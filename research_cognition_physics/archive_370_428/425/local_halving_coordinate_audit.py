"""Round 425: consistent local halving, NSS, and a nonabelian chart witness.

The general implications are proved in research_note_425.md. Finite tests
check certificates and witnesses, not arbitrary topological groups.
"""
import argparse
from fractions import Fraction as F
import io
import json
from pathlib import Path
import platform
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'local_halving_coordinate_audit_results.json'
OBS = {}
I = np.eye(2, dtype=complex)
SIGMA = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
RADIUS = 0.75  # strictly smaller than pi/2
ETA = 0.8


def su2_exp(x):
    x = np.asarray(x, dtype=float)
    radius = np.linalg.norm(x)
    if radius == 0:
        return I.copy()
    return np.cos(radius) * I + 1j * np.sin(radius) * np.einsum('i,ijk->jk', x / radius, SIGMA)


def su2_log(g):
    """Principal vector logarithm on SU(2) away from -I."""
    scalar = float(np.trace(g).real / 2)
    vector = np.array([(np.trace(s @ g) / (2j)).real for s in SIGMA])
    length = np.linalg.norm(vector)
    if length < 1e-14:
        if scalar < 0:
            raise ValueError('The logarithm is not unique at -I')
        return vector
    return np.arctan2(length, scalar) * vector / length


def half(g):
    return su2_exp(su2_log(g) / 2)


def cost(g):
    return float(np.linalg.norm(su2_log(g)))


def effect(g):
    x = su2_log(g)
    r = np.linalg.norm(x)
    if r == 0:
        raise ValueError('Direction is undefined at the identity')
    return (I + ETA * np.einsum('i,ijk->jk', x / r, SIGMA)) / 2


def direction_probabilities(g):
    e = effect(g)
    return np.array([np.trace((I + s) @ e / 2).real for s in SIGMA])


def reconstruct(g):
    """Witness only: assumes exact cost and three calibrated probabilities."""
    return cost(g) * (2 * direction_probabilities(g) - 1) / ETA


def escape_depth(b, bound, q):
    assert 0 < b <= bound and 0 < q < 1
    depth = 0
    while q ** depth * bound >= b:
        depth += 1
    return depth


def torus_cost(angles):
    total = F(0)
    for j, angle in enumerate(angles, 1):
        t = angle % 1
        total += F(1, 2 ** j) * min(t, 1 - t)
    return total


class Audit(unittest.TestCase):
    def test_01_exact_finite_escape_certificates(self):
        certificates = []
        for q in (F(1, 2), F(2, 3), F(3, 4)):
            for b in (F(1, 3), F(1, 8), F(7, 64)):
                k = escape_depth(b, F(1), q)
                self.assertLess(q ** k, b)
                self.assertGreaterEqual(q ** (k - 1), b)
                certificates.append(dict(q=str(q), initial_cost=str(b), depth=k,
                                         maximal_replay_count=2 ** k))
        # Explicit additive witness with a strict local unit budget.
        b = F(7, 64)
        first_exit = next(k for k in range(10) if 2 ** k * b >= 1)
        self.assertLessEqual(first_exit, escape_depth(b, F(1), F(1, 2)))
        OBS['exact_escape_certificates'] = certificates
        OBS['additive_witness_first_exit_depth'] = first_exit

    def test_02_su2_local_root_and_replay_consistency(self):
        rng = np.random.default_rng(42502)
        square_error = replay_error = budget_error = identity_error = 0.0
        for _ in range(120):
            direction = rng.normal(size=3)
            direction /= np.linalg.norm(direction)
            x = rng.uniform(0.001, 0.99 * RADIUS) * direction
            g = su2_exp(x)
            h = half(g)
            square_error = max(square_error, np.linalg.norm(h @ h - g))
            budget_error = max(budget_error, abs(cost(h) - cost(g) / 2))
            identity_error = max(identity_error, np.linalg.norm(g.conj().T @ g - I))
            # This smaller input ensures that g and g^2 both remain in U.
            a = su2_exp(x / 2)
            replay_error = max(replay_error, np.linalg.norm(half(a @ a) - a))
        self.assertLess(max(square_error, replay_error, budget_error, identity_error), 3e-14)
        OBS['local_square_root_error'] = float(square_error)
        OBS['repeat_then_half_error'] = float(replay_error)
        OBS['local_halving_budget_error'] = float(budget_error)
        OBS['su2_unitarity_error'] = float(identity_error)

    def test_03_halving_is_not_a_group_homomorphism(self):
        g = su2_exp([0.24, 0, 0])
        h = su2_exp([0, 0.28, 0])
        self.assertLess(max(cost(g), cost(h), cost(g @ h)), RADIUS)
        defect = np.linalg.norm(half(g @ h) - half(g) @ half(h))
        commutator = np.linalg.norm(g @ h - h @ g)
        self.assertGreater(defect, 0.02)
        self.assertGreater(commutator, 0.15)
        self.assertTrue(np.allclose((-I) @ (-I), I))
        OBS['local_root_product_defect'] = float(defect)
        OBS['endpoint_noncommutativity'] = float(commutator)
        OBS['central_order_two_element_excluded_from_local_domain'] = True

    def test_04_complete_local_qubit_direction_witness(self):
        rng = np.random.default_rng(42504)
        covariance_error = inverse_error = budget_error = 0.0
        for _ in range(80):
            v = rng.normal(size=3)
            v *= 0.35 / np.linalg.norm(v)
            g = su2_exp(v)
            u = su2_exp(rng.normal(size=3))
            transformed = u @ g @ u.conj().T
            covariance_error = max(covariance_error,
                                   np.linalg.norm(effect(transformed) - u @ effect(g) @ u.conj().T))
            budget_error = max(budget_error, abs(cost(transformed) - cost(g)))
            inverse_error = max(inverse_error, np.linalg.norm(effect(g) + effect(g.conj().T) - I))
            self.assertGreaterEqual(np.linalg.eigvalsh(effect(g))[0], 0.099999999)
        differences = []
        for axis in np.eye(3):
            g = su2_exp(0.3 * axis)
            difference = effect(g) - effect(g.conj().T)
            differences.append([np.trace(s @ difference).real / 2 for s in SIGMA])
        rank = np.linalg.matrix_rank(differences)
        self.assertEqual(rank, 3)
        self.assertLess(max(covariance_error, inverse_error, budget_error), 3e-14)
        OBS['qubit_covariance_error'] = float(covariance_error)
        OBS['inverse_effect_complement_error'] = float(inverse_error)
        OBS['conjugation_budget_error'] = float(budget_error)
        OBS['effect_difference_rank'] = int(rank)
        OBS['inverse_operator_gap'] = float(np.linalg.norm(effect(su2_exp([0.3, 0, 0])) -
                                                           effect(su2_exp([-0.3, 0, 0])), 2))

    def test_05_local_chart_and_nonadditive_overlap(self):
        rng = np.random.default_rng(42505)
        reconstruction_error = overlap_error = 0.0
        a, b = su2_exp([0.09, -0.04, 0.02]), su2_exp([-0.02, 0.07, 0.03])
        for _ in range(80):
            x = rng.normal(size=3)
            x *= rng.uniform(0.01, 0.3) / np.linalg.norm(x)
            endpoint = a @ su2_exp(x)
            relative_b = b.conj().T @ endpoint
            self.assertLess(cost(relative_b), RADIUS)
            decoded = reconstruct(relative_b)
            reconstruction_error = max(reconstruction_error, np.linalg.norm(decoded - su2_log(relative_b)))
            overlap_error = max(overlap_error, np.linalg.norm(b @ su2_exp(decoded) - endpoint))
        x, y = np.array([0.24, 0, 0]), np.array([0, 0.28, 0])
        additive_defect = np.linalg.norm(su2_log(su2_exp(x) @ su2_exp(y)) - x - y)
        self.assertLess(max(reconstruction_error, overlap_error), 3e-14)
        self.assertGreater(additive_defect, 0.06)
        OBS['cost_plus_probabilities_coordinate_error'] = float(reconstruction_error)
        OBS['nonabelian_chart_overlap_error'] = float(overlap_error)
        OBS['additive_coordinate_rule_defect'] = float(additive_defect)
        OBS['readout_requires_exact_cost_and_three_probabilities'] = True

    def test_06_compact_infinite_torus_contraction_is_not_halving(self):
        # Finite-support identities represent exact points in the infinite product.
        angles = (F(1, 3), F(3, 5), F(1, 7))
        self.assertEqual(torus_cost((F(0),) + angles), torus_cost(angles) / 2)
        certificates = []
        for n in (1, 4, 10, 20):
            torsion = (F(0),) * n + (F(1, 2),)
            squared = tuple(2 * x % 1 for x in torsion)
            self.assertEqual(torus_cost(squared), 0)
            self.assertEqual(torus_cost(torsion), F(1, 2 ** (n + 2)))
            certificates.append(dict(fixed_initial_coordinates=n,
                                     whole_tail_subgroup_budget_upper_bound=str(F(1, 2 ** (n + 1))),
                                     nonidentity_order_two_cost=str(torus_cost(torsion))))
        OBS['infinite_torus_shift_exact_ratio'] = '1/2'
        OBS['infinite_torus_small_subgroup_certificates'] = certificates

    def test_07_halving_does_not_select_dimension_or_supply_compactness(self):
        errors = []
        for n in (1, 2, 3, 4, 8, 32):
            vector = np.arange(1, n + 1, dtype=float)
            vector *= 0.3 / np.linalg.norm(vector)
            errors.append(abs(np.linalg.norm(vector / 2) - np.linalg.norm(vector) / 2))
        self.assertEqual(max(errors), 0)
        # Truncated orthogonal family checks the formula used in the l2 proof.
        points = 0.2 * np.eye(32)
        distances = [np.linalg.norm(points[i] - points[j]) for i in range(32) for j in range(i)]
        self.assertTrue(np.allclose(distances, 0.2 * np.sqrt(2)))
        OBS['finite_additive_dimensions_with_same_halving_contract'] = [1, 2, 3, 4, 8, 32]
        OBS['hilbert_orthogonal_family_separation'] = float(min(distances))
        OBS['infinite_dimension_and_noncompactness_proved_analytically'] = True

    def test_08_additive_error_leaves_a_nonzero_resolution_floor(self):
        q, epsilon, bound = F(1, 2), F(1, 100), F(1)
        floor = epsilon / (1 - q)
        current = bound
        for k in range(1, 15):
            current = q * current + epsilon
            exact = q ** k * bound + epsilon * (1 - q ** k) / (1 - q)
            self.assertEqual(current, exact)
            self.assertGreater(current, floor)
        self.assertEqual(floor, q * floor + epsilon)
        OBS['approximate_budget_recurrence_floor'] = str(floor)
        OBS['finite_error_claims_exact_NSS'] = False


def run():
    output = io.StringIO()
    tests = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=425, baseline_round=424, status='verified_conditional_local_chart_bridge',
                tests_run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors),
                python=platform.python_version(), numpy=np.__version__, observations=OBS,
                scope=dict(local_halving_implies_NSS=True,
                           locally_compact_NSS_Lie_theorem_reused=True,
                           global_contraction_automorphism_required=False,
                           halving_preserves_group_products_required=False,
                           actual_endpoint_torsor_still_input=True,
                           local_compactness_still_input=True,
                           compatible_local_halving_and_cost_still_input=True,
                           complete_local_qubit_shell_still_separate_input=True,
                           conditional_three_dimensional_local_atlas=True,
                           additive_displacement_coordinates_inherited=False,
                           full_round424_contract_claimed_logically_weakened=False,
                           budget_numerical_access_in_witness=True,
                           witness_claimed_complete_cognitive_realization=False,
                           three_dimensional_space_unconditionally_derived=False,
                           full_cognitive_countermodel_completed=False,
                           full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        saved = json.loads(TARGET.read_text(encoding='utf-8'))
        assert saved == result
    else:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
