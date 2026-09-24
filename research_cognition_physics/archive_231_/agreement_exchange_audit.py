"""Round 429: agreement-preserving reversible pair rules.

Schur--Weyl duality and partial swap are established results, explicitly cited
in the note. The tests audit the premise-to-rule bridge and its two distinct
agreement domains. No space, autonomous scheduler or universal processor is
constructed here.
"""
import argparse
import io
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

TARGET = Path(__file__).with_name('agreement_exchange_audit_results.json')
OBS = {}


def swap(d):
    out = np.zeros((d*d, d*d), dtype=complex)
    for i in range(d):
        for j in range(d):
            out[j*d+i, i*d+j] = 1
    return out


def unitary(h, time):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*time*values)) @ vectors.conj().T


def state(d, rng):
    a = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
    out = a @ a.conj().T
    return out/np.trace(out)


def marginals(joint, d):
    tensor = joint.reshape(d, d, d, d)
    return np.trace(tensor, axis1=1, axis2=3), np.trace(tensor, axis1=0, axis2=2)


def distance(a, b):
    difference = a-b
    return float(np.sum(np.abs(np.linalg.eigvalsh((difference+difference.conj().T)/2)))/2)


def su_generators(d):
    out = []
    for j in range(d):
        for k in range(j+1, d):
            x = np.zeros((d, d), dtype=complex)
            x[j, k] = x[k, j] = 1
            y = np.zeros_like(x)
            y[j, k], y[k, j] = -1j, 1j
            out.extend([x, y])
    for j in range(d-1):
        z = np.zeros((d, d), dtype=complex)
        z[j, j], z[j+1, j+1] = 1, -1
        out.append(z)
    return out


def modular_rank(matrix, modulus=1009):
    # Gaussian integers mapped to F_1009 with i -> 469; 469^2 = -1.
    assert 469**2 % modulus == modulus-1
    assert np.all(matrix.real == np.rint(matrix.real))
    assert np.all(matrix.imag == np.rint(matrix.imag))
    a = (np.rint(matrix.real).astype(np.int64)+469*np.rint(matrix.imag).astype(np.int64)) % modulus
    rank = 0
    pivots = []
    for col in range(a.shape[1]):
        candidates = np.flatnonzero(a[rank:, col])
        if not len(candidates):
            continue
        row = rank+int(candidates[0])
        a[[rank, row]] = a[[row, rank]]
        pivot = int(a[rank, col])
        pivots.append(pivot)
        a[rank] = (a[rank]*pow(pivot, -1, modulus)) % modulus
        if rank+1 < a.shape[0]:
            a[rank+1:] = (a[rank+1:]-a[rank+1:, col, None]*a[rank]) % modulus
        rank += 1
        if rank == a.shape[0]:
            break
    return rank, pivots


def commutator_constraints(generators):
    n = generators[0].shape[0]
    identity = np.eye(n)
    # Row-major vectorization of Q H - H Q.
    return np.vstack([np.kron(q, identity)-np.kron(identity, q.T) for q in generators])


class Audit(unittest.TestCase):
    def test_01_finite_exact_commutant_certificates(self):
        rows = []
        for d in (2, 3):
            identity = np.eye(d)
            local = su_generators(d)
            common = [np.kron(x, identity)+np.kron(identity, x) for x in local]
            independent = [np.kron(x, identity) for x in local]+[np.kron(identity, x) for x in local]
            c = commutator_constraints(common)
            separate = commutator_constraints(independent)
            common_rank, pivots = modular_rank(c)
            separate_rank, other_pivots = modular_rank(separate)
            self.assertEqual(common_rank, d**4-2)
            self.assertEqual(separate_rank, d**4-1)
            self.assertTrue(np.array_equal(c @ np.eye(d*d).reshape(-1), np.zeros(c.shape[0])))
            self.assertTrue(np.array_equal(c @ swap(d).reshape(-1), np.zeros(c.shape[0])))
            rows.append(dict(dimension=d, matrix_unknowns=d**4,
                             common_rank_mod_1009=common_rank,
                             independent_rank_mod_1009=separate_rank,
                             common_nonzero_pivots=pivots, independent_nonzero_pivots=other_pivots))
        OBS['finite_commutant_certificates'] = rows

    def test_02_equal_product_agreement_and_common_covariance(self):
        rng = np.random.default_rng(42902)
        product_error = covariance_error = group_error = 0.
        for d in (2, 3):
            s = swap(d)
            h = .41*np.eye(d*d)+.73*s
            w = unitary(h, .62)
            group_error = max(group_error, float(np.linalg.norm(unitary(h, .21) @ unitary(h, .41)-w)))
            for _ in range(5):
                rho = state(d, rng)
                before = np.kron(rho, rho)
                after = w @ before @ w.conj().T
                product_error = max(product_error, float(np.linalg.norm(after-before)))
                a, b = marginals(after, d)
                self.assertLess(distance(a, rho), 1e-12)
                self.assertLess(distance(b, rho), 1e-12)
                q = state(d, rng)
                v = unitary(q, 2.3)
                lift = np.kron(v, v)
                covariance_error = max(covariance_error, float(np.linalg.norm(w @ lift-lift @ w)))
        self.assertLess(max(product_error, covariance_error, group_error), 1e-12)
        OBS['agreement_rule'] = dict(equal_product_error=product_error,
            common_basis_covariance_error=covariance_error, homogeneous_group_error=group_error,
            all_dimensions_proof_is_analytic=True)

    def test_03_internal_state_response_and_unknown_reference(self):
        rng = np.random.default_rng(42903)
        worst = 0.
        for d in (2, 3):
            reference = 3
            theta = .61
            c, s = math.cos(theta), math.sin(theta)
            w = unitary(swap(d), theta)
            rho = state(d*reference, rng)
            sigma = state(d, rng)
            # Ordering D,P,R. Sigma is prepared independently of the entire D,R.
            initial = np.einsum('arbs,pq->aprbqs', rho.reshape(d, reference, d, reference), sigma)
            initial = initial.reshape(d*d*reference, d*d*reference)
            lifted = np.kron(w, np.eye(reference))
            final = lifted @ initial @ lifted.conj().T
            tensor = final.reshape(d, d, reference, d, d, reference)
            actual = np.trace(tensor, axis1=1, axis2=4).reshape(d*reference, d*reference)
            rho_r = np.trace(rho.reshape(d, reference, d, reference), axis1=0, axis2=2)
            sigma_lift = np.kron(sigma, np.eye(reference))
            expected = c*c*rho+s*s*np.kron(sigma, rho_r)-1j*c*s*(sigma_lift @ rho-rho @ sigma_lift)
            worst = max(worst, float(np.linalg.norm(actual-expected)))
            self.assertLess(np.linalg.norm(actual-expected), 1e-12)
            self.assertLess(np.linalg.norm(lifted.conj().T @ final @ lifted-initial), 1e-12)
        # One fixed H and time, differing internal target preparations.
        data = np.diag([1., 0.])
        targets = [np.diag([1., 0.]), np.diag([0., 1.])]
        w = unitary(swap(2), math.pi/4)
        effect = np.diag([0., 1.])
        probabilities = []
        for target in targets:
            out, _ = marginals(w @ np.kron(data, target) @ w.conj().T, 2)
            probabilities.append(float(np.trace(effect @ out).real))
        self.assertAlmostEqual(probabilities[1]-probabilities[0], .5)
        OBS['internal_target_response'] = dict(reference_formula_error=worst,
            same_rule_same_time_readout_probabilities=probabilities,
            full_joint_reference_recoverable_by_inverse=True, source_preparation_is_explicit_input=True)

    def test_04_pure_agreement_does_not_control_mixed_agreement(self):
        # All coefficients of U and rho are binary-exact rationals.
        d = 3
        a = np.zeros(d*d)
        a[1], a[3] = 1., -1.
        h = np.outer(a, a)/2
        w = np.eye(d*d)-2*h
        rng = np.random.default_rng(42904)
        error = 0.
        for _ in range(12):
            psi = rng.normal(size=d)+1j*rng.normal(size=d)
            psi /= np.linalg.norm(psi)
            identical = np.kron(psi, psi)
            error = max(error, float(np.linalg.norm(w @ identical-identical)))
        self.assertLess(error, 1e-13)
        rho = np.array([[.5, 0., 0.], [0., .25, .25], [0., .25, .25]])
        target = np.array([[.5, 0., 0.], [0., .25, .125], [0., .125, .25]])
        left, right = marginals(w @ np.kron(rho, rho) @ w.T, d)
        self.assertTrue(np.array_equal(left, target))
        self.assertTrue(np.array_equal(right, target))
        effect = np.array([[0., 0., 0.], [0., .5, .5], [0., .5, .5]])
        gap = float(np.trace(effect @ (rho-left)))
        self.assertEqual(gap, .125)
        self.assertAlmostEqual(distance(rho, left), .125)
        OBS['pure_only_boundary'] = dict(dimension=3, sampled_pure_agreement_error=error,
            mixed_source=rho.tolist(), both_output_marginals=target.tolist(),
            exact_readout_probability_gap='1/8', pure_state_all_quantifier_proved_by_antisymmetry=True)

    def test_05_equal_marginals_with_hidden_correlation(self):
        psi = np.array([0., 1., 1j, 0.])/math.sqrt(2)
        before = np.outer(psi, psi.conj())
        a, b = marginals(before, 2)
        self.assertTrue(np.allclose(a, np.eye(2)/2))
        self.assertTrue(np.allclose(a, b))
        w = unitary(swap(2), math.pi/4)
        after = w @ before @ w.conj().T
        target = np.diag([0., 1., 0., 0.])
        self.assertLess(np.linalg.norm(after-target), 1e-12)
        left, right = marginals(after, 2)
        self.assertAlmostEqual(distance(left, a), .5)
        self.assertAlmostEqual(distance(right, b), .5)
        self.assertLess(np.linalg.norm(w.conj().T @ after @ w-before), 1e-12)
        OBS['correlated_agreement_boundary'] = dict(initial_marginals='I/2, I/2',
            final_joint_state='|01><01|', each_marginal_trace_distance=.5,
            full_joint_information_erased=False,
            all_correlated_equal_marginals_would_force_trivial_generator=True)

    def test_06_finite_interaction_is_not_an_exact_programmed_unitary(self):
        rows = []
        for d in (2, 3):
            sigma = np.diag(np.arange(1, d+1, dtype=float))
            sigma /= np.trace(sigma)
            first, second = np.zeros((d, d)), np.zeros((d, d))
            first[0, 0], second[1, 1] = 1, 1
            for theta in (.2, .7, math.pi/2):
                w = unitary(swap(d), theta)
                a, _ = marginals(w @ np.kron(first, sigma) @ w.conj().T, d)
                b, _ = marginals(w @ np.kron(second, sigma) @ w.conj().T, d)
                observed = distance(a, b)
                self.assertAlmostEqual(observed, math.cos(theta)**2)
                rows.append(dict(dimension=d, angle=theta, input_distance=1., output_distance=observed))
        OBS['finite_control_boundary'] = dict(examples=rows,
            nontrivial_single_interaction_gives_unitary_data_channel=False,
            infinitesimal_commutator_not_promoted_to_finite_exact_control=True,
            no_free_repeated_copy_assumption=True)


def run():
    output = io.StringIO()
    tests = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=429, baseline_round=428, status='agreement_selected_pair_rule_with_explicit_domain',
        tests_run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(equal_independent_states_preservation_is_new_candidate_input=True,
            same_dimension_pair_and_identification_required=True,
            reversible_pair_rule_classified=True,
            fixed_exchange_rule_responds_to_internal_state=True,
            mixed_agreement_not_inferred_from_pure_agreement=True,
            arbitrary_correlated_agreement_would_force_triviality=True,
            swap_or_schur_weyl_claimed_as_new_discovery=False,
            nonzero_coupling_or_pair_selection_derived=False,
            exact_universal_autonomous_control_completed=False,
            full_cognitive_implementation_completed=False,
            spatial_dimension_selected=False,
            three_dimensional_space_unconditionally_derived=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert json.loads(TARGET.read_text(encoding='utf-8')) == result
    else:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
