"""Round 432: which exchange weights can change internal relations?

Central class-sum facts are standard. The audit connects their precise scope to
429--431 and checks an expectation-feedback shortcut against affine evolution.
"""
import argparse
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import agreement_exchange_audit as exact
import exchange_relation_audit as relation
import exchange_reader_audit as reader

TARGET = Path(__file__).with_name('exchange_participation_audit_results.json')
OBS = {}


def swap(n, d, first, second):
    out = np.zeros((d**n, d**n), complex)
    for col in range(d**n):
        digits = []
        value = col
        for _ in range(n):
            digits.append(value % d)
            value //= d
        digits.reverse()
        digits[first], digits[second] = digits[second], digits[first]
        row = 0
        for value in digits:
            row = row*d+value
        out[row, col] = 1
    return out


def family(n, d=2):
    return {edge: swap(n, d, *edge) for edge in itertools.combinations(range(n), 2)}


def power(rho, n):
    out = np.array([[1.]], complex)
    for _ in range(n):
        out = np.kron(out, rho)
    return out


class Audit(unittest.TestCase):
    def close(self, a, b, tolerance=5e-12):
        self.assertLess(float(np.linalg.norm(a-b)), tolerance)

    def test_01_central_sum_and_exact_kernel_ranks(self):
        certificates = []
        for n, d in ((3, 2), (4, 2), (5, 2), (3, 3)):
            fs = family(n, d)
            generators = list(fs.values())
            central = sum(generators)
            for s in generators:
                self.close(central@s, s@central)
            columns = []
            for h in generators:
                columns.append(np.concatenate([(h@s-s@h).ravel() for s in generators]))
            matrix = np.column_stack(columns)
            rank, pivots = exact.modular_rank(matrix)
            self.assertEqual(rank, len(generators)-1)
            self.close(matrix@np.ones(len(generators)), np.zeros(matrix.shape[0]))
            certificates.append(dict(subject_count=n, local_dimension=d, weights=len(generators),
                rank_mod_1009=rank, kernel_dimension=1, nonzero_pivots=pivots))
        OBS['uniformity_certificates'] = certificates

    def test_02_one_excitation_reconstructs_weights(self):
        rows = []
        for n, d in ((4, 2), (3, 3)):
            fs = family(n, d)
            weights = {edge: j-2 for j, edge in enumerate(fs)}
            h = sum(weights[edge]*s for edge, s in fs.items())
            indices = [d**(n-1-k) for k in range(n)]
            reduced = h[np.ix_(indices, indices)]
            for edge, value in weights.items():
                self.assertEqual(reduced[edge], value)
            self.close(reduced@np.ones(n), sum(weights.values())*np.ones(n))
            rows.append(dict(subject_count=n, local_dimension=d,
                weights=[int(x) for x in weights.values()],
                off_diagonal_weight_reconstruction_exact=True))
        OBS['necessity_witness'] = rows

    def test_03_uniform_completion_freezes_reader_setup(self):
        _, effect, plus, minus, _ = reader.setup()
        central = sum(family(5).values())
        for rho in (plus, minus):
            self.close(central@rho, rho@central)
        errors = []
        for t in (.17, 1.5, 2.3):
            u = relation.evolve(central, t)
            for rho in (plus, minus):
                after = u@rho@u.conj().T
                errors.append(float(np.linalg.norm(after-rho)))
                self.close(after, rho)
                self.assertAlmostEqual(np.trace(effect@after).real, 1)
        OBS['reader_uniform_completion'] = dict(all_ten_pairs_have_same_strength=True,
            initial_states_unchanged_from_431=True, maximum_state_change=max(errors),
            reader_probability_for_both_inputs=1., reader_contrast=0.)

    def test_04_pair_current_from_incident_weight_differences(self):
        n = 4
        fs = family(n)
        weights = {edge: j-2 for j, edge in enumerate(fs)}
        h = sum(weights[edge]*s for edge, s in fs.items())
        def weight(i, j):
            return weights[tuple(sorted((i, j)))]
        errors = []
        for (a, b), sab in fs.items():
            expected = np.zeros_like(h)
            for k in range(n):
                if k not in (a, b):
                    sak = fs[tuple(sorted((a, k)))]
                    expected += (weight(a, k)-weight(b, k))*1j*(sak@sab-sab@sak)
            actual = 1j*(h@sab-sab@h)
            errors.append(float(np.linalg.norm(actual-expected)))
            self.close(actual, expected)
        OBS['pair_relation_current'] = dict(subject_count=n,
            maximum_operator_identity_error=max(errors),
            current_depends_on_weight_differences_and_three_body_chirality=True)

    def test_05_equal_product_condition_does_not_select_weights(self):
        rng = np.random.default_rng(43205)
        fs = family(4)
        rho = relation.density(2, rng)
        product = power(rho, 4)
        examples = []
        for weights in ([1]*6, [1, 0, 0, 1, 0, 1], [0, 2, -1, 3, 1, -2]):
            h = sum(value*s for value, s in zip(weights, fs.values()))
            self.close(h@product, product@h)
            u = relation.evolve(h, .37)
            self.close(u@product@u.conj().T, product)
            examples.append(dict(weights=weights, same_independent_product_preserved=True))
        OBS['same_state_condition'] = dict(examples=examples,
            arbitrary_weights_preserve_rho_tensor_n_analytically=True,
            coupling_pattern_selected_by_agreement_alone=False)

    def test_06_expectation_weight_feedback_is_not_affine(self):
        a, b, c, p, x, y, z = relation.operators()
        rho_x, rho_y = (p+x)/4, (p+y)/4
        mixed = (rho_x+rho_y)/2
        def proposed_derivative(rho):
            h = sum(float(np.trace(s@rho).real)*s for s in (a, b, c))
            return -1j*(h@rho-rho@h)
        self.close(proposed_derivative(rho_x), np.zeros((8, 8)))
        self.close(proposed_derivative(rho_y), np.zeros((8, 8)))
        self.close(proposed_derivative(mixed), 3*z/16)
        # Analytic solution for mixed initial r=(1/2,1/2,0):
        # r_x=(1/2)cos(3t/2), r_y=1/2, r_z=(1/2)sin(3t/2).
        residuals = []
        for t in (.17, .64, math.pi/3):
            r = np.array([.5*math.cos(1.5*t), .5, .5*math.sin(1.5*t)])
            dr = np.array([-.75*math.sin(1.5*t), 0., .75*math.cos(1.5*t)])
            rho = (p+r[0]*x+r[1]*y+r[2]*z)/4
            predicted = (dr[0]*x+dr[1]*y+dr[2]*z)/4
            residuals.append(float(np.linalg.norm(proposed_derivative(rho)-predicted)))
            self.close(proposed_derivative(rho), predicted)
        direct = (p+.5*y+.5*z)/4
        branch_average = mixed
        effect = (p+z)/2
        gap = float(np.trace(effect@(direct-branch_average)).real)
        self.assertAlmostEqual(gap, .25)
        OBS['expectation_feedback_boundary'] = dict(rule='J_ij(rho)=Tr(rho SWAP_ij)',
            both_pure_logical_components_stationary=True,
            mixture_changes=True, time='pi/3', effect_probability_affinity_defect=gap,
            maximum_analytic_ODE_residual=max(residuals),
            proposed_rule_is_fixed_CPTP_evolution=False,
            mean_field_or_explicit_controller_alternatives_refuted=False)


def run():
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=432, baseline_round=431, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(central_exchange_sum_is_established_algebra=True,
            uniform_weights_iff_all_relations_constant_under_model=True,
            arbitrary_initial_state_cannot_change_invariant_reads_under_uniform_rule=True,
            equal_product_condition_does_not_select_couplings=True,
            expectation_weight_shortcut_affinity_defect_proved=True,
            conditional_pair_rule_429_refuted=False,
            arbitrary_continuous_internal_evolution_refuted=False,
            full_cognitive_implementation_completed=False,
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
