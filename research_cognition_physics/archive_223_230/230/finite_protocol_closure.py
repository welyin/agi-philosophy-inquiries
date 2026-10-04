"""Round 230: finite adaptive experiments and the scope of operational density.

Checks telescoping error bounds with retained outcome histories and entangled
references, the optimal discrimination identity, and rare-event conditioning.
The general theorem is analytic. Sample perturbations validate formulas; they
are not an implementation of an unspecified actual dense gate set.
"""
import argparse
import json
from pathlib import Path
import platform
import unittest

import numpy as np

from reversible_dynamics_bridge import dagger, sample_generator, unitary
from tensor_process_bridge import maximally_entangled


def trace_norm(x):
    return float(np.sum(np.linalg.svd(x, compute_uv=False)))


def adaptive_step(stage, history):
    """Two efficient branches, with different instruments on different histories."""
    code = sum((value+1)*3**i for i, value in enumerate(history))
    k = sample_generator(2, 23000+100*stage+code)
    u = unitary(k, .13+.04*stage)
    gamma = .15+.06*((stage+code) % 7)
    damping = [np.diag([1., np.sqrt(1-gamma)]), np.array([[0, np.sqrt(gamma)], [0, 0]])]
    return [operator@u for operator in damping]


def perturbing_unitary(stage, history, epsilon):
    code = sum((value+1)*3**i for i, value in enumerate(history))
    k = sample_generator(2, 23100+100*stage+code)
    k /= np.linalg.norm(k, 2)
    return unitary(k, epsilon)


def run_protocol(depth, epsilon):
    """Block dictionaries retain the full classical outcome history and reference."""
    initial = maximally_entangled(2)
    exact = {(): initial}
    nearby = {(): initial}
    budgets = []
    for stage in range(depth):
        next_exact, next_nearby = {}, {}
        step_budget = 0.
        for history in exact:
            instrument = adaptive_step(stage, history)
            v = perturbing_unitary(stage, history, epsilon)
            step_budget = max(step_budget, 2*float(np.linalg.norm(v-np.eye(2), 2)))
            for result, operator in enumerate(instrument):
                key = history+(result,)
                total = np.kron(operator, np.eye(2))
                perturbed = np.kron(v@operator, np.eye(2))
                next_exact[key] = total@exact[history]@dagger(total)
                next_nearby[key] = perturbed@nearby[history]@dagger(perturbed)
        exact, nearby = next_exact, next_nearby
        budgets.append(step_budget)
    return exact, nearby, budgets


def experiment_certificate(depth, epsilon):
    exact, nearby, budgets = run_protocol(depth, epsilon)
    flagged_norm = sum(trace_norm(exact[h]-nearby[h]) for h in exact)
    # A concrete final Z measurement on S, retaining its result and the history.
    record_tv = 0.
    for h in exact:
        difference = (exact[h]-nearby[h]).reshape(2, 2, 2, 2)
        for value in range(2):
            record_tv += abs(float(np.trace(difference[value, :, value, :]).real)) / 2
    ideal_trace = sum(float(np.trace(block).real) for block in exact.values())
    nearby_trace = sum(float(np.trace(block).real) for block in nearby.values())
    return {'depth': depth, 'histories': len(exact), 'perturbation': epsilon,
            'per_stage_diamond_upper_bounds': budgets,
            'analytic_flagged_norm_upper_bound': sum(budgets),
            'observed_flagged_trace_norm_on_bell_input': flagged_norm,
            'observed_final_record_total_variation': record_tv,
            'analytic_record_total_variation_upper_bound': min(1., sum(budgets)/2),
            'normalization_error': max(abs(ideal_trace-1), abs(nearby_trace-1))}


def rare_event(p):
    zero = np.diag([1., 0.])
    one = np.diag([0., 1.])
    exact = [p*zero, (1-p)*zero]
    nearby = [p*one, (1-p)*zero]
    unconditional = sum(trace_norm(x-y) for x, y in zip(exact, nearby))
    return {'success_probability': p,
            'exact_replacement_instrument_diamond_distance': unconditional,
            'conditional_trace_distance': trace_norm(exact[0]/p-nearby[0]/p)/2,
            'expected_independent_trials_until_success': 1/p}


def quantum_discrimination_certificate(epsilon):
    zero = np.diag([1., 0.]).astype(complex)
    y = np.array([[0, -1j], [1j, 0]])
    v = unitary(y, epsilon)
    alternative = v@zero@dagger(v)
    difference = zero-alternative
    values, vectors = np.linalg.eigh(difference)
    projector = vectors[:, values > 1e-12]@dagger(vectors[:, values > 1e-12])
    success = .5*(np.trace(projector@zero)+np.trace((np.eye(2)-projector)@alternative)).real
    return {'observed_optimal_success': float(success),
            'helstrom_value': .5+trace_norm(difference)/4,
            'trace_distance': trace_norm(difference)/2}


class FiniteProtocolTests(unittest.TestCase):
    def test_adaptive_instruments_preserve_total_probability(self):
        for stage in range(4):
            for history in ((), (0,), (1, 0), (1, 1, 1)):
                operators = adaptive_step(stage, history)
                np.testing.assert_allclose(sum(dagger(k)@k for k in operators), np.eye(2), atol=2e-14)

    def test_flagged_entangled_outputs_obey_telescoping_bound(self):
        for depth in (1, 2, 4):
            for epsilon in (.04, .004):
                c = experiment_certificate(depth, epsilon)
                self.assertLessEqual(c['observed_flagged_trace_norm_on_bell_input'], c['analytic_flagged_norm_upper_bound']+1e-12)
                self.assertLess(c['normalization_error'], 2e-13)

    def test_terminal_measurement_contracts_to_total_variation(self):
        for depth in (1, 2, 4):
            c = experiment_certificate(depth, .04)
            self.assertLessEqual(c['observed_final_record_total_variation'], c['observed_flagged_trace_norm_on_bell_input']/2+1e-12)

    def test_exact_match_has_no_artificial_error(self):
        c = experiment_certificate(4, 0.)
        self.assertLess(c['observed_flagged_trace_norm_on_bell_input'], 1e-13)

    def test_postselection_requires_a_probability_denominator(self):
        for p in (.1, .001, .00001):
            c = rare_event(p)
            self.assertAlmostEqual(c['exact_replacement_instrument_diamond_distance'], 2*p)
            self.assertAlmostEqual(c['conditional_trace_distance'], 1.)

    def test_positive_part_measurement_attains_discrimination_formula(self):
        for epsilon in (.02, .2, .7):
            c = quantum_discrimination_certificate(epsilon)
            self.assertAlmostEqual(c['observed_optimal_success'], c['helstrom_value'])

    def test_repeated_rare_event_records_can_accumulate_evidence(self):
        p, trials = .01, 50
        # Correct distinguishability needs the number of attempts in the budget.
        cumulative = 1-(1-p)**trials
        self.assertGreater(cumulative, p)
        self.assertLessEqual(cumulative, trials*p)


def report():
    return {'round': 230, 'date': '2026-09-21',
            'analytic_result': 'For a finite adaptive protocol, sum of instrument diamond errors E bounds final flagged trace norm by E and record total variation by min(1,E/2). No finite fixed protocol has a positive robust separation from the dense baseline theory.',
            'quantifiers': 'For every fixed finite quantum protocol and positive tolerance there exists an allowed approximating protocol. Not one universal approximation for all protocols, not exact equality, not equality of efficient resource costs.',
            'hypotheses': ['F+U+C+P', 'round228 density with classical records and external references', 'finite maximum protocol depth and finite outcomes', 'comparison with the same existing input/output system types'],
            'adaptive_protocol_certificates': [experiment_certificate(d,e) for d in (1,2,4) for e in (.04,.004)],
            'rare_event_certificates': [rare_event(p) for p in (.1,.001,.00001)],
            'discrimination_certificates': [quantum_discrimination_certificate(e) for e in (.02,.2,.7)],
            'numerical_scope': 'Finite examples with an explicit Bell reference; the reported trace norms are not numerically computed diamond norms. The diamond bounds and all-reference statements are analytic.',
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
            'unresolved_but_not_blocking_closure': ['Does bare F+U+C+P imply a nonconstant reversible operation path?', 'Exact gate synthesis cost and natural Hamiltonian selection are not fixed.'],
            'next': 'Stage closure complete; await the user direction conjectures before beginning Hamiltonians, fields or gravity.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FiniteProtocolTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': tests.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        target = Path(__file__).with_name('finite_protocol_closure_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing result differs; review before replacing it.')
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
