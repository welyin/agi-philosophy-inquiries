"""Rounds 201/202: a standard real-qubit comparison for delayed self-prediction.

The quantum memory and Born readout are explicit comparison assumptions, not
derived from cognitive principles. This task asks for one selected bit once.
"""

import argparse
import json
from math import sqrt
from pathlib import Path
import unittest

import numpy as np
import internal_prediction_memory as classical


I = np.eye(2)
Z = np.diag([1., -1.])
X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
OBSERVABLES = (Z, X)


def code_states(weight=.5):
    if not 0 <= weight <= 1:
        raise ValueError("Query weight must be a probability.")
    norm = sqrt(weight**2+(1-weight)**2)
    return {s: (I+(-1)**s[0]*weight/norm*Z+(-1)**s[1]*(1-weight)/norm*X)/2
            for s in classical.INPUTS}


def effect(query, outcome):
    return (I+(-1)**outcome*OBSERVABLES[query])/2


def mean_success(weight=.5):
    return sum(query_weight*np.trace(rho @ effect(query, state[query])).real/4
               for state, rho in code_states(weight).items()
               for query, query_weight in enumerate((weight, 1-weight)))


def forecast_data(weight=.5):
    states = code_states(weight)
    result = {}
    for query, query_weight in enumerate((weight, 1-weight)):
        for outcome in (0, 1):
            event_weights = {state: np.trace(rho @ effect(query, outcome)).real/4
                            for state, rho in states.items()}
            mass = sum(event_weights.values())
            posterior = sum(p for state, p in event_weights.items() if state[query])/mass
            result[(query, outcome)] = (query_weight*mass, posterior, event_weights)
    return result


def mean_brier_for_accuracy_code(weight=.5):
    return sum(mass*posterior*(1-posterior)
               for mass, posterior, _ in forecast_data(weight).values())


def optimized_preparations_score(first_observable, second_observable, weight):
    # Binary effect differences B_j obey -I <= B_j <= I. For fixed B's the best
    # preparation for each input is the largest-eigenvalue state of the sum.
    return .5+sum(np.linalg.eigvalsh(weight*(-1)**s[0]*first_observable
                                     +(1-weight)*(-1)**s[1]*second_observable)[-1]
                  for s in classical.INPUTS)/8


class RealPredictionCodeTests(unittest.TestCase):
    def test_real_states_and_effects_are_valid_and_need_no_y_component(self):
        for weight in (0., .2, .5, .8, 1.):
            for rho in code_states(weight).values():
                self.assertTrue(np.isrealobj(rho))
                np.testing.assert_allclose(rho @ rho, rho, atol=1e-14)
                self.assertAlmostEqual(np.trace(rho), 1)
                self.assertGreaterEqual(np.linalg.eigvalsh(rho).min(), -1e-14)
                self.assertAlmostEqual(abs(np.trace(rho @ Y)), 0)
            for query in (0, 1):
                np.testing.assert_array_equal(effect(query, 0)+effect(query, 1), I)
        np.testing.assert_allclose(sum(code_states().values()), 2*I, atol=1e-14)

    def test_equal_weight_code_has_uniform_success_for_all_input_query_pairs(self):
        expected = (1+1/sqrt(2))/2
        for state, rho in code_states().items():
            for query in (0, 1):
                self.assertAlmostEqual(np.trace(rho @ effect(query, state[query])), expected)
        self.assertGreater(mean_success(), .75)

    def test_weighted_real_construction_attains_general_qubit_upper_bound(self):
        for weight in np.linspace(0, 1, 21):
            bound = .5+.5*sqrt(weight**2+(1-weight)**2)
            self.assertAlmostEqual(mean_success(weight), bound)
            self.assertAlmostEqual(optimized_preparations_score(Z, X, weight), bound)

    def test_general_complex_biased_binary_effects_respect_analytic_bound(self):
        rng = np.random.default_rng(201)
        for _ in range(180):
            observables = []
            for _ in (0, 1):
                scalar = rng.uniform(-1, 1)
                vector = rng.normal(size=3)
                vector *= (1-abs(scalar))*rng.uniform(0, 1)/np.linalg.norm(vector)
                observables.append(scalar*I+vector[0]*X+vector[1]*Y+vector[2]*Z)
            weight = rng.uniform(0, 1)
            score = optimized_preparations_score(*observables, weight)
            self.assertLessEqual(score, .5+.5*sqrt(weight**2+(1-weight)**2)+1e-14)

    def test_both_bits_retrieval_differs_from_one_later_selected_bit(self):
        states = code_states()
        povm = {state: rho/2 for state, rho in states.items()}
        np.testing.assert_allclose(sum(povm.values()), I, atol=1e-14)
        full_pair_success = sum(np.trace(povm[state] @ rho).real/4 for state, rho in states.items())
        self.assertAlmostEqual(full_pair_success, .5)
        sequential_success = 0.
        for state, rho in states.items():
            first = effect(0, state[0])
            updated = first @ rho @ first
            sequential_success += np.trace(updated @ effect(1, state[1])).real/4
        self.assertAlmostEqual(sequential_success, mean_success()/2)
        self.assertLess(sequential_success, .5)

    def test_honest_forecasts_tie_classical_brier_for_the_symmetric_code(self):
        self.assertAlmostEqual(mean_brier_for_accuracy_code(), 1/8)
        self.assertAlmostEqual(mean_brier_for_accuracy_code(), float(classical.optimize_codes(2)[1]))
        for (query, _), (_, posterior, events) in forecast_data().items():
            squared_error = sum(p*(posterior-state[query])**2 for state, p in events.items())
            self.assertAlmostEqual(squared_error/sum(events.values()), posterior*(1-posterior))

    def test_accuracy_advantage_does_not_order_probability_forecast_loss(self):
        weight = .8
        self.assertGreater(mean_success(weight), .5+.5*weight)
        self.assertGreater(mean_brier_for_accuracy_code(weight), (1-weight)/4)
        # The classical strategy is itself also a valid orthogonal rebit code.
        # Therefore this is not a proof that optimal quantum Brier loss is worse.


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checked = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(RealPredictionCodeTests))
    if not checked.wasSuccessful():
        raise SystemExit(1)
    result = {
        "rounds": [201, 202],
        "scope": "Single isolated standard dimension-2 real/complex quantum memory, one later-selected bit; no entanglement or source-correlated side channel",
        "quantum_memory_and_born_rule_are_assumed_for_comparison": True,
        "uniform_classical_bit_optimal_mean_success": .75,
        "uniform_rebit_and_qubit_optimal_mean_success": mean_success(),
        "weighted_dimension_two_optimum": "1/2 + sqrt(w^2+(1-w)^2)/2",
        "y_component_needed": False,
        "exact_both_bits_from_one_copy_possible": False,
        "optimal_uniform_full_pair_success_dimension_two": .5,
        "symmetric_accuracy_code_brier": mean_brier_for_accuracy_code(),
        "symmetric_classical_bit_optimal_brier": .125,
        "global_optimal_quantum_brier_proved": False,
        "cognitive_necessity_of_one_bit_budget_or_quantum_structure_derived": False,
        "weighted_samples": [{"w": w, "classical_success": .5+max(w, 1-w)/2,
                               "rebit_success": mean_success(w),
                               "classical_optimal_brier": min(w, 1-w)/4,
                               "brier_of_accuracy_optimal_rebit_code": mean_brier_for_accuracy_code(w)}
                              for w in (.2, .5, .8)],
        "automated_tests": checked.testsRun,
    }
    if args.write_results:
        Path(__file__).with_name("real_prediction_code_results.json").write_text(
            json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"tests": checked.testsRun, "rounds": result["rounds"]}))


if __name__ == "__main__":
    main()
