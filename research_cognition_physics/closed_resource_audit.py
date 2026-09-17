"""Round 205: resource closure, using the actual round-122/123 record law.

Capacity and entropy are different quantities; neither is an energy price.
Numerical entropy values are diagnostics. The support and rate bounds are
analytic, under the fixed diagnostic setting and fresh-target assumptions.
"""
import argparse
from fractions import Fraction as F
from itertools import product
import json
import math
from pathlib import Path
import unittest

import numpy as np

from operational_effect_closure import majority_certificate
from reference_history_audit import exact_posterior, signed_outcome
from reusable_network_reference import resources, sector_tables


RECORDS = tuple(product((-1, 1), range(4), (-1, 1)))


def binary_entropy(p):
    if not 0 <= p <= 1:
        raise ValueError("Probability outside [0,1].")
    return -math.fsum(t * math.log2(t) for t in (p, 1-p) if t > 0)


def entropy(probabilities):
    return -math.fsum(p * math.log2(p) for p in probabilities if p > 0)


def diagnostic_visibility(count=5):
    gamma = majority_certificate(count)["effective_visibility_diagnostic"]
    return gamma**3 / math.sqrt(2)


def analytic_record(sector, record, visibility):
    return (1 + sector * visibility * signed_outcome(*record)) / 16


def actual_history_probability(history, g, count=5):
    tables = sector_tables(count)
    return math.fsum((1+s*g)/2 * math.prod(tables[s][(1, 4)+r] for r in history)
                     for s in (-1, 1))


def history_entropy(n, g, visibility):
    """O(n) diagnostic via the sufficient count; no enumeration of 16**n logs."""
    if not isinstance(n, int) or n < 0 or not -1 <= g <= 1 or not 0 <= visibility < 1:
        raise ValueError("Require integer n>=0, |g|<=1 and 0<=visibility<1.")
    prior = (1+g)/2
    p = (1+visibility)/2
    rate = 3 + binary_entropy(p)
    if n == 0 or visibility == 0 or abs(g) == 1:
        mutual = 0.
        mass = 1.
    else:
        terms, masses = [], []
        for k in range(n+1):
            choose = math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
            plus = math.log(prior)+choose+k*math.log(p)+(n-k)*math.log1p(-p)
            minus = math.log1p(-prior)+choose+k*math.log1p(-p)+(n-k)*math.log(p)
            peak = max(plus, minus)
            a, b = math.exp(plus-peak), math.exp(minus-peak)
            weight = math.exp(peak)*(a+b)
            masses.append(weight)
            terms.append(weight*binary_entropy(a/(a+b)))
        mass = math.fsum(masses)
        mutual = binary_entropy(prior)-math.fsum(terms)
    return {"trials": n, "sector_prior_bias": g,
            "conditional_entropy_bits_per_trial": rate,
            "sector_information_bits": mutual,
            "sector_information_upper_bits": binary_entropy(prior),
            "record_entropy_bits": n*rate+mutual,
            "count_distribution_mass_diagnostic": mass,
            "fixed_capacity_exact_archive_bits": 4*n,
            "classical_count_summary_bits_at_known_n": n.bit_length(),
            "count_summary_preserves_only_declared_future_predictions": True}


def pointer_cycle(seed, pointer, log):
    """Three reversible XOR gates; on blank pointer/log, copy seed to log.

    This is a classical resource example, not a replacement for old noisy
    quantum pointer instruments. The seed and log remain inside the whole.
    """
    pointer ^= seed
    log ^= pointer
    pointer ^= seed
    return seed, pointer, log


def report():
    gamma = majority_certificate(5)["effective_visibility_diagnostic"]
    v = diagnostic_visibility()
    return {"round": 205,
            "resource_accounting_requirement_accepted_by_user": True,
            "physical_axiom_selecting_quantum_theory_adopted": False,
            "record_interface": {"setting": [1, 4], "outcomes_per_trial": 16,
                "raw_readouts_are_not_all_archived_in_this_bound": True,
                "fresh_targets_conditionally_independent_given_fixed_sector": True,
                "gamma_diagnostic": gamma, "visibility_diagnostic": v,
                "smallest_single_record_probability_diagnostic": (1-v)/16,
                "exact_full_archive_dimension_lower_bound": "16**N",
                "entropy_identity_bits": "H(R^N)=N*(3+h2((1+v)/2))+I(S;R^N)",
                "mutual_information_upper_bits": "h2((1+g)/2)",
                "count_summary_is_round_123_result_reused": True},
            "examples": [history_entropy(n, gamma, v) for n in (1, 10, 1000)],
            "legacy_protocol_budget_for_1000": resources(1000),
            "resource_closure_status": {
                "reference_reuse_proved": True,
                "whole_apparatus_cyclic_return_proved": False,
                "complete_autonomous_physical_implementation_supplied": False,
                "missing_physical_accounts": ["target source preparation and recovery",
                    "pointer preparation, readout dilation and reset",
                    "controller, setting program, clock and communication hardware",
                    "storage and noise environment",
                    "Hamiltonians, work store, heat reservoir and operating time"],
                "resource_counts_may_share_physical_carriers_do_not_add_blindly": True,
                "marginal_reference_return_does_not_imply_uncorrelated_return": True},
            "scope": {"fixed_capacity_and_full_recoverable_history_required_for_4N_bound": True,
                "all_cognition_requires_full_history_retention": False,
                "entropy_equals_heat": False,
                "closed_universe_assumed_finite_dimensional": False,
                "real_complex_matched_record_laws_have_same_archive_bound": True,
                "ordinary_real_quantum_excluded": False,
                "complex_quantum_necessity_derived": False}}


class ClosedResourceAuditTests(unittest.TestCase):
    def test_entropy_uses_actual_old_matrix_record_law(self):
        for count in (1, 3, 5):
            v = diagnostic_visibility(count)
            tables = sector_tables(count)
            for s in (-1, 1):
                actual = [tables[s][(1, 4)+r] for r in RECORDS]
                predicted = [analytic_record(s, r, v) for r in RECORDS]
                np.testing.assert_allclose(actual, predicted, atol=2e-16, rtol=0)
                self.assertAlmostEqual(entropy(actual), 3+binary_entropy((1+v)/2), places=13)

    def test_enumerated_joint_histories_match_entropy_identity(self):
        v = diagnostic_visibility()
        for n, g in product((1, 2, 3), (-1., 0., .61, 1.)):
            probabilities = [actual_history_probability(h, g) for h in product(RECORDS, repeat=n)]
            self.assertAlmostEqual(math.fsum(probabilities), 1., places=13)
            self.assertAlmostEqual(entropy(probabilities), history_entropy(n, g, v)["record_entropy_bits"], places=12)

    def test_all_histories_have_positive_probability_even_with_known_sector(self):
        v = diagnostic_visibility()
        for g in (-1., -.2, 0., .9, 1.):
            probabilities = [actual_history_probability(h, g) for h in product(RECORDS, repeat=2)]
            self.assertEqual(len(probabilities), 16**2)
            self.assertGreaterEqual(min(probabilities), ((1-v)/16)**2-1e-16)

    def test_long_batch_bounds_and_binomial_normalization(self):
        for n, g in product((10, 1000), (0., .61, .999997222423702)):
            result = history_entropy(n, g, diagnostic_visibility())
            self.assertAlmostEqual(result["count_distribution_mass_diagnostic"], 1., places=11)
            self.assertGreaterEqual(result["sector_information_bits"], -1e-12)
            self.assertLessEqual(result["sector_information_bits"], result["sector_information_upper_bits"]+1e-12)
            self.assertLess(result["record_entropy_bits"], 4*n)

    def test_zero_information_boundaries(self):
        for n, g, v in ((0, .3, .7), (4, 1., .7), (4, -1., .7), (4, .3, 0.)):
            result = history_entropy(n, g, v)
            self.assertEqual(result["sector_information_bits"], 0.)
        self.assertEqual(history_entropy(4, .3, 0.)["record_entropy_bits"], 16.)

    def test_shared_sector_information_is_bounded_while_new_record_entropy_grows(self):
        result = history_entropy(1000, 0., diagnostic_visibility())
        self.assertAlmostEqual(result["sector_information_bits"], 1., places=13)
        self.assertGreater(result["record_entropy_bits"], 3600.)

    def test_old_count_is_sufficient_for_fixed_future_prediction(self):
        g, v = F(3, 10), F(2, 3)
        for history in product((-1, 1), repeat=5):
            plus = (1+g)/2*math.prod((1+v*w)/2 for w in history)
            minus = (1-g)/2*math.prod((1-v*w)/2 for w in history)
            posterior = (plus-minus)/(plus+minus)
            self.assertEqual(posterior, exact_posterior(g, v, 5, history.count(1)))

    def test_exact_classical_predictor_has_n_plus_one_distinct_values(self):
        for n in (1, 7, 20):
            values = {exact_posterior(F(3, 10), F(2, 3), n, k) for k in range(n+1)}
            self.assertEqual(len(values), n+1)
        self.assertEqual(history_entropy(1000, .3, .7)["classical_count_summary_bits_at_known_n"], 10)

    def test_pointer_can_cycle_reversibly_when_record_and_source_remain(self):
        inputs = tuple(product((0, 1), repeat=3))
        outputs = {pointer_cycle(*state) for state in inputs}
        self.assertEqual(outputs, set(inputs))
        for seed in (0, 1):
            self.assertEqual(pointer_cycle(seed, 0, 0), (seed, 0, seed))

    def test_local_marginal_entropies_do_not_make_an_additive_resource_account(self):
        # Fair seed, blank pointer, blank log -> correlated seed/log, blank pointer.
        distribution = {pointer_cycle(seed, 0, 0): .5 for seed in (0, 1)}
        joint = entropy(distribution.values())
        marginals = []
        for index in range(3):
            probabilities = [sum(p for state, p in distribution.items() if state[index] == bit) for bit in (0, 1)]
            marginals.append(entropy(probabilities))
        self.assertEqual(joint, 1.)
        self.assertEqual(marginals, [1., 0., 1.])
        self.assertEqual(sum(marginals)-joint, 1.)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ClosedResourceAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result["automated_checks"] = {"run": checks.testsRun, "failures": 0, "errors": 0}
    if args.write_results:
        Path(__file__).with_name("closed_resource_audit_results.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
