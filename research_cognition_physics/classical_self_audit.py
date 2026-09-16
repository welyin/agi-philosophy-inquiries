"""Round 118: exact finite classical self-audit, with no hidden-world oracle."""
import argparse
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import unittest


@dataclass(frozen=True)
class Prediction:
    prior: F
    requested_action: int
    executed_action: int
    noise: F
    predicted_one: F
    model_version: str = "declared-binary-sensor-v1"
    source: str = "inference-after-declared-intervention"


@dataclass(frozen=True)
class Record:
    prediction: Prediction
    observation: int
    posterior: F
    observation_source: str = "sensor"


def probability(value):
    value = F(value)
    if not 0 <= value <= 1:
        raise ValueError("Probability outside [0,1].")
    return value


def commit(prior, noise=F(1,4), allow_toggle=True):
    prior, noise = probability(prior), probability(noise)
    requested = int(prior < F(1,2))
    action = requested if allow_toggle else 0
    propagated = 1-prior if action else prior
    return Prediction(prior, requested, action, noise,
                      noise + (1-2*noise)*propagated)


def close_record(prediction, observation):
    if observation not in (0,1):
        raise ValueError("Binary observation required.")
    propagated = 1-prediction.prior if prediction.executed_action else prediction.prior
    likelihood = 1-prediction.noise if observation else prediction.noise
    denominator = prediction.predicted_one if observation else 1-prediction.predicted_one
    if denominator == 0:
        raise ValueError("Observed event impossible under declared model: revise model.")
    return Record(prediction, observation, propagated*likelihood/denominator)


def audit(record):
    p = record.prediction
    if p.requested_action not in (0,1) or p.executed_action not in (0,p.requested_action):
        return False
    rebuilt = commit(p.prior, p.noise, allow_toggle=bool(p.executed_action))
    return (p == rebuilt and record == close_record(p, record.observation))


def joint_posterior(prior, likelihood):
    """Finite joint Bayesian update; correlations are retained, not factorized."""
    weights = {state: mass*likelihood[state] for state,mass in prior.items()}
    total = sum(weights.values())
    if total <= 0:
        raise ValueError("Impossible observation.")
    return {state: mass/total for state,mass in weights.items()}


def enumerate_histories(depth=8, actual_noise=F(1,4), declared_noise=F(1,4)):
    """Evaluator keeps world masses; the controller sees only its committed belief."""
    if depth < 0:
        raise ValueError("Nonnegative horizon required.")
    actual_noise, declared_noise = probability(actual_noise), probability(declared_noise)
    nodes = {(): ((F(2,3),F(1,3)), F(1,3))}
    audits = forecasts = 0
    max_forecast_gap = F(0)
    max_posterior_gap = F(0)
    expected_mistakes = F(0)
    maximum_bits = 0
    for step in range(depth):
        following = {}
        for history,(masses,belief) in nodes.items():
            pending = commit(belief, declared_noise, allow_toggle=(step % 3 != 1))
            propagated = masses[::-1] if pending.executed_action else masses
            mass = sum(propagated)
            true_one = (propagated[0]*actual_noise + propagated[1]*(1-actual_noise))/mass
            max_forecast_gap = max(max_forecast_gap, abs(true_one-pending.predicted_one))
            guess = int(pending.predicted_one >= F(1,2))
            expected_mistakes += mass*(1-true_one if guess else true_one)
            forecasts += 1
            for observation in (0,1):
                next_masses = tuple(propagated[x]*(1-actual_noise if x==observation else actual_noise)
                                    for x in (0,1))
                if sum(next_masses) == 0:
                    continue
                record = close_record(pending, observation)
                if not audit(record):
                    raise AssertionError("Internal audit failed.")
                audits += 1
                true_posterior = next_masses[1]/sum(next_masses)
                max_posterior_gap = max(max_posterior_gap, abs(true_posterior-record.posterior))
                maximum_bits = max(maximum_bits, record.posterior.numerator.bit_length(),
                                   record.posterior.denominator.bit_length())
                following[history+(observation,)] = (next_masses, record.posterior)
        nodes = following
    return {"horizon":depth, "terminal_histories":len(nodes), "forecasts":forecasts,
            "audited_records":audits, "terminal_mass":str(sum(sum(v[0]) for v in nodes.values())),
            "max_forecast_gap":str(max_forecast_gap), "max_posterior_gap":str(max_posterior_gap),
            "expected_incorrect_binary_guesses":str(expected_mistakes),
            "maximum_posterior_integer_bits":maximum_bits,
            "controller_has_hidden_world_access":False}


def storage_bound(horizon):
    """Specific prior 1/3 and noise 1/4. Lengths, program and model ID are separate."""
    if horizon < 0:
        raise ValueError("Nonnegative horizon required.")
    bits = (3*4**horizon).bit_length()
    return {"horizon":horizon, "one_integer_bits_upper":bits,
            "current_fraction_bits_upper":2*bits,
            "all_posterior_fraction_bits_upper":sum(2*(3*4**t).bit_length() for t in range(horizon+1)),
            "actions_requests_observations_bits":3*horizon,
            "extra_forecast_fraction_bits_upper":sum(2*(3*4**(t+1)).bit_length() for t in range(horizon))}


class ClassicalSelfAuditTests(unittest.TestCase):
    def test_every_history_is_exactly_calibrated_under_the_declared_model(self):
        report=enumerate_histories()
        self.assertEqual(report["terminal_histories"],256)
        self.assertEqual(report["audited_records"],510)
        self.assertEqual(report["terminal_mass"],"1")
        self.assertEqual(report["max_forecast_gap"],"0")
        self.assertEqual(report["max_posterior_gap"],"0")

    def test_veto_is_logged_before_forecast(self):
        self.assertEqual(commit(F(1,3)).predicted_one,F(7,12))
        veto=commit(F(1,3),allow_toggle=False)
        self.assertEqual((veto.requested_action,veto.executed_action),(1,0))
        self.assertEqual(veto.predicted_one,F(5,12))

    def test_tampered_prediction_or_posterior_fails_internal_audit(self):
        record=close_record(commit(F(1,3)),1)
        self.assertTrue(audit(record))
        self.assertFalse(audit(replace(record,posterior=F(1,2))))
        self.assertFalse(audit(replace(record,prediction=replace(record.prediction,predicted_one=F(1,2)))))
        self.assertFalse(audit(replace(record,observation_source="prediction")))

    def test_honest_wrong_model_is_not_external_truth(self):
        report=enumerate_histories(actual_noise=F(1,2))
        self.assertEqual(report["audited_records"],510)
        self.assertGreater(F(report["max_forecast_gap"]),0)
        self.assertGreater(F(report["max_posterior_gap"]),0)

    def test_correct_probability_does_not_promise_correct_single_outcome(self):
        self.assertGreater(F(enumerate_histories()["expected_incorrect_binary_guesses"]),0)

    def test_joint_model_preserves_relations_and_is_independent_of_grouping(self):
        states=list(product((0,1),repeat=3))
        prior={s:F(1,4) if sum(s)%2==0 else F(0) for s in states}
        likelihood={s:F(3,4) if s[0]==1 else F(1,4) for s in states}
        result=joint_posterior(prior,likelihood)
        regroup={(s[0],(s[1],s[2])):m for s,m in prior.items()}
        regroup_like={(s[0],(s[1],s[2])):m for s,m in likelihood.items()}
        updated=joint_posterior(regroup,regroup_like)
        for s in states:
            self.assertEqual(result[s],updated[(s[0],(s[1],s[2]))])
        for i in range(3):
            self.assertEqual(sum(m for s,m in prior.items() if s[i]),F(1,2))
        self.assertEqual(sum(m for s,m in prior.items() if sum(s)%2),0)

    def test_exact_memory_is_finite_and_charged(self):
        report=enumerate_histories()
        self.assertLessEqual(report["maximum_posterior_integer_bits"],storage_bound(8)["one_integer_bits_upper"])
        self.assertEqual(storage_bound(8)["current_fraction_bits_upper"],36)
        self.assertEqual(storage_bound(8)["all_posterior_fraction_bits_upper"],180)

    def test_impossible_events_require_model_revision_and_bad_inputs_fail(self):
        with self.assertRaises(ValueError): close_record(commit(F(1),noise=0),0)
        with self.assertRaises(ValueError): commit(F(2))
        with self.assertRaises(ValueError): joint_posterior({0:F(1)},{0:F(0)})
        with self.assertRaises(ValueError): enumerate_histories(-1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ClassicalSelfAuditTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":118,"matched_model":enumerate_histories(),
            "honest_but_misspecified_model":enumerate_histories(actual_noise=F(1,2)),
            "classical_storage_upper":storage_bound(8),
            "preparation_noncontextuality_derived":False,"Y_readout_needed":False,
            "scope":"four operational audit requirements; no claim to exhaust philosophical self-cognition",
            "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("classical_self_audit_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
