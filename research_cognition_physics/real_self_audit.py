"""Round 119: classical audit ledgers for the existing real instrument."""
import argparse
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product
import json
import math
from pathlib import Path
import unittest
import numpy as np

from quantum_interface_audit import branch_kraus, rotation_unitary, two_rebit_states
from operational_effect_closure import compiled_copy_gate


def trace_distance(a,b):
    return float(np.abs(np.linalg.eigvalsh(a-b)).sum()/2)


def instrument(history, outcome):
    """One old copy gate, two local rotations, then an original local readout."""
    step=len(history)
    angle=.17*step+(.23 if history and history[-1] else -.31)
    local=np.kron(rotation_unitary(angle),rotation_unitary(-.13*step))
    gate=compiled_copy_gate() @ local
    eta=(.3,.6,.9)[step%3]
    operators=[np.kron(k,np.eye(2)) @ gate for k in branch_kraus(angle,outcome,eta)]
    if max(np.max(np.abs(k.imag)) for k in operators)>1e-14:
        raise AssertionError("Unexpected complex physical operation.")
    return tuple(k.real for k in operators)


def branch(rho, history, outcome):
    return sum(k @ rho @ k.T for k in instrument(history,outcome))


@dataclass(frozen=True)
class QuantumPrediction:
    history: tuple
    model_entries: tuple
    predicted_one: float
    model_version: str = "real_rotation-original-instrument"
    capabilities: str = "joint real model; no direct local Y; no unknown-state oracle"


def commit(rho,history):
    return QuantumPrediction(tuple(history),tuple(map(float,rho.ravel())),
                             float(np.trace(branch(rho,history,1))))


def audit(pending,outcome,posterior):
    rho=np.array(pending.model_entries).reshape(4,4)
    rebuilt=commit(rho,pending.history)
    raw=branch(rho,pending.history,outcome)
    return (pending.model_version==rebuilt.model_version
            and pending.capabilities==rebuilt.capabilities
            and abs(pending.predicted_one-rebuilt.predicted_one)<1e-13
            and np.allclose(posterior,raw/np.trace(raw),atol=1e-12,rtol=0))


def run_tree(depth=6):
    # Both the true initial preparation and the declared predictive model are known here.
    initial=two_rebit_states()[0].real
    nodes={(): (initial,initial)}
    error=0.
    audits=0
    minimum_eigenvalue=0.
    for _ in range(depth):
        new={}
        for history,(raw,model) in nodes.items():
            pending=commit(model,history)  # happens before the observation
            for outcome in (0,1):
                actual=branch(raw,history,outcome)
                predicted_raw=branch(model,history,outcome)
                probability=float(np.trace(predicted_raw))
                updated=predicted_raw/probability
                actual_probability=float(np.trace(actual)/np.trace(raw))
                error=max(error,abs(probability-actual_probability),
                          float(np.max(np.abs(updated-actual/np.trace(actual)))))
                if not audit(pending,outcome,updated):
                    raise AssertionError("Invalid audit.")
                audits+=1
                minimum_eigenvalue=min(minimum_eigenvalue,float(np.linalg.eigvalsh(updated).min()))
                new[history+(outcome,)]=(actual,updated)
        nodes=new
    return {"horizon":depth,"terminal_histories":len(nodes),"audited_records":audits,
            "terminal_probability":float(sum(np.trace(v[0]) for v in nodes.values())),
            "max_numerical_consistency_residual":error,
            "minimum_numerical_eigenvalue":minimum_eigenvalue,
            "resources_per_executed_path":{"pair_rotations":depth,
                                          "local_rotations":5*depth,
                                          "original_readouts":depth,
                                          "classical_outcome_bits":depth},
            "numeric_values_are_diagnostics_not_interval_certificates":True}


def rational_model(rho,bits):
    """Real PSD rational representative with trace distance <= d^2 * 2^-bits.
    Certificate assumes exact entry rounding; float input is the diagnostic target.
    """
    if bits<1 or int(bits)!=bits:
        raise ValueError("Positive integer precision required.")
    if not np.allclose(rho,rho.T,atol=1e-14) or np.max(np.abs(np.imag(rho)))>1e-14:
        raise ValueError("Real symmetric model required.")
    d=len(rho)
    scale=2**bits
    rounded=[[F(round(float(np.real(rho[i,j]))*scale),scale) for j in range(d)] for i in range(d)]
    shift=F(d,2*scale)
    for i in range(d): rounded[i][i]+=shift
    normalizer=sum(rounded[i][i] for i in range(d))
    return tuple(tuple(value/normalizer for value in row) for row in rounded)


def model_storage(d,bits):
    if d<1 or bits<1: raise ValueError("Positive dimensions and precision required.")
    return {"independent_real_entries":d*(d+1)//2,
            "rounded_model_payload_bits":d*(d+1)//2*(bits+2),
            "trace_distance_upper_exact":str(F(d*d,2**bits)),
            "metadata_and_arithmetic_work_space_not_included":True}


def conditional_distance_bound(unnormalized_trace_norm_error, probability):
    error,prob=F(unnormalized_trace_norm_error),F(probability)
    if error<0 or not 0<prob<=1: raise ValueError("Invalid error or event probability.")
    return min(F(1),error/prob)


def rare_record_example(p):
    p=F(p)
    if not 0<p<1: raise ValueError("Use 0<p<1.")
    # First register is an already classical, orthogonal history flag.
    rho=np.diag([float(1-p),0,float(p),0])
    sigma=np.diag([float(1-p),0,0,float(p)])
    return rho,sigma,rho[2:,2:]/float(p),sigma[2:,2:]/float(p)


class RealSelfAuditTests(unittest.TestCase):
    def test_all_original_instrument_records_have_consistent_real_audits(self):
        report=run_tree()
        self.assertEqual(report["audited_records"],126)
        self.assertEqual(report["terminal_histories"],64)
        self.assertAlmostEqual(report["terminal_probability"],1,places=13)
        self.assertLess(report["max_numerical_consistency_residual"],2e-14)
        self.assertGreater(report["minimum_numerical_eigenvalue"],-2e-14)

    def test_adaptive_instruments_preserve_probability_and_are_real(self):
        for history in ((),(0,),(1,0),(1,1,0)):
            operators=[k for outcome in (0,1) for k in instrument(history,outcome)]
            np.testing.assert_allclose(sum(k.T@k for k in operators),np.eye(4),atol=2e-15,rtol=0)

    def test_model_register_tampering_is_detected(self):
        rho=two_rebit_states()[0].real
        pending=commit(rho,())
        raw=branch(rho,(),1)
        self.assertTrue(audit(pending,1,raw/np.trace(raw)))
        self.assertFalse(audit(replace(pending,predicted_one=.99),1,raw/np.trace(raw)))
        self.assertFalse(audit(pending,1,np.eye(4)/4))

    def test_rational_positive_models_have_the_proved_error_budget(self):
        rng=np.random.default_rng(119)
        for d,bits in product((2,4,8),(8,16)):
            a=rng.normal(size=(d,d))
            rho=a@a.T
            rho/=np.trace(rho)
            rational=rational_model(rho,bits)
            self.assertEqual(sum(rational[i][i] for i in range(d)),F(1))
            sigma=np.array(rational,dtype=float)
            self.assertGreaterEqual(np.linalg.eigvalsh(sigma).min(),-1e-14)
            self.assertLessEqual(trace_distance(rho,sigma),float(F(d*d,2**bits)))
        self.assertEqual(model_storage(4,20)["rounded_model_payload_bits"],220)

    def test_complete_record_distance_is_bounded_by_input_model_distance(self):
        rho=two_rebit_states()[0].real
        sigma=np.array(rational_model(rho,12),dtype=float)
        nodes=[(rho,sigma)]
        for depth in range(5):
            nodes=[(branch(a,h,o),branch(b,h,o))
                   for index,(a,b) in enumerate(nodes)
                   for o in (0,1)
                   for h in [tuple(map(int,format(index,f"0{depth}b"))) if depth else ()]]
        tv=sum(abs(np.trace(a)-np.trace(b)) for a,b in nodes)/2
        self.assertLessEqual(tv,trace_distance(rho,sigma)+1e-14)

    def test_small_unconditional_error_can_hide_large_conditional_error(self):
        for p in (F(1,100),F(1,10**8)):
            rho,sigma,a,b=rare_record_example(p)
            self.assertAlmostEqual(trace_distance(rho,sigma),float(p),places=15)
            self.assertAlmostEqual(trace_distance(a,b),1)
            self.assertEqual(conditional_distance_bound(2*p,p),1)

    def test_conditional_normalization_bound_on_random_positive_matrices(self):
        rng=np.random.default_rng(120)
        for _ in range(20):
            a=rng.normal(size=(4,4));a=a@a.T/30
            b=rng.normal(size=(4,4));b=b@b.T/30
            p,q=np.trace(a),np.trace(b)
            error=np.abs(np.linalg.eigvalsh(a-b)).sum()
            self.assertLessEqual(trace_distance(a/p,b/q),error/p+1e-14)

    def test_joint_relation_survives_classical_description_without_local_Y(self):
        plus,minus=(x.real for x in two_rebit_states())
        np.testing.assert_allclose(plus.reshape(2,2,2,2).trace(axis1=1,axis2=3),np.eye(2)/2)
        np.testing.assert_allclose(minus.reshape(2,2,2,2).trace(axis1=1,axis2=3),np.eye(2)/2)
        self.assertEqual(trace_distance(plus,minus),1)
        self.assertNotEqual(commit(plus,()).model_entries,commit(minus,()).model_entries)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RealSelfAuditTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":119,"original_real_instrument_tree":run_tree(),
            "finite_real_model_example":model_storage(4,20),
            "rare_record_example":{"probability_exact":"1/100000000","initial_trace_distance_exact":"1/100000000",
                                   "conditional_trace_distance_exact":"1"},
            "unknown_state_readout_used":False,"native_local_Y_used":False,
            "general_claim":"append a finite classical audit to a specified computable real model; not exact finite tomography",
            "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("real_self_audit_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
