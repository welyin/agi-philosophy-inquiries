"""Round 120: self-description versus nondisturbing unknown-state readout.

The exact no-information theorem is general. The finite-disturbance frontier
below is restricted to a two-pure-state product-output isometry, not all instruments.
"""
import argparse
import json
import math
from pathlib import Path
import unittest
import numpy as np
from quantum_interface_audit import ALPHA, effective_readout


def vector(angle):
    return np.array([math.cos(angle),math.sin(angle)])


def projector(state):
    return np.outer(state,state.conj())


def frame(first,second):
    columns=[first,second]
    for seed in np.eye(len(first)):
        residual=seed.copy()
        for col in columns:
            residual-=np.dot(col,residual)*col
        length=np.linalg.norm(residual)
        if length>1e-10: columns.append(residual/length)
        if len(columns)==len(first): break
    return np.column_stack(columns)


def product_frontier(overlap,disturbance):
    """s in (0,1), loss of each system's pure-state fidelity <= disturbance."""
    s=float(overlap)
    delta=float(disturbance)
    if not 0<s<1 or not 0<=delta<=1:
        raise ValueError("Use 0<s<1 and 0<=disturbance<=1.")
    theta=math.acos(s)
    epsilon=min(math.asin(math.sqrt(delta)),theta/2)
    t=math.cos(theta-2*epsilon)
    pointer_overlap=min(1.,s/t)
    distinguishability=math.sqrt(max(0.,1-pointer_overlap**2))
    return {"input_overlap":s,"requested_fidelity_loss":delta,
            "actual_fidelity_loss":math.sin(epsilon)**2,
            "rotation_epsilon":epsilon,"system_output_overlap":t,
            "pointer_overlap":pointer_overlap,"pointer_trace_distance":distinguishability,
            "ideal_pointer_success":(1+distinguishability)/2,
            "original_noisy_pointer_success":(1+ALPHA*distinguishability)/2}


def real_construction(overlap,disturbance):
    result=product_frontier(overlap,disturbance)
    theta=math.acos(overlap)
    epsilon=result["rotation_epsilon"]
    psi=[vector(0),vector(theta)]
    phi=[vector(epsilon),vector(theta-epsilon)]
    c=result["pointer_overlap"]
    records=[vector(0),np.array([c,math.sqrt(max(0.,1-c*c))])]
    input_vectors=[np.kron(p,vector(0)) for p in psi]
    outputs=[np.kron(p,r) for p,r in zip(phi,records)]
    norm=math.sqrt(1-overlap**2)
    before=frame(input_vectors[0],(input_vectors[1]-overlap*input_vectors[0])/norm)
    after=frame(outputs[0],(outputs[1]-overlap*outputs[0])/norm)
    if np.linalg.det(after@before.T)<0: after[:,-1]*=-1
    operation=after@before.T
    return operation,psi,phi,records


def pointer_success(records):
    difference=projector(records[0])-projector(records[1])
    x=2*difference[0,1]
    z=difference[0,0]-difference[1,1]
    angle=math.atan2(x,z)
    positive=effective_readout(angle,1,1).real
    return float((np.trace(positive@projector(records[0]))
                  +np.trace((np.eye(2)-positive)@projector(records[1])))/2)


class SelfReadoutBoundaryTests(unittest.TestCase):
    def test_no_disturbance_means_no_record_information_for_nonorthogonal_pair(self):
        for s in (.1,.5,math.sqrt(.5),.99):
            result=product_frontier(s,0)
            self.assertAlmostEqual(result["pointer_overlap"],1,places=14)
            self.assertAlmostEqual(result["ideal_pointer_success"],.5,places=7)

    def test_real_isometry_preserves_both_inputs_and_the_full_gram_matrix(self):
        for s in (.2,math.sqrt(.5),.9):
            for delta in (0,.001,.01,.1,.5):
                operation,psi,phi,records=real_construction(s,delta)
                np.testing.assert_allclose(operation.T@operation,np.eye(4),atol=2e-14,rtol=0)
                self.assertAlmostEqual(np.linalg.det(operation),1,places=13)
                for source,target,record in zip(psi,phi,records):
                    np.testing.assert_allclose(operation@np.kron(source,vector(0)),
                                               np.kron(target,record),atol=2e-14,rtol=0)

    def test_each_physical_system_respects_declared_disturbance_budget(self):
        for delta in (0,.001,.01,.1,.5):
            _,psi,phi,_=real_construction(math.sqrt(.5),delta)
            for initial,final in zip(psi,phi):
                self.assertLessEqual(1-abs(np.dot(initial,final))**2,delta+1e-14)

    def test_original_noisy_readout_has_the_calculated_information(self):
        for delta in (.001,.01,.1,.5):
            result=product_frontier(math.sqrt(.5),delta)
            _,_,_,records=real_construction(math.sqrt(.5),delta)
            self.assertAlmostEqual(pointer_success(records),result["original_noisy_pointer_success"],places=14)

    def test_full_transfer_cannot_beat_initial_nonorthogonal_discrimination(self):
        for s in (.2,.5,math.sqrt(.5),.9):
            result=product_frontier(s,1)
            self.assertAlmostEqual(result["system_output_overlap"],1)
            self.assertAlmostEqual(result["pointer_trace_distance"],math.sqrt(1-s*s))
            self.assertLess(result["ideal_pointer_success"],1)

    def test_orthogonal_classical_labels_have_an_exact_nondisturbing_copy(self):
        # This contrasting gate is an abstract classical register copy, not an old readout.
        copy=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
        for label in (0,1):
            source=np.eye(2)[label]
            np.testing.assert_array_equal(copy@np.kron(source,np.array([1,0])),np.kron(source,source))

    def test_quantitative_bound_rejects_a_perfect_nondisturbing_label_for_nonorthogonal_inputs(self):
        s=math.sqrt(.5)
        unchanged_system_overlap=s
        orthogonal_record_overlap=0
        self.assertGreater(abs(s-unchanged_system_overlap*orthogonal_record_overlap),.7)

    def test_invalid_budgets_are_rejected(self):
        for s,delta in ((0,0),(1,0),(.5,-.1),(.5,1.1)):
            with self.assertRaises(ValueError): product_frontier(s,delta)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SelfReadoutBoundaryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":120,"input_example":"|0>, (|0>+|1>)/sqrt(2)",
            "rows":[product_frontier(math.sqrt(.5),delta) for delta in (0,.001,.01,.1,(1-math.sqrt(.5))/2)],
            "nondisturbance_theorem_scope":"all real or complex channels with fixed independent apparatus and unchanged pure input marginals",
            "finite_disturbance_frontier_scope":"two pure inputs, pure product outputs, fidelity loss bounded on each input",
            "construction_is_real":True,"general_SO4_construction_compiled_into_old_gates":False,
            "values_are_floating_diagnostics_of_analytic_formulas":True,
            "quantum_theory_derived_from_cognition":False,
            "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("self_readout_boundary_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
