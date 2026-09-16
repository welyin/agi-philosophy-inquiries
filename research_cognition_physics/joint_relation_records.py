"""Round 98: complete noisy records and the orthogonal-observable joint-read bound."""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from bilocal_record_tomography import embed_operator
from certified_intervals import Interval as I, SCALE
from complex_whole_real_interfaces import random_state
from incompatible_relation_writes import (G,H,sequence_isometry,written_state,
    forget_memories,old_channel,calibration_state)
from independent_source_alignment import keep_systems
from one_bit_real_network import visibility_interval
from operational_effect_closure import compiled_copy_gate, majority_effect_by_circuit
from quantum_interface_audit import (ALPHA, PAULI_X, PAULI_Z, branch_kraus,
    effective_readout, rotation_unitary)
from weak_relation_tradeoff import read_axis


def memory_axis(c):
    return c*PAULI_X-math.sqrt((1-c)*(1+c))*PAULI_Z


def validated_count(count):
    if not isinstance(count,int) or count<1 or count%2==0:
        raise ValueError("Use a positive odd read count.")


def pointer_effect(c,count,sign):
    validated_count(count)
    if sign not in (-1,1): raise ValueError("Use sign +/-1.")
    axis=read_axis(2*math.acos(c))
    if count==1:
        return effective_readout(axis,int(sign==1),1)
    rotate=rotation_unitary(axis)
    plus=rotate@majority_effect_by_circuit(count)@rotate.conj().T
    return plus if sign==1 else np.eye(2)-plus


def compiled_joint_effect(c,d,m,n,r,t,reverse=False):
    v=sequence_isometry(c,d,reverse)
    effect=np.kron(np.kron(pointer_effect(c,m,r),pointer_effect(d,n,t)),np.eye(4))
    return v.conj().T@effect@v


def analytic_joint_effect(c,d,gamma,eta,r,t,reverse=False):
    v,w=math.sqrt((1-c)*(1+c)),math.sqrt((1-d)*(1+d))
    x,y=(gamma*d*v,eta*w) if reverse else (gamma*v,eta*c*w)
    return (np.eye(4)+r*x*G+t*y*H)/4


def read_memory_tree(state,site,c,count):
    return axis_read_memory_tree(state,site,read_axis(2*math.acos(c)),count)


def axis_read_memory_tree(state,site,axis,count):
    """Actual finite Kraus branches. A protected read uses one reset pointer."""
    validated_count(count)
    systems=int(round(math.log2(len(state))))
    if count==1:
        return {(r,):sum((full:=embed_operator(k,(site,),systems))@state@full.conj().T
                          for k in branch_kraus(axis,r,1.)) for r in (0,1)}
    rotate=embed_operator(rotation_unitary(-axis),(site,),systems)
    leaves={():rotate@state@rotate.conj().T}
    copier=embed_operator(compiled_copy_gate(),(0,site+1),systems+1)
    noise={r:[np.kron(k,np.eye(len(state))) for k in branch_kraus(0,r,1.)] for r in (0,1)}
    for _ in range(count):
        following={}
        for history,branch in leaves.items():
            copied=copier@np.kron(np.diag([1.,0]),branch)@copier.conj().T
            for r in (0,1):
                measured=sum(k@copied@k.conj().T for k in noise[r])
                following[history+(r,)]=keep_systems(measured,tuple(range(1,systems+1)),systems+1)
        leaves=following
    return leaves


def complete_record_tree(state,c,d,m,n):
    after=written_state(state,c,d)
    return {(first,second):branch2
            for first,branch in read_memory_tree(after,0,c,m).items()
            for second,branch2 in read_memory_tree(branch,1,d,n).items()}


def complete_record_probability(state,c,d,first,second):
    v,w=math.sqrt((1-c)*(1+c)),math.sqrt((1-d)*(1+d))
    g,h=np.trace(G@state).real,np.trace(H@state).real
    return sum((1+r*v*g+t*c*w*h)/4
        *math.prod((1+(2*x-1)*r*ALPHA)/2 for x in first)
        *math.prod((1+(2*x-1)*t*ALPHA)/2 for x in second)
        for r,t in product((-1,1),repeat=2))


def interval_fields(value):
    return {"lower_exact":str(Fraction(value.lo,SCALE)),
            "upper_exact":str(Fraction(value.hi,SCALE)),
            "diagnostic":sum(value.floats())/2}


def record_certificate(c,d,m,n):
    c,d=Fraction(c),Fraction(d)
    if not 0<=c<=1 or not 0<=d<=1: raise ValueError("Overlaps must be in [0,1].")
    validated_count(m); validated_count(n)
    v,w=(1-I.exact(c*c)).sqrt(),(1-I.exact(d*d)).sqrt()
    x=visibility_interval(m)*v
    y=visibility_interval(n)*I.exact(c)*w
    return {"c_exact":str(c),"d_exact":str(d),"G_read_count":m,"H_read_count":n,
        "G_calibration_success":interval_fields((1+x)/2),
        "H_calibration_success":interval_fields((1+y)/2),
        "ideal_joint_sharpness_squared_sum_exact":str(1-c*c*d*d),
        "finite_joint_sharpness_squared_sum":interval_fields(x*x+y*y),
        "interpretation":"Two separate eigenstate calibration ensembles; not simultaneous hidden G,H values"}


class JointRelationRecordsTests(unittest.TestCase):
    def test_actual_pointer_circuits_give_all_four_effects_in_both_orders(self):
        for c,d in ((.8,.8),(.3,.6),(0,0),(1,.4)):
            for m,n in ((1,1),(3,1),(1,3)):
                gamma=sum(visibility_interval(m).floats())/2
                eta=sum(visibility_interval(n).floats())/2
                for reverse in (False,True):
                    for r,t in product((-1,1),repeat=2):
                        np.testing.assert_allclose(compiled_joint_effect(c,d,m,n,r,t,reverse),
                            analytic_joint_effect(c,d,gamma,eta,r,t,reverse),atol=2e-15)

    def test_full_actual_record_tree_matches_latent_pair_formula(self):
        state=random_state(np.random.default_rng(98),4)
        for m,n in ((1,1),(3,1),(3,3)):
            leaves=complete_record_tree(state,.8,.6,m,n)
            self.assertEqual(len(leaves),2**(m+n))
            for (first,second),branch in leaves.items():
                self.assertAlmostEqual(np.trace(branch).real,
                    complete_record_probability(state,.8,.6,first,second),places=14)
            np.testing.assert_allclose(sum(forget_memories(b) for b in leaves.values()),
                old_channel(state,.8,.6),atol=2e-15)

    def test_each_calibration_success_matches_certified_formula(self):
        row=record_certificate(Fraction(4,5),Fraction(4,5),1,1)
        for op,index,key in ((G,0,"G_calibration_success"),(H,1,"H_calibration_success")):
            success=0.
            for sign in (-1,1):
                leaves=complete_record_tree(calibration_state(op,sign),.8,.8,1,1)
                success+=sum(np.trace(b).real for history,b in leaves.items()
                             if 2*history[index][0]-1==sign)/2
            self.assertAlmostEqual(success,row[key]["diagnostic"],places=14)

    def test_effects_are_normalized_and_have_no_product_record_term(self):
        effects={(r,t):analytic_joint_effect(.8,.6,.9,.7,r,t) for r,t in product((-1,1),repeat=2)}
        np.testing.assert_allclose(sum(effects.values()),np.eye(4),atol=1e-16)
        np.testing.assert_allclose(sum(r*t*f for (r,t),f in effects.items()),0,atol=1e-16)
        for effect in effects.values(): self.assertGreater(np.linalg.eigvalsh(effect).min(),0)

    def test_entire_ideal_disk_is_realizable_by_sequential_strengths(self):
        for x,y in ((0,1),(.6,.8),(.3,.4),(1,0),(0,0)):
            c=math.sqrt(max(0,1-x*x))
            d=math.sqrt(max(0,1-(y/c)**2)) if c else 1.
            np.testing.assert_allclose(analytic_joint_effect(c,d,1,1,1,1),
                (np.eye(4)+x*G+y*H)/4,atol=3e-16)

    def test_universal_positivity_proof_ingredients_for_arbitrary_complex_effects(self):
        rng=np.random.default_rng(198)
        for _ in range(20):
            raw=[random_state(rng,4) for _ in range(4)]
            vals,vecs=np.linalg.eigh(sum(raw))
            inverse=(vecs/np.sqrt(vals))@vecs.conj().T
            effects=[inverse@f@inverse for f in raw]
            coords=[]
            for (r,t),f in zip(product((-1,1),repeat=2),effects):
                p=np.trace(f).real/4
                u,v=np.trace(G@f).real/4,np.trace(H@f).real/4
                self.assertLessEqual(math.hypot(u,v),p+1e-15)
                coords.append([r*u,t*v])
            self.assertLessEqual(np.linalg.norm(np.sum(coords,axis=0)),1+1e-15)

    def test_exact_joint_bound_rejects_two_perfect_original_sign_records(self):
        self.assertLess(np.linalg.eigvalsh((np.eye(4)+G+H)/4).min(),0)
        for x in (0,.3,.6,1):
            y=math.sqrt(1-x*x)
            for r,t in product((-1,1),repeat=2):
                effect=(np.eye(4)+r*x*G+t*y*H)/4
                self.assertGreaterEqual(np.linalg.eigvalsh(effect).min(),-1e-15)

    def test_finite_certificates_have_strict_joint_margin_and_improve_under_rereading(self):
        earlier=None
        for m in (1,3,5):
            row=record_certificate(Fraction(4,5),Fraction(4,5),m,m)
            self.assertLess(Fraction(row["finite_joint_sharpness_squared_sum"]["upper_exact"]),
                            Fraction(row["ideal_joint_sharpness_squared_sum_exact"]))
            current=Fraction(row["H_calibration_success"]["lower_exact"])
            if earlier is not None: self.assertGreater(current,earlier)
            earlier=Fraction(row["H_calibration_success"]["upper_exact"])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(JointRelationRecordsTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":98,"forward_effect":"F_rt=(I+r*gamma_m*sqrt(1-c^2)*G+t*gamma_n*c*sqrt(1-d^2)*H)/4",
        "universal_joint_unbiased_marginal_condition":"x^2+y^2<=1 (necessary and sufficient)",
        "ideal_sequential_disk_covered":True,"original_joint_hidden_eigenvalues_assumed":False,
        "examples":[record_certificate(Fraction(4,5),Fraction(4,5),m,m) for m in (1,3,5)],
        "universal_equal_resource_optimality_claimed":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("joint_relation_records_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
