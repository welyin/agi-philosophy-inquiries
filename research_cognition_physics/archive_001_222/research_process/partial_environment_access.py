"""Round 103: exact query bounds with only part of a redundant environment."""
import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

from bilocal_record_tomography import embed_operator
from certified_intervals import Interval as I
from coherent_memory_compression import vector
from complex_whole_real_interfaces import random_state
from deferred_relation_queries import memory_channel
from incompatible_relation_writes import G,H,calibration_state,written_state,forget_memories,old_channel
from independent_source_alignment import keep_systems
from joint_relation_records import memory_axis,interval_fields
from memory_dephasing_threshold import dephasing_unitary,dephase_memory,valid_overlap
from weak_relation_tradeoff import trace_distance


def partition_overlaps(overlaps,accessible):
    overlaps=tuple(overlaps)
    accessible=tuple(accessible)
    for x in overlaps: valid_overlap(x)
    if len(set(accessible))!=len(accessible) or any(
            not isinstance(i,int) or i<0 or i>=len(overlaps) for i in accessible):
        raise ValueError("Accessible indices must be distinct and in range.")
    return (math.prod(overlaps[i] for i in accessible),
            math.prod(x for i,x in enumerate(overlaps) if i not in accessible))


def fragment_unitary(c,overlaps,rest_systems):
    """Wires E_0,...,E_(n-1),M,rest. Actual round-100 gates."""
    overlaps=tuple(overlaps)
    partition_overlaps(overlaps,())
    n=len(overlaps)
    total=np.eye(2**(n+rest_systems),dtype=complex)
    for i,lam in enumerate(overlaps):
        total=embed_operator(dephasing_unitary(c,lam),(i,n),n+rest_systems)@total
    return total


def branching_isometry(c,overlaps,rest_dimension):
    plus=np.ones(1)
    minus=np.ones(1)
    for lam in overlaps:
        plus=np.kron(plus,[1.,0.])
        minus=np.kron(minus,vector(lam))
    axis=np.kron(memory_axis(c),np.eye(rest_dimension//2))
    return (np.kron(plus[:,None],(np.eye(rest_dimension)+axis)/2)
            +np.kron(minus[:,None],(np.eye(rest_dimension)-axis)/2))


def fragment_state(state,c,overlaps):
    n=len(overlaps)
    unitary=fragment_unitary(c,overlaps,int(round(math.log2(len(state)))))
    vacuum=np.zeros((2**n,2**n)); vacuum[0,0]=1
    return unitary@np.kron(vacuum,state)@unitary.conj().T


def accessible_state(state,c,overlaps,accessible):
    partition_overlaps(overlaps,accessible)
    n=len(overlaps)
    rest=int(round(math.log2(len(state))))
    return keep_systems(fragment_state(state,c,overlaps),
                        tuple(accessible)+tuple(range(n,n+rest)),n+rest)


def accessible_certificate(c,d,overlaps,accessible):
    c,d=Fraction(c),Fraction(d)
    valid_overlap(c);valid_overlap(d)
    overlaps=tuple(map(Fraction,overlaps));accessible=tuple(accessible)
    a,b=partition_overlaps(overlaps,accessible)
    w=(1-I.exact(d*d)).sqrt()
    return {"c_exact":str(c),"d_exact":str(d),
        "fragment_overlaps_exact":list(map(str,overlaps)),"accessible_indices":list(accessible),
        "accessible_overlap_exact":str(a),"inaccessible_overlap_exact":str(b),
        "all_environment_discarded_overlap_exact":str(a*b),
        "optimal_accessible_H_distance":interval_fields(w*I.exact(max(c,b))),
        "H_distance_without_environment":interval_fields(w*I.exact(max(c,a*b))),
        "full_old_system_available_to_recovery":False,
        "optimality_scope":"All measurements on accessible fragments and M,N; no old AB or inaccessible fragments",
        "pre_environment_memory_and_reference_recoverable_by_coherent_inverse":b==1}


class PartialEnvironmentAccessTests(unittest.TestCase):
    def test_actual_gates_equal_branching_isometry_with_complex_reference(self):
        rho=random_state(np.random.default_rng(103),8)
        for overlaps in ((),(.7,),(.2,.9),(0.,.6,1.)):
            v=branching_isometry(.8,overlaps,8)
            np.testing.assert_allclose(v.conj().T@v,np.eye(8),atol=1e-15)
            np.testing.assert_allclose(fragment_state(rho,.8,overlaps),v@rho@v.conj().T,atol=2e-15)

    def test_every_partial_trace_is_accessible_isometry_after_inaccessible_dephasing(self):
        rho=random_state(np.random.default_rng(203),4)
        overlaps=(.2,.7,.95)
        for count in range(4):
            for subset in combinations(range(3),count):
                _,b=partition_overlaps(overlaps,subset)
                v=branching_isometry(.8,[overlaps[i] for i in subset],4)
                np.testing.assert_allclose(accessible_state(rho,.8,overlaps,subset),
                    v@dephase_memory(rho,.8,b)@v.conj().T,atol=2e-15)

    def test_H_spectral_bound_is_exact_for_different_subsets_and_endpoints(self):
        for c,d,overlaps in ((.8,.8,(.6,.9)),(0.,0.,(0.,.5)),(1.,.4,(.2,1.))):
            for subset in ((),(0,),(1,),(1,0)):
                pair=[accessible_state(memory_channel(calibration_state(H,s),c,d),c,overlaps,subset)
                      for s in (-1,1)]
                _,b=partition_overlaps(overlaps,subset)
                self.assertAlmostEqual(trace_distance(*pair),math.sqrt(1-d*d)*max(c,b),places=13)

    def test_coherent_inverse_attains_remaining_channel_and_full_access_recovers_every_state(self):
        rho=random_state(np.random.default_rng(303),8)
        overlaps=(.3,.8,.9);subset=(2,0)
        accessible=accessible_state(rho,.8,overlaps,subset)
        u=fragment_unitary(.8,[overlaps[i] for i in subset],3)
        vacuum=np.diag([1.,0.,0.,0.])
        np.testing.assert_allclose(u.conj().T@accessible@u,
            np.kron(vacuum,dephase_memory(rho,.8,.8)),atol=2e-15)
        u=fragment_unitary(.8,overlaps,3)
        np.testing.assert_allclose(u.conj().T@fragment_state(rho,.8,overlaps)@u,
            np.kron(np.diag([1.]+[0.]*7),rho),atol=2e-15)

    def test_one_inaccessible_perfect_record_blocks_H_gain(self):
        overlaps=(.9,0.,.95)
        for subset in ((),(0,),(2,),(0,2)):
            row=accessible_certificate(Fraction(4,5),Fraction(4,5),overlaps,subset)
            self.assertEqual(Fraction(row["inaccessible_overlap_exact"]),0)
            self.assertAlmostEqual(row["optimal_accessible_H_distance"]["diagnostic"],.48)
        row=accessible_certificate(Fraction(4,5),Fraction(4,5),overlaps,(0,1,2))
        self.assertAlmostEqual(row["optimal_accessible_H_distance"]["diagnostic"],.6)

    def test_G_information_and_old_unconditional_reference_are_preserved(self):
        pair=[accessible_state(memory_channel(calibration_state(G,s),.8,.6),.8,(.3,.7),(1,))
              for s in (-1,1)]
        self.assertAlmostEqual(trace_distance(*pair),.6,places=13)
        rho=random_state(np.random.default_rng(403),8)
        written=written_state(rho,.8,.6)
        reduced=accessible_state(written,.8,(.3,.7),())
        np.testing.assert_allclose(forget_memories(reduced),old_channel(rho,.8,.6),atol=1e-15)

    def test_exact_equal_strength_access_threshold_at_96_of_100(self):
        c=Fraction(4,5);lam=Fraction(19,20)
        self.assertLess(lam**5,c)
        self.assertGreater(lam**4,c)
        previous=accessible_certificate(c,c,[lam]*100,range(95))
        nextrow=accessible_certificate(c,c,[lam]*100,range(96))
        self.assertAlmostEqual(previous["optimal_accessible_H_distance"]["diagnostic"],.48)
        self.assertGreater(nextrow["optimal_accessible_H_distance"]["diagnostic"],.48)

    def test_invalid_access_and_physical_parameter_domains_are_rejected(self):
        for indices in ((0,0),(-1,),(2,)):
            with self.assertRaises(ValueError): partition_overlaps((.8,.9),indices)
        with self.assertRaises(ValueError): partition_overlaps((1.1,),())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(PartialEnvironmentAccessTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    f=Fraction
    report={"round":103,"model":"Independent pure fragments record the same memory axis; full whole retained",
        "universal_accessible_H_distance":"sqrt(1-d^2)*max(c,product_inaccessible_overlaps)",
        "proof":"Partial trace equals accessible isometry applied to residual dephasing; trace norm isometry invariance",
        "examples":[accessible_certificate(f(4,5),f(4,5),[f(19,20)]*100,range(k))
                    for k in (0,95,96,99,100)],
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("partial_environment_access_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
