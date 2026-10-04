"""Round 93: a sharp bound for instruments preserving every state in each sector.

The class allows complex Kraus amplitudes, coarse records, and arbitrary
spectators. It does not cover arbitrary instruments that disturb sector interiors.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE
from complex_whole_real_interfaces import random_state
from one_bit_real_network import visibility_interval
from persistent_relation_memory import sector_projector
from quantum_interface_audit import ALPHA, branch_kraus
from role_symmetry_and_swap import pauli_word
from weak_relation_tradeoff import memory_vectors, read_axis, trace_distance


def qnd_kraus(a,b):
    a,b=np.asarray(a,dtype=complex),np.asarray(b,dtype=complex)
    if a.ndim!=1 or a.shape!=b.shape or len(a)==0:
        raise ValueError("Use two nonempty amplitude vectors of equal size.")
    if not np.isclose(np.vdot(a,a),1) or not np.isclose(np.vdot(b,b),1):
        raise ValueError("Each amplitude vector must have unit norm.")
    return tuple(x*sector_projector(-1)+y*sector_projector(1) for x,y in zip(a,b))


def coherence(a,b):
    return np.dot(a,np.conj(b))


def complex_coherence_channel(state,value):
    spectator=np.eye(len(state)//4)
    p,q=np.kron(sector_projector(-1),spectator),np.kron(sector_projector(1),spectator)
    return p@state@p+q@state@q+value*p@state@q+np.conj(value)*q@state@p


def record_statistics(a,b,groups=None):
    if groups is None: groups=[(i,) for i in range(len(a))]
    if sorted(i for group in groups for i in group)!=list(range(len(a))):
        raise ValueError("Groups must partition all fine outcomes.")
    p=np.array([sum(abs(a[i])**2 for i in group) for group in groups])
    q=np.array([sum(abs(b[i])**2 for i in group) for group in groups])
    return {"record_tv":float(np.abs(p-q).sum()/2),
            "classical_overlap":float(np.sqrt(p*q).sum()),
            "success":float(np.maximum(p,q).sum()/2)}


def actual_noisy_amplitudes(angle):
    m0,m1=memory_vectors(angle)
    a,b,groups=[],[],[]
    for result in (0,1):
        group=[]
        for k in branch_kraus(read_axis(angle),result,1.):
            for row in range(2):
                group.append(len(a))
                a.append((k@m0)[row])
                b.append((k@m1)[row])
        groups.append(tuple(group))
    return np.array(a),np.array(b),groups


def ideal_minimum_disturbance(target):
    target=Fraction(target)
    if not Fraction(1,2)<=target<=1:
        raise ValueError("Use a target success in [1/2,1].")
    g=2*target-1
    return (1-(1-I.exact(g*g)).sqrt())/2


def finite_target_certificate(target,count):
    target=Fraction(target)
    if not Fraction(1,2)<target<1:
        raise ValueError("Use a target success strictly between 1/2 and 1.")
    gamma=visibility_interval(count)
    g=2*target-1
    if gamma.lo<=I.exact(g).hi:
        raise ValueError("This readout budget does not certify the target.")
    value=(1-(1-(g/gamma)**2).sqrt())/2
    ideal=ideal_minimum_disturbance(target)
    return {"target_success_exact":str(target),"old_protected_readouts":count,
            "ideal_qnd_minimum_disturbance_lower_exact":str(Fraction(ideal.lo,SCALE)),
            "ideal_qnd_minimum_disturbance_upper_exact":str(Fraction(ideal.hi,SCALE)),
            "finite_family_required_disturbance_lower_exact":str(Fraction(value.lo,SCALE)),
            "finite_family_required_disturbance_upper_exact":str(Fraction(value.hi,SCALE)),
            "finite_family_required_disturbance_diagnostic":sum(value.floats())/2}


class NondemolitionInformationBoundTests(unittest.TestCase):
    def test_general_complex_amplitudes_define_a_complete_cp_instrument(self):
        rng=np.random.default_rng(93)
        for size in (2,4,7):
            a=rng.normal(size=size)+1j*rng.normal(size=size)
            b=rng.normal(size=size)+1j*rng.normal(size=size)
            a/=np.linalg.norm(a); b/=np.linalg.norm(b)
            ks=qnd_kraus(a,b)
            np.testing.assert_allclose(sum(k.conj().T@k for k in ks),np.eye(4),atol=4e-16)

    def test_every_sector_interior_state_and_its_correlated_spectator_are_preserved(self):
        a=np.array([1,2j,3.])/np.sqrt(14)
        b=np.array([2j,-1,1.])/np.sqrt(6)
        ks=qnd_kraus(a,b)
        for sign in (-1,1):
            p=np.kron(sector_projector(sign),np.eye(2))
            state=p@random_state(np.random.default_rng(193+sign),8)@p
            state/=np.trace(state).real
            output=sum(np.kron(k,np.eye(2))@state@np.kron(k,np.eye(2)).conj().T for k in ks)
            np.testing.assert_allclose(output,state,atol=2e-16)

    def test_complex_coherence_formula_and_exact_worst_disturbance(self):
        rng=np.random.default_rng(293)
        witness=np.diag([1.,0,0,0])
        for _ in range(12):
            a=rng.normal(size=4)+1j*rng.normal(size=4)
            b=rng.normal(size=4)+1j*rng.normal(size=4)
            a/=np.linalg.norm(a); b/=np.linalg.norm(b)
            c=coherence(a,b)
            ks=qnd_kraus(a,b)
            state=random_state(rng,4)
            actual=sum(k@state@k.conj().T for k in ks)
            np.testing.assert_allclose(actual,complex_coherence_channel(state,c),atol=2e-16)
            delta=abs(1-c)/2
            self.assertAlmostEqual(trace_distance(witness,complex_coherence_channel(witness,c)),delta)
            spectator_state=random_state(rng,8)
            self.assertLessEqual(trace_distance(spectator_state,complex_coherence_channel(spectator_state,c)),delta+1e-15)

    def test_classical_records_obey_both_cauchy_schwarz_bounds_even_when_coarsened(self):
        rng=np.random.default_rng(393)
        for _ in range(30):
            a=rng.normal(size=6)+1j*rng.normal(size=6)
            b=rng.normal(size=6)+1j*rng.normal(size=6)
            a/=np.linalg.norm(a); b/=np.linalg.norm(b)
            c=coherence(a,b)
            for groups in (None,[(0,1),(2,3,4),(5,)]):
                row=record_statistics(a,b,groups)
                self.assertLessEqual(abs(c),row["classical_overlap"]+1e-15)
                self.assertLessEqual(row["record_tv"]**2+row["classical_overlap"]**2,1+1e-15)
                self.assertLessEqual(row["record_tv"]**2+(1-abs(1-c))**2,1+1e-15)

    def test_real_two_outcome_instruments_attain_the_ideal_budget_bound(self):
        for delta in (0,.01,.1,.25,.4,.5):
            g=2*math.sqrt(delta*(1-delta))
            a=np.sqrt([(1+g)/2,(1-g)/2])
            b=a[::-1].copy()
            row=record_statistics(a,b)
            self.assertAlmostEqual(abs(1-coherence(a,b))/2,delta,places=14)
            self.assertAlmostEqual(row["success"],.5+math.sqrt(delta*(1-delta)),places=14)

    def test_actual_old_noisy_instrument_keeps_the_same_coherence_with_less_record_information(self):
        angle=2*math.atan2(3,4)
        a,b,groups=actual_noisy_amplitudes(angle)
        self.assertAlmostEqual(coherence(a,b).real,.8)
        self.assertAlmostEqual(coherence(a,b).imag,0)
        row=record_statistics(a,b,groups)
        self.assertAlmostEqual(row["record_tv"],.6*ALPHA)
        self.assertAlmostEqual(row["success"],(1+.6*ALPHA)/2)
        self.assertLess(row["record_tv"]**2+abs(coherence(a,b))**2,1)

    def test_preserving_all_pure_rays_in_each_sector_leaves_exactly_two_scalar_blocks(self):
        values,vectors=np.linalg.eigh(pauli_word("YY").real)
        rays=[]
        for sign in (-1,1):
            basis=vectors[:,np.isclose(values,sign)]
            rays.extend([basis[:,0],basis[:,1],(basis[:,0]+basis[:,1])/np.sqrt(2)])
        constraints=np.vstack([np.kron(v.reshape(1,-1),np.eye(4)-np.outer(v,v)) for v in rays])
        self.assertEqual(16-np.linalg.matrix_rank(constraints),2)
        for p in (sector_projector(-1),sector_projector(1)):
            np.testing.assert_allclose(constraints@p.reshape(-1,order="F"),0,atol=3e-16)

    def test_finite_target_certificates_exceed_the_ideal_lower_bound_and_improve_with_readout_budget(self):
        previous=Fraction(1)
        for count in (1,3,5):
            row=finite_target_certificate(Fraction(99,100),count)
            self.assertGreater(Fraction(row["ideal_qnd_minimum_disturbance_lower_exact"]),Fraction(2,5))
            self.assertGreater(Fraction(row["finite_family_required_disturbance_lower_exact"]),
                               Fraction(row["ideal_qnd_minimum_disturbance_upper_exact"]))
            self.assertLess(Fraction(row["finite_family_required_disturbance_upper_exact"]),previous)
            previous=Fraction(row["finite_family_required_disturbance_lower_exact"])
        value=ideal_minimum_disturbance(Fraction(4,5))
        self.assertLessEqual(Fraction(value.lo,SCALE),Fraction(1,10))
        self.assertGreaterEqual(Fraction(value.hi,SCALE),Fraction(1,10))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(NondemolitionInformationBoundTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":93,
        "scope":"Finite complex CPTP instruments whose sum fixes every state inside either YY sector",
        "complete_kraus_class":"K_rj=a_rj Pi_minus+b_rj Pi_plus",
        "coherence":"C=sum_rj a_rj conjugate(b_rj)",
        "worst_old_state_trace_distance_including_spectators":"|1-C|/2",
        "record_coherence_bound":"TV(p_minus,p_plus)^2+|C|^2<=1",
        "sharp_budget_success_for_0_le_Delta_le_half":"P<=1/2+sqrt(Delta(1-Delta))",
        "sharp_budget_success_for_Delta_ge_half":"P<=1",
        "ideal_bound_attained_by_real_two_outcome_instruments":True,
        "old_finite_noisy_scheme_approaches_but_does_not_attain_ideal_for_nontrivial_strength":True,
        "ninety_nine_percent_examples":[finite_target_certificate(Fraction(99,100),m) for m in (1,3,5)],
        "optimality_for_all_non_qnd_instruments_claimed":False,
        "optimality_for_all_finite_noisy_compilers_claimed":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("nondemolition_information_bound_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
