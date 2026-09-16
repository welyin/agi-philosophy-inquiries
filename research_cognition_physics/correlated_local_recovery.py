"""Round 107: prepare the correlated source and recover with local real messages."""
import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from bilocal_record_tomography import embed_operator
from certified_intervals import Interval as I
from complex_whole_real_interfaces import random_state
from correlated_environment_response import correlated_environment,correlated_accessible
from deferred_relation_queries import memory_channel
from environment_query_recovery import interval_abs,interval_max
from fragment_query_compiler import query_branches
from incompatible_relation_writes import H,calibration_state,written_state,forget_memories,old_channel
from independent_source_alignment import keep_systems
from joint_relation_records import memory_axis,axis_read_memory_tree,interval_fields,validated_count
from one_bit_real_network import visibility_interval
from operational_effect_closure import compiled_copy_gate
from quantum_interface_audit import rotation_unitary
from role_symmetry_and_swap import pauli_word,pauli_rotation
from weak_relation_tradeoff import forget_front_pointer,trace_distance


def source_correlator():
    first=embed_operator(pauli_rotation("YX",math.pi/2),(0,1),3)
    before=embed_operator(rotation_unitary(math.pi/2),(1,),3)
    pair=embed_operator(pauli_rotation("XY",math.pi/2),(1,2),3)
    after=before.conj().T
    return after@pair@before@first


def source_preparation(kappa):
    """Six pure slots E0,E1,E2,R0,R1,R2. R purifies E and stays in the whole."""
    if not -1<=kappa<=1: raise ValueError("Use kappa in [-1,1].")
    u=np.eye(64,dtype=complex)
    for i,angle in enumerate((math.acos(kappa),math.pi/2,math.pi/2)):
        u=embed_operator(rotation_unitary(angle),(i,),6)@u
        u=embed_operator(compiled_copy_gate(),(i+3,i),6)@u
    return np.kron(source_correlator(),np.eye(8))@u


def signed_dephase(state,c,t):
    if not -1-1e-14<=t<=1+1e-14: raise ValueError("Signed coherence must have modulus <=1.")
    axis=np.kron(memory_axis(c),np.eye(len(state)//2))
    return (1+t)*state/2+(1-t)*axis@state@axis/2


def direct_eraser_axis(lam0):
    # F=s0 Z-lam0 X maximizes the real part of the conditional response.
    return -math.asin(lam0)


def correlated_messages(state,c,kappa,overlaps,count,undo=False):
    joint=correlated_accessible(state,c,correlated_environment(kappa),overlaps,(0,),undo)
    axis=math.pi/2 if undo else direct_eraser_axis(overlaps[0])
    records={r:np.zeros_like(state,dtype=complex) for r in (-1,1)}
    for history,branch in axis_read_memory_tree(joint,0,axis,count).items():
        r=1 if sum(history)>count//2 else -1
        records[r]+=forget_front_pointer(branch)
    return records


def conditional_coherence(kappa,overlaps,gamma,r,undo=False):
    l0,l1,l2=overlaps
    q=kappa*math.sqrt((1-l1*l1)*(1-l2*l2))
    # After inverse the measured axis is +X, so the correlation term is negative.
    return l1*l2-r*gamma*q if undo else l0*l1*l2+r*gamma*q


def finite_records(old,c,d,kappa,overlaps,m,q,undo=False):
    gamma=sum(visibility_interval(m).floats())/2
    branches=correlated_messages(written_state(old,c,d),c,kappa,overlaps,m,undo)
    output={}
    for r,branch in branches.items():
        t=conditional_coherence(kappa,overlaps,gamma,r,undo)
        action=("joint",1 if t>=0 else -1) if abs(t)>c else ("second_only",1)
        for history,value in query_branches(branch,c,d,q,*action).items():
            output[(r,history)]=value
    return output


def prepared_budget(m,q,undo=False):
    validated_count(m);validated_count(q)
    em=m if m>1 else 0;qm=q if q>1 else 0
    slots=8+int(m>1)+int(q>1)
    return {"scope":"Fixed six-slot purified-source compiler, original MN writes, three couplings, worst joint query branch",
        "pair_rotations":12+int(undo)+em+1+qm,
        "extra_local_rotations":35+5*int(undo)+(3*m+1 if m>1 else 0)+5+3*qm,
        "global_extra_quantum_slots_peak":slots,"initial_pure_preparations":slots,
        "retained_source_purifier_slots":3,"environment_readouts":m,"final_H_readouts":q,
        "total_readouts":m+q,"pointer_resets":m-1+q-1,"return_classical_message_bits":1,
        "coherent_return_E0_M_pair_gates":int(undo),
        "source_requires_role_reversed_XY":True,"transport_or_routing_cost_included":False}


def finite_certificate(c,d,kappa,overlaps,m,q,undo=False):
    c,d,kappa=map(Fraction,(c,d,kappa));overlaps=tuple(map(Fraction,overlaps))
    if len(overlaps)!=3 or any(not 0<=x<=1 for x in (c,d,*overlaps)) or not -1<=kappa<=1:
        raise ValueError("Invalid source or write parameters.")
    l0,l1,l2=overlaps
    base=I.exact(l1*l2 if undo else l0*l1*l2)
    coupling=I.exact(kappa)*((1-I.exact(l1*l1))*(1-I.exact(l2*l2))).sqrt()
    gamma=visibility_interval(m);eta=visibility_interval(q)
    total=I.exact(0)
    for r in (-1,1):
        t=base+(-r if undo else r)*gamma*coupling
        total+=interval_max(I.exact(c),interval_abs(t))/2
    w=(1-I.exact(d*d)).sqrt()
    return {"kappa_exact":str(kappa),"overlaps_exact":list(map(str,overlaps)),
        "inverse_accessible_coupling":undo,
        "finite_H_success":interval_fields((1+eta*w*total)/2),
        "fixed_message_H_factor":interval_fields(total),"resource_budget":prepared_budget(m,q,undo),
        "optimality_scope":"Ideal one-way environment-first real POVM, then arbitrary MN decision; finite construction uses original noisy reads",
        "arbitrary_LOCC_or_complex_environment_POVMs_optimized":False}


class CorrelatedLocalRecoveryTests(unittest.TestCase):
    def test_two_real_pair_gates_encode_Z_as_XYY(self):
        t=source_correlator()
        np.testing.assert_allclose(t@pauli_word("ZII")@t.conj().T,pauli_word("XYY"),atol=1e-15)
        np.testing.assert_allclose(t.imag,0,atol=1e-15)

    def test_six_slot_unitary_prepares_exact_source_and_keeps_purification(self):
        for kappa in (-1.,0.,.4,1.):
            u=source_preparation(kappa);psi=u[:,0]
            np.testing.assert_allclose(u.conj().T@u,np.eye(64),atol=2e-15)
            np.testing.assert_allclose(keep_systems(np.outer(psi,psi.conj()),(0,1,2),6),
                                       correlated_environment(kappa),atol=2e-15)
            self.assertAlmostEqual(np.vdot(psi,psi).real,1,places=13)

    def test_actual_direct_and_inverse_messages_equal_signed_dephasing_on_complex_reference(self):
        state=random_state(np.random.default_rng(107),8)
        for kappa,overlaps in ((0.,(.9,.9,.9)),(1.,(.9,.9,.9)),(1.,(0.,0.,0.)),(-.4,(.6,.7,.8))):
            for undo in (False,True):
                for m in (1,3):
                    gamma=sum(visibility_interval(m).floats())/2
                    for r,branch in correlated_messages(state,.8,kappa,overlaps,m,undo).items():
                        np.testing.assert_allclose(branch,signed_dephase(state,.8,
                            conditional_coherence(kappa,overlaps,gamma,r,undo))/2,atol=2e-15)

    def test_actual_finite_queries_attain_certificates_without_postselection(self):
        f=Fraction
        for lam,m,q,undo in ((f(9,10),1,1,False),(f(9,10),1,1,True),(f(0),3,1,False),(f(0),1,3,False)):
            cert=finite_certificate(f(4,5),f(4,5),1,[lam]*3,m,q,undo)
            success=0.
            for s in (-1,1):
                records=finite_records(calibration_state(H,s),.8,.8,1,[float(lam)]*3,m,q,undo)
                self.assertAlmostEqual(sum(np.trace(rho).real for _,rho in records.values()),1,places=13)
                success+=sum(np.trace(rho).real for guess,rho in records.values() if guess==s)/2
            self.assertAlmostEqual(success,cert["finite_H_success"]["diagnostic"],places=13)

    def test_real_local_axis_optimum_agrees_with_actual_projective_angle_scan(self):
        # Numerical check of the analytic convex endpoint proof, not a global numerical proof.
        c=.8;d=.8;overlaps=(.9,.9,.9)
        pair=[correlated_accessible(memory_channel(calibration_state(H,s),c,d),c,
                correlated_environment(1),overlaps,(0,)) for s in (-1,1)]
        zeta=.9**3;q=1-.9**2
        optimum=.6*(max(c,abs(zeta+q))+max(c,abs(zeta-q)))/2
        for angle in [direct_eraser_axis(.9),*np.linspace(-math.pi,math.pi,19)]:
            axis=rotation_unitary(angle)@pauli_word("Z")@rotation_unitary(angle).conj().T
            value=0.
            for r in (-1,1):
                p=(np.eye(2)+r*axis)/2
                branches=[forget_front_pointer(np.kron(p,np.eye(4))@rho) for rho in pair]
                value+=trace_distance(*branches)
            self.assertLessEqual(value,optimum+2e-15)
        self.assertAlmostEqual(optimum,.5157,places=13)

    def test_all_local_operations_preserve_old_unconditional_complex_reference(self):
        old=random_state(np.random.default_rng(207),8)
        records=finite_records(old,.8,.6,1,(.9,.9,.9),1,1)
        np.testing.assert_allclose(sum(forget_memories(rho) for _,rho in records.values()),
                                   old_channel(old,.8,.6),atol=2e-15)

    def test_zero_overlap_direct_messages_reach_full_ideal_query_with_correlation(self):
        for kappa,expected in ((0.,.8),(1.,1.)):
            ts=[conditional_coherence(kappa,(0.,0.,0.),1,r) for r in (-1,1)]
            self.assertAlmostEqual(sum(max(.8,abs(t)) for t in ts)/2,expected)
        f=Fraction
        result=finite_certificate(f(4,5),f(4,5),1,[f(0)]*3,1,1)
        self.assertGreater(Fraction(result["finite_H_success"]["lower_exact"]),f(793,1000))

    def test_resources_include_source_purifier_and_distinguish_coherent_return(self):
        self.assertEqual([prepared_budget(1,1)[key] for key in
            ("pair_rotations","extra_local_rotations","global_extra_quantum_slots_peak")],[13,40,8])
        self.assertEqual(prepared_budget(1,1,True)["coherent_return_E0_M_pair_gates"],1)
        self.assertEqual(prepared_budget(3,3)["retained_source_purifier_slots"],3)
        self.assertEqual(prepared_budget(3,3)["global_extra_quantum_slots_peak"],10)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(CorrelatedLocalRecoveryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    f=Fraction
    report={"round":107,"direct_real_eraser_axis":"sqrt(1-lambda0^2)*Z-lambda0*X",
        "conditional_coherence":"zeta+r*gamma*kappa*sqrt((1-lambda1^2)*(1-lambda2^2))",
        "finite_examples":[finite_certificate(f(4,5),f(4,5),k,[lam]*3,m,q,undo)
            for lam,k,m,q,undo in ((f(9,10),0,1,1,False),(f(9,10),1,1,1,False),
                (f(9,10),1,1,1,True),(f(0),0,1,1,False),(f(0),1,1,1,False),(f(0),1,3,1,False))],
        "source_preparation_is_explicit":True,"all_purifiers_retained_in_the_whole":True,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("correlated_local_recovery_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
