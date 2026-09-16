"""Round 110: deterministic real preparation and finite recovery with nonzero records."""
import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I
from complex_whole_real_interfaces import random_state
from correlated_environment_response import correlated_accessible
from correlated_local_recovery import signed_dephase
from environment_query_recovery import interval_abs,interval_max
from fragment_query_compiler import query_branches
from incompatible_relation_writes import H,calibration_state,written_state,forget_memories,old_channel
from independent_source_alignment import keep_systems
from joint_relation_records import axis_read_memory_tree,interval_fields,validated_count
from nonzero_fragment_records import biased_environment,correlation_capacity,validate_biases
from one_bit_real_network import visibility_interval
from role_symmetry_and_swap import pauli_rotation,pauli_word
from weak_relation_tradeoff import forget_front_pointer


def native_pulses(word,angle):
    """Compile the needed odd-Y rotations to local Y and pair XY/YX pulses."""
    sites=[i for i,x in enumerate(word) if x!="I"]
    if len(sites)<=2:
        if word.count("Y")!=1: raise ValueError("Expected one Y.")
        if "Z" not in word: return [(word,angle)]
        z=word.index("Z")
        local="".join("Y" if i==z else "I" for i in range(len(word)))
        return [(local,math.pi/2),(word.replace("Z","X"),angle),(local,-math.pi/2)]
    routes={"ZYX":("YIX","XYI"),"ZXY":("YXI","XIY"),"ZZY":("YZI","XIY")}
    if word not in routes: raise ValueError("Unsupported three-body rotation.")
    outer,inner=routes[word]
    return native_pulses(outer,-math.pi/2)+native_pulses(inner,-angle)+native_pulses(outer,math.pi/2)


def mapped_pulses(word,angle,sites,total=6):
    result=[]
    for local,a in native_pulses(word,angle):
        full=["I"]*total
        for i,site in enumerate(sites): full[site]=local[i]
        result.append(("".join(full),a))
    return result


def phase_pulses(kind,sites,total=6):
    families={
        "seed":(("IYX",1),("IXY",1),("ZYX",-1),("ZXY",-1)),
        "twirl":(("IYX",1),("IXY",-1),("ZYX",-1),("ZXY",1)),
        "cz":(("IIY",1),("ZIY",-1),("IZY",-1),("ZZY",1))}
    return [pulse for word,sign in families[kind]
            for pulse in mapped_pulses(word,sign*math.pi/2,sites,total)]


def pulse_unitary(pulses,systems):
    out=np.eye(2**systems,dtype=complex)
    for word,angle in pulses: out=pauli_rotation(word,angle)@out
    return out


def biased_source_pulses(biases,tau):
    validate_biases(biases)
    if not -1<=tau<=1: raise ValueError("Use tau in [-1,1].")
    pulses=[]
    for i,r in enumerate(biases):
        pulses+=mapped_pulses("Y",math.acos(r),(i,))
    pulses+=phase_pulses("seed",(0,1,2))
    for control,other in ((3,1),(4,2)):
        pulses+=mapped_pulses("Y",math.pi/2,(control,))
        pulses+=phase_pulses("twirl",(control,0,other))
    pulses+=mapped_pulses("Y",math.acos(tau),(5,))
    pulses+=phase_pulses("cz",(5,0,1))
    return pulses


def biased_source_preparation(biases,tau):
    # E0,E1,E2,R0,R1,R2. No measurement, filter, or discarded failed branch.
    return pulse_unitary(biased_source_pulses(biases,tau),6)


def biased_environment_messages(state,c,biases,kappa,overlaps,count):
    accessible=correlated_accessible(state,c,biased_environment(biases,kappa),overlaps,(0,),True)
    result={r:np.zeros_like(state,dtype=complex) for r in (-1,1)}
    for history,branch in axis_read_memory_tree(accessible,0,math.pi/2,count).items():
        r=1 if sum(history)>count//2 else -1
        result[r]+=forget_front_pointer(branch)
    return result


def finite_biased_records(old,c,d,biases,kappa,overlaps,m,q):
    gamma=sum(visibility_interval(m).floats())/2
    b=overlaps[1]*overlaps[2]
    k=kappa*math.sqrt((1-overlaps[1]**2)*(1-overlaps[2]**2))
    output={}
    for r,branch in biased_environment_messages(written_state(old,c,d),c,biases,kappa,overlaps,m).items():
        t=b-r*gamma*k
        action=("joint",1 if t>=0 else -1) if abs(t)>c else ("second_only",1)
        for history,value in query_branches(branch,c,d,q,*action).items():
            output[(r,history)]=value
    return output


def biased_budget(m,q):
    validated_count(m);validated_count(q)
    em=m if m>1 else 0;qm=q if q>1 else 0
    slots=8+int(m>1)+int(q>1)
    return {"scope":"Fixed compiler including source, old MN writes, three couplings, E0-M inverse, worst joint query branch",
        "source_pair_rotations":29,"source_local_rotations":15,
        "pair_rotations":38+em+qm,
        "extra_local_rotations":46+(3*m+1 if m>1 else 0)+3*qm,
        "global_extra_quantum_slots_peak":slots,"initial_pure_preparations":slots,
        "source_purifier_slots_retained":3,"environment_readouts":m,"final_H_readouts":q,
        "pointer_resets":m-1+q-1,"return_classical_message_bits":1,
        "coherent_return_E0_M_pair_gates":1,"postselection":False,
        "source_uses_role_reversed_XY":True,"transport_and_routing_included":False}


def finite_biased_certificate(c,d,biases,kappa,overlaps,m,q):
    c,d,kappa=map(Fraction,(c,d,kappa))
    biases=tuple(map(Fraction,biases));overlaps=tuple(map(Fraction,overlaps))
    validate_biases(biases)
    if len(overlaps)!=3 or any(not 0<=x<=1 for x in (c,d,*overlaps)):
        raise ValueError("Invalid parameters.")
    if kappa*kappa>math.prod(1-r*r for r in biases): raise ValueError("Unphysical source.")
    b=I.exact(overlaps[1]*overlaps[2])
    k=I.exact(kappa)*((1-I.exact(overlaps[1]**2))*(1-I.exact(overlaps[2]**2))).sqrt()
    gamma=visibility_interval(m);eta=visibility_interval(q)
    factor=sum((interval_max(I.exact(c),interval_abs(b-r*gamma*k))/2 for r in (-1,1)),I.exact(0))
    baseline=I.exact(max(c,overlaps[1]*overlaps[2]))
    w=(1-I.exact(d*d)).sqrt()
    return {"biases_exact":list(map(str,biases)),"kappa_exact":str(kappa),
        "gate_cosines_exact":list(map(str,overlaps)),
        "finite_H_success":interval_fields((1+eta*w*factor)/2),
        "independent_source_same_recovery_success":interval_fields((1+eta*w*baseline)/2),
        "strict_success_gain":interval_fields(eta*w*(factor-baseline)/2),
        "resource_budget":biased_budget(m,q),"finite_same_resource_global_optimality_proved":False}


class BiasedSourceCompilerTests(unittest.TestCase):
    def test_three_body_conjugation_compilers_match_exact_target_rotations(self):
        for word in ("ZYX","ZXY","ZZY"):
            for angle in (.3,-math.pi/2):
                pulses=native_pulses(word,angle)
                np.testing.assert_allclose(pulse_unitary(pulses,3),pauli_rotation(word,angle),atol=1e-15)

    def test_phase_gadgets_act_exactly_on_all_basis_states_and_spectator(self):
        seed=np.ones(8);seed[[4,7]]=-1
        twirl=np.ones(8);twirl[[5,6]]=-1
        cz=np.ones(8);cz[[6,7]]=-1
        for kind,diagonal in (("seed",seed),("twirl",twirl),("cz",cz)):
            u=pulse_unitary(phase_pulses(kind,(0,1,2),3),3)
            np.testing.assert_allclose(u,np.diag(diagonal),atol=2e-15)

    def test_full_six_slot_preparation_yields_biased_source_without_filtering(self):
        for biases,tau in (((.6,)*3,0.),((.6,)*3,1.),((.2,.4,.8),-.7),((1.,.4,.6),.5)):
            u=biased_source_preparation(biases,tau);psi=u[:,0]
            self.assertAlmostEqual(np.vdot(psi,psi).real,1,places=13)
            actual=keep_systems(np.outer(psi,psi.conj()),(0,1,2),6)
            np.testing.assert_allclose(actual,biased_environment(biases,tau*correlation_capacity(biases)),atol=2e-15)
            np.testing.assert_allclose(u.conj().T@u,np.eye(64),atol=3e-15)

    def test_every_source_pulse_is_an_original_real_local_or_role_symmetric_pair_rotation(self):
        pulses=biased_source_pulses((.6,)*3,.4)
        counts=[sum(ch!="I" for ch in word) for word,angle in pulses]
        self.assertEqual(counts.count(1),15)
        self.assertEqual(counts.count(2),29)
        for word,angle in pulses:
            self.assertEqual(word.count("Y"),1)
            self.assertTrue(set(word)<=set("IXY"))

    def test_actual_noisy_environment_instrument_is_independent_of_bias_at_fixed_kappa(self):
        state=random_state(np.random.default_rng(110),8)
        for biases in ((.6,.4,.2),(.2,.3,.4)):
            for m in (1,3):
                gamma=sum(visibility_interval(m).floats())/2
                for r,branch in biased_environment_messages(state,.8,biases,.3,(.9,)*3,m).items():
                    np.testing.assert_allclose(branch,signed_dephase(state,.8,.81-r*gamma*.3*.19)/2,atol=2e-15)

    def test_actual_finite_H_queries_attain_strict_certified_gain(self):
        f=Fraction
        for m,q in ((1,1),(3,1),(1,3)):
            cert=finite_biased_certificate(f(4,5),f(4,5),[f(3,5)]*3,f(64,125),[f(9,10)]*3,m,q)
            success=0.
            for s in (-1,1):
                records=finite_biased_records(calibration_state(H,s),.8,.8,(.6,)*3,.512,(.9,)*3,m,q)
                self.assertAlmostEqual(sum(np.trace(rho).real for _,rho in records.values()),1,places=13)
                success+=sum(np.trace(rho).real for guess,rho in records.values() if guess==s)/2
            self.assertAlmostEqual(success,cert["finite_H_success"]["diagnostic"],places=13)
            self.assertGreater(Fraction(cert["strict_success_gain"]["lower_exact"]),f(12,1000))

    def test_old_unconditional_complex_reference_channel_is_preserved(self):
        old=random_state(np.random.default_rng(210),8)
        records=finite_biased_records(old,.8,.6,(.6,)*3,.512,(.9,)*3,1,1)
        np.testing.assert_allclose(sum(forget_memories(rho) for _,rho in records.values()),old_channel(old,.8,.6),atol=2e-15)

    def test_budget_keeps_purifiers_and_coherent_return_explicit(self):
        b=biased_budget(1,1)
        self.assertEqual((b["pair_rotations"],b["extra_local_rotations"],b["global_extra_quantum_slots_peak"]),(38,46,8))
        self.assertEqual(b["coherent_return_E0_M_pair_gates"],1)
        self.assertFalse(b["postselection"])
        self.assertEqual(biased_budget(3,3)["global_extra_quantum_slots_peak"],10)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(BiasedSourceCompilerTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    f=Fraction
    report={"round":110,"source_parameterization":"kappa=tau*product_i sqrt(1-r_i^2)",
        "source_preparation":"Pure biased product, real sign gate, two coherent parity records, one coherent phase record",
        "source_pair_rotations":29,"source_local_rotations":15,"source_pure_slots":6,
        "messages_after_actual_E0_M_inverse":"J_r(omega)=Lambda_(B-r*gamma*kappa*s1*s2)(omega)/2",
        "finite_examples":[finite_biased_certificate(f(4,5),f(4,5),[f(3,5)]*3,k,[f(9,10)]*3,m,q)
                          for k,m,q in ((f(0),1,1),(f(64,125),1,1),(f(64,125),3,1),(f(64,125),3,3))],
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("biased_source_compiler_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
