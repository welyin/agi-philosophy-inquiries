"""Round 91: an exactly compiled weak write and its information-disturbance curve."""

import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from bell_network_statistics import bell_basis
from bilocal_record_tomography import embed_operator
from certified_intervals import SCALE
from complex_whole_real_interfaces import flagged_complex_whole, random_state
from one_bit_real_network import visibility_interval
from operational_effect_closure import compiled_copy_gate
from persistent_relation_memory import sector_projector, record_tree
from quantum_interface_audit import ALPHA, IDENTITY, PAULI_X, PAULI_Z, rotation_unitary, effective_readout
from role_symmetry_and_swap import pauli_word


def weak_write_unitary(angle):
    if not 0 <= angle <= math.pi:
        raise ValueError("Use an angle in [0,pi].")
    decoder=np.kron(IDENTITY,bell_basis())
    return decoder@embed_operator(compiled_copy_gate(angle),(0,2),3)@decoder.conj().T


def weak_write_isometry(angle):
    return weak_write_unitary(angle)@np.kron(np.array([[1.],[0.]]),np.eye(4))


def memory_vectors(angle):
    return np.array([1.,0.]),np.array([math.cos(angle/2),math.sin(angle/2)])


def memory_state(sign,angle):
    vector=memory_vectors(angle)[int(sign==1)]
    return np.outer(vector,vector)


def write_state(state,angle):
    spectator=len(state)//4
    isometry=np.kron(weak_write_isometry(angle),np.eye(spectator))
    return isometry@state@isometry.conj().T


def forget_front_pointer(state):
    d=len(state)//2
    return np.einsum("abad->bd",state.reshape(2,d,2,d))


def coherence_channel(state,c):
    if not -1 <= c <= 1:
        raise ValueError("Use a coherence factor in [-1,1].")
    g=np.kron(pauli_word("YY"),np.eye(len(state)//4))
    return (1+c)*state/2+(1-c)*g@state@g/2


def trace_distance(a,b):
    return float(np.abs(np.linalg.eigvalsh(a-b)).sum()/2)


def read_axis(angle):
    return math.pi/2+angle/2


def tradeoff_certificate(c,v,count):
    c,v=Fraction(c),Fraction(v)
    if not 0 <= c <= 1 or not 0 <= v <= 1 or c*c+v*v!=1:
        raise ValueError("Use nonnegative exact c,v with c^2+v^2=1.")
    success=(1+v*visibility_interval(count))/2
    return {"overlap_exact":str(c),"memory_distinguishability_exact":str(v),
            "worst_old_state_trace_distance_exact":str((1-c)/2),
            "old_readouts":count,
            "success_lower_exact":str(Fraction(success.lo,SCALE)),
            "success_upper_exact":str(Fraction(success.hi,SCALE)),
            "success_diagnostic":sum(success.floats())/2,
            "single_memory_ideal_ceiling_exact":str((1+v)/2)}


class WeakRelationTradeoffTests(unittest.TestCase):
    def test_variable_old_gate_matches_the_analytic_isometry_and_is_reversible(self):
        for angle in (0,.2,.9,2*math.atan2(3,4),math.pi):
            m0,m1=memory_vectors(angle)
            expected=np.kron(m0[:,None],sector_projector(-1))+np.kron(m1[:,None],sector_projector(1))
            np.testing.assert_allclose(weak_write_isometry(angle),expected,atol=5e-16)
            u=weak_write_unitary(angle)
            np.testing.assert_allclose(u.conj().T@u,np.eye(8),atol=9e-16)

    def test_actual_reduced_output_matches_the_channel_on_correlated_complex_inputs(self):
        rng=np.random.default_rng(91)
        for d in (4,8,12):
            state=random_state(rng,d)
            for angle in (.4,1.3,math.pi):
                actual=forget_front_pointer(write_state(state,angle))
                np.testing.assert_allclose(actual,coherence_channel(state,math.cos(angle/2)),atol=2e-16)

    def test_old_promised_whole_remains_unchanged_even_with_nonorthogonal_memory(self):
        for angle in (.2,1.4,math.pi):
            for sign in (-1,1):
                old=flagged_complex_whole(sign)
                np.testing.assert_allclose(write_state(old,angle),
                                           np.kron(memory_state(sign,angle),old),atol=3e-16)

    def test_worst_state_distance_is_attained_by_an_old_real_product_state(self):
        rng=np.random.default_rng(191)
        witness=np.diag([1.,0,0,0])
        for c in (0,.2,.8,1):
            delta=(1-c)/2
            self.assertAlmostEqual(trace_distance(witness,coherence_channel(witness,c)),delta)
            self.assertAlmostEqual(np.trace(coherence_channel(witness,c)@pauli_word("IZ")).real,c)
            for d in (4,8,12):
                state=random_state(rng,d)
                self.assertLessEqual(trace_distance(state,coherence_channel(state,c)),delta+1e-15)

    def test_memory_difference_and_actual_old_readout_match_the_optimal_real_axis(self):
        for angle in (0,.2,1.1,2.7,math.pi):
            c,v=math.cos(angle/2),math.sin(angle/2)
            difference=memory_state(1,angle)-memory_state(-1,angle)
            np.testing.assert_allclose(difference,v*(c*PAULI_X-v*PAULI_Z),atol=3e-16)
            self.assertAlmostEqual(trace_distance(memory_state(1,angle),memory_state(-1,angle)),v)
            table=np.array([[np.trace(memory_state(s,angle)@effective_readout(read_axis(angle),r,1)).real
                             for r in (0,1)] for s in (-1,1)])
            self.assertAlmostEqual(np.max(table,axis=0).sum()/2,(1+ALPHA*v)/2)

    def test_three_actual_protected_reads_give_the_certified_success(self):
        angle=2*math.atan2(3,4)
        turn=rotation_unitary(-read_axis(angle))
        records=[]
        for s in (-1,1):
            state=turn@memory_state(s,angle)@turn.conj().T
            records.append({h:np.trace(branch).real for h,branch in record_tree(state,3).items()})
        actual=sum(max(records[0][h],records[1][h]) for h in records[0])/2
        expected=tradeoff_certificate(Fraction(4,5),Fraction(3,5),3)
        self.assertAlmostEqual(actual,expected["success_diagnostic"],places=14)

    def test_reading_only_memory_cannot_change_the_unconditional_old_reduced_state(self):
        from quantum_interface_audit import branch_kraus
        state=random_state(np.random.default_rng(291),8)
        written=write_state(state,1.2)
        output=np.zeros_like(written)
        for r in (0,1):
            for k in branch_kraus(read_axis(1.2),r,1.):
                full=np.kron(k,np.eye(8))
                output+=full@written@full.conj().T
        np.testing.assert_allclose(forget_front_pointer(output),forget_front_pointer(written),atol=2e-16)

    def test_exact_certificate_exhibits_the_tradeoff_and_nonzero_readout_gap(self):
        for c,v in ((Fraction(4,5),Fraction(3,5)),(Fraction(3,5),Fraction(4,5)),
                    (Fraction(12,13),Fraction(5,13))):
            delta=(1-c)/2
            self.assertEqual(v*v+(1-2*delta)**2,1)
            for count in (1,3,5):
                row=tradeoff_certificate(c,v,count)
                self.assertLess(Fraction(row["success_upper_exact"]),(1+v)/2)
                self.assertGreater(Fraction(row["success_lower_exact"]),Fraction(1,2))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(WeakRelationTradeoffTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":91,
        "write_isometry":"|0> Pi_minus + (cos(theta/2)|0>+sin(theta/2)|1>) Pi_plus",
        "unconditional_old_channel":"(1+c)/2 id +(1-c)/2 Ad_YY, c=cos(theta/2)",
        "worst_old_state_trace_distance_including_spectators":"Delta=(1-c)/2",
        "single_memory_ideal_success":"(1+sqrt(1-c^2))/2",
        "finite_protected_read_success":"(1+gamma_m sqrt(1-c^2))/2",
        "finite_protocol_curve":"((2 P-1)/gamma_m)^2+(1-2 Delta)^2=1",
        "unconditional_reading_only_memory_adds_old_state_disturbance":False,
        "write_pair_interactions":3,"write_local_rotations":3,
        "protected_read_pair_interactions_per_read":1,
        "persistent_memory_rebits":1,"reusable_read_pointer_rebits":1,
        "exact_c_four_fifths_examples":[tradeoff_certificate(Fraction(4,5),Fraction(3,5),m) for m in (1,3,5)],
        "universal_optimality_outside_this_protocol_proved":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("weak_relation_tradeoff_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
