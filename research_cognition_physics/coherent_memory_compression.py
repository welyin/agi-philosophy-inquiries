"""Round 94: six real pulses compress two memories of the same binary relation."""

import argparse
import json
import math
import unittest
from pathlib import Path

import numpy as np

from bilocal_record_tomography import embed_operator
from complex_whole_real_interfaces import flagged_complex_whole, random_state
from persistent_relation_memory import sector_projector
from role_symmetry_and_swap import pauli_rotation
from weak_relation_tradeoff import (weak_write_unitary, weak_write_isometry,
    coherence_channel, trace_distance)


def vector(c):
    if not 0 <= c <= 1:
        raise ValueError("Use a real overlap in [0,1].")
    return np.array([c, math.sqrt((1-c)*(1+c))])


def compressor_pulses(a, b):
    """Chronological Pauli rotations on accumulator A and scratch B.

    XY is an explicitly available role-reversed interaction (round 50).
    This compiler does not claim to work with only a fixed-direction YX gate.
    """
    sa, sb = vector(a)[1], vector(b)[1]
    phi = -2*math.atan2(sb, b)
    beta = math.atan2(a*sb, sa) if sa or a*sb else 0.
    return [("YI", math.pi/2), ("XY", -phi/2),
            ("YI", -math.pi/2), ("IY", phi/2),
            ("XY", -beta), ("YX", beta)]


def compressor(a, b):
    result = np.eye(4, dtype=complex)
    for word, angle in compressor_pulses(a, b):
        result = pauli_rotation(word, angle)@result
    return result


def plane_rotation(i, j, angle):
    """Independent two-coordinate expression, not a primitive control."""
    out = np.eye(4)
    out[i, i] = out[j, j] = math.cos(angle)
    out[i, j], out[j, i] = -math.sin(angle), math.sin(angle)
    return out


def geometric_compressor(a, b):
    sa, sb = vector(a)[1], vector(b)[1]
    beta = math.atan2(a*sb, sa) if sa or a*sb else 0.
    return plane_rotation(1, 2, beta)@plane_rotation(2, 3, -math.atan2(sb, b))


def stream_unitary(overlaps):
    """Four wires A(accumulator), B(scratch), old C,D; scratch is reused coherently."""
    overlaps = tuple(overlaps)
    if not overlaps:
        raise ValueError("At least one write is required.")
    for c in overlaps:
        vector(c)
    total = embed_operator(weak_write_unitary(2*math.acos(overlaps[0])), (0, 2, 3), 4)
    a = overlaps[0]
    for b in overlaps[1:]:
        write = embed_operator(weak_write_unitary(2*math.acos(b)), (1, 2, 3), 4)
        collect = np.kron(compressor(a, b), np.eye(4))
        total = collect@write@total
        a *= b
    return total


def stream_isometry(overlaps):
    return stream_unitary(overlaps)@np.kron(np.eye(4)[:, :1], np.eye(4))


def ideal_stream_isometry(overlaps):
    c = math.prod(overlaps)
    m0 = np.array([1., 0.])
    return (np.kron(np.kron(m0, m0)[:, None], sector_projector(-1))
            + np.kron(np.kron(vector(c), m0)[:, None], sector_projector(1)))


def trace_memories(state):
    d = len(state)//4
    return np.einsum("abad->bd", state.reshape(4, d, 4, d))


class CoherentMemoryCompressionTests(unittest.TestCase):
    def test_six_actual_pulses_equal_two_plane_rotations_including_endpoints(self):
        for a in (0, .1, .8, 1):
            for b in (0, .4, .8, 1):
                np.testing.assert_allclose(compressor(a,b), geometric_compressor(a,b), atol=7e-16)
                self.assertEqual(sum(w.count("I")==0 for w,_ in compressor_pulses(a,b)), 3)

    def test_gate_is_real_proper_and_maps_both_candidates_with_same_pure_garbage(self):
        for a,b in ((0,0),(0,.7),(.2,.8),(.8,.8),(1,.5),(1,1)):
            u = compressor(a,b)
            np.testing.assert_allclose(u.imag, 0, atol=1e-16)
            np.testing.assert_allclose(u.T@u, np.eye(4), atol=8e-16)
            self.assertAlmostEqual(np.linalg.det(u).real, 1)
            np.testing.assert_allclose(u@np.eye(4)[:,0], np.eye(4)[:,0], atol=4e-16)
            np.testing.assert_allclose(u@np.kron(vector(a),vector(b)),
                                       np.kron(vector(a*b),[1.,0]), atol=7e-16)

    def test_reused_scratch_stream_matches_full_operator_isometry(self):
        for cs in ((1,),(.8,),(.8,)*3,(.3,1,.7,0),(.999,)*15):
            np.testing.assert_allclose(stream_isometry(cs), ideal_stream_isometry(cs), atol=8e-15)

    def test_complex_correlated_spectators_and_the_old_channel_are_preserved(self):
        rng = np.random.default_rng(94)
        for d in (4,8,12):
            rho = random_state(rng,d)
            cs = (.8,.7,.9)
            w = np.kron(stream_isometry(cs),np.eye(d//4))
            ideal = np.kron(ideal_stream_isometry(cs),np.eye(d//4))
            actual = w@rho@w.conj().T
            np.testing.assert_allclose(actual, ideal@rho@ideal.conj().T, atol=7e-16)
            np.testing.assert_allclose(trace_memories(actual),coherence_channel(rho,math.prod(cs)), atol=8e-16)

    def test_flagged_old_whole_factorizes_and_remains_identical_in_each_sector(self):
        cs = (.8,)*3
        w = np.kron(stream_isometry(cs),np.eye(2))
        for sign in (-1,1):
            rho = flagged_complex_whole(sign)
            m = np.kron(vector(math.prod(cs)) if sign==1 else [1.,0],[1.,0])
            np.testing.assert_allclose(w@rho@w.conj().T,np.kron(np.outer(m,m),rho), atol=8e-16)

    def test_the_complete_stream_is_reversible_on_unknown_old_inputs(self):
        u = stream_unitary((.8,.4,.9,.2))
        np.testing.assert_allclose(u.conj().T@u,np.eye(16),atol=9e-15)
        rho = random_state(np.random.default_rng(194),4)
        initial = np.kron(np.diag([1.,0,0,0]),rho)
        np.testing.assert_allclose(u.conj().T@(u@initial@u.conj().T)@u,initial,atol=4e-15)

    def test_one_stronger_write_has_exactly_the_same_effective_isometry(self):
        for cs in ((.8,)*3,(.2,.7,.95)):
            # Select scratch=0 only after its deterministic factorization was checked.
            w = stream_isometry(cs).reshape(2,2,4,4)[:,0].reshape(8,4)
            np.testing.assert_allclose(w,weak_write_isometry(2*math.acos(math.prod(cs))),atol=3e-15)

    def test_old_disturbance_bound_is_still_attained_after_compression(self):
        rho = np.diag([1.,0,0,0])
        w = stream_isometry((.8,)*3)
        self.assertAlmostEqual(trace_distance(trace_memories(w@rho@w.conj().T),rho),.244)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(CoherentMemoryCompressionTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":94,"compiler":"3 pair rotations + 3 local rotations per merge",
        "extra_control_assumption":"Explicit XY role reversal of YX, as in round 50",
        "compressor_pulses_a_b_four_fifths":compressor_pulses(.8,.8),
        "overlap_update":"a,b -> a*b","persistent_memory_qubits":1,"reusable_scratch_qubits":1,
        "scratch_returned_pure_and_uncorrelated_for_arbitrary_old_inputs":True,
        "entire_stream_unitary_before_readout":True,
        "n_write_pair_pulses":"6*n-3","n_write_local_pulses":"6*n-3",
        "same_isometry_as_one_stronger_write":True,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("coherent_memory_compression_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
