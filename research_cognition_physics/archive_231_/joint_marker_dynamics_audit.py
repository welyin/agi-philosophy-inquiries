"""Round 448: common single-marker carrier and state-uniform record protection.
All operators have explicit raw-register extensions. No geometry is derived.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import io
import itertools as it
import json
from pathlib import Path
import platform
import unittest

import numpy as np
import reciprocal_exchange_dynamics_audit as old
import internal_partner_response_audit as previous
import whole_subject_exchange_audit as transport

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_marker_dynamics_audit_results.json'
OBS = {}


def parity(word):
    return sum(word[i] > word[j] for i in range(len(word)) for j in range(i+1, len(word))) % 2


def raw(perm, bits):
    return previous.physical_word(perm, tuple(1+b for b in bits))


def packet_swap(word, i, j):
    out = list(word)
    out[2*i:2*i+2], out[2*j:2*j+2] = word[2*j:2*j+2], word[2*i:2*i+2]
    return tuple(out)


def single_marker(word):
    return all(sum(v != -1 for v in word[2*i:2*i+2]) == 1
               for i in range(len(word)//2))


def record_energy(word):
    marker = []
    for i in range(len(word)//2):
        a, d = word[2*i:2*i+2]
        marker.append(a if a >= 0 and d == -1 else
                      (d if a == -1 and d >= 0 else -1))
    return old.defects(tuple(marker))


@lru_cache(None)
def model():
    n = 4
    perms = tuple(it.permutations(range(n)))
    bits = tuple(it.product(range(2), repeat=n))
    words = tuple(raw(p, b) for p in perms for b in bits)
    index = {w: i for i, w in enumerate(words)}
    internal = np.zeros((384, 384), dtype=np.int64)
    packet = np.zeros_like(internal)
    for col, w in enumerate(words):
        for i in range(n):
            internal[index[old.swap(w, (2*i, 2*i+1))], col] += 1
        for i, j in it.combinations(range(n), 2):
            packet[index[packet_swap(w, i, j)], col] += 1
    active = np.array([previous.raw_energy(w) for w in words], dtype=np.int64)
    record = np.array([record_energy(w) for w in words], dtype=np.int64)
    gamma = old.matrices(perms)[1]
    c = transport.permutation_matrix(np.concatenate([
        k*16+transport.data_perm(p, 2) for k, p in enumerate(perms)]))
    match = np.array([old.defects(p) == 0 for p in perms for b in bits])
    z = np.eye(384, dtype=int)[:, match]
    odd = np.array([parity(p) == 1 for p in perms for b in bits])
    all_one = np.array([b == (1, 1, 1, 1) for p in perms for b in bits])
    return perms, bits, words, active, record, internal, packet, gamma, c, z, odd, all_one


class Audit(unittest.TestCase):
    def test_01_raw_common_carrier_and_primitive_closure(self):
        perms, bits, words, active, record, internal, packet, gamma, c, z, odd, all_one = model()
        lookup = {w: i for i, w in enumerate(words)}
        for p in perms:
            for b in bits:
                w = raw(p, b)
                col = lookup[w]
                expected_e = 4-2*sum(p[i] == j and p[j] == i and b[i] == b[j] == 0
                                    for i, j in it.combinations(range(4), 2))
                self.assertEqual(active[col], expected_e)
                self.assertEqual(record[col], old.defects(p))
        bad, total = 0, 0
        for w in words:
            for i, j in it.combinations(range(4), 2):
                self.assertIn(packet_swap(w, i, j), lookup)
                for left, right in it.product((0, 1), repeat=2):
                    total += 1
                    bad += not single_marker(old.swap(w, (2*i+left, 2*j+right)))
        self.assertEqual(bad*2, total)
        np.testing.assert_array_equal(internal.T, internal)
        np.testing.assert_array_equal(packet.T, packet)
        self.assertEqual(z.shape, (384, 48))
        self.assertGreater(np.linalg.norm(packet@z-z@(z.T@packet@z)), 1)
        OBS['carrier'] = dict(single_marker_permutation_dimension=384,
            reciprocal_matching_logic_dimension=48,
            individual_cross_register_swaps_checked=total,
            cross_register_swaps_leaving_single_marker_code=bad,
            complete_packet_and_internal_swaps_exactly_closed=True,
            arbitrary_joint_marker_logic_reference_retained=True,
            reciprocal_matching_subspace_not_invariant=True,
            single_marker_closure_does_not_imply_reciprocal_relation_closure=True)

    def test_02_exact_controlled_frame_and_physical_readout(self):
        perms, bits, words, active, record, internal, packet, gamma, c, z, odd, all_one = model()
        x = np.array([[0, 1], [1, 0]], dtype=int)
        xsum = sum((transport.site_operator(x, i, 4) for i in range(4)),
                   np.zeros((16, 16), dtype=int))
        np.testing.assert_array_equal(c.T@packet@c, np.kron(gamma, np.eye(16, dtype=int)))
        np.testing.assert_array_equal(c.T@internal@c, np.kron(np.eye(24, dtype=int), xsum))
        # If p(i)=j and p(j)=i, the unordered pair of bit indices is unchanged.
        np.testing.assert_array_equal(c.T@(active[:, None]*c), np.diag(active))
        np.testing.assert_array_equal(c.T@c, np.eye(384, dtype=int))
        for site in range(4):
            actual = np.diag([int(w[2*site] == -1) for w in words])
            expected = np.diag([b[p[site]] for p in perms for b in bits])
            np.testing.assert_array_equal(c.T@actual@c, expected)
        # Turning off packet exchange recovers 447 exactly, not merely its spectrum.
        p447, loc447, w447, e447, v447, _ = previous.build(4)
        sel = [j for j, q in enumerate(loc447) if all(a in (1, 2) for a in q)]
        old_indices = [81*k+j for k in range(3) for j in sel]
        old_restriction = (np.diag(e447)+v447)[np.ix_(old_indices, old_indices)]
        np.testing.assert_array_equal(z.T@(np.diag(active)+internal)@z, old_restriction)
        OBS['frame'] = dict(packet_exchange_removed_from_logic_by_exact_conjugacy=True,
            active_potential_remains_marker_conditioned_logic_interaction=True,
            physical_blank_readout_depends_on_current_marker=True,
            unchanged_447_fixed_matching_model_recovered_at_g_zero=True,
            all_joint_unknown_inputs_covered_by_integer_intertwining=True)

    def test_03_free_exchange_and_uniform_rational_certificate(self):
        perms, bits, words, active, record, internal, packet, gamma, c, z, odd, all_one = model()
        identity = np.eye(24, dtype=np.int64)
        square = gamma@gamma
        fourth = square@square
        np.testing.assert_array_equal(gamma@(square-4*identity)@(square-36*identity),
                                      np.zeros((24, 24), dtype=int))
        np.testing.assert_array_equal(np.diag(square), np.full(24, 6))
        np.testing.assert_array_equal(np.diag(fourth), np.full(24, 120))
        # Weights of |Gamma|=6,2,0 in any basis vector are 1/12,3/4,1/6.
        weight_6 = F(120-4*6, 36*(36-4))
        weight_2 = F(36*6-120, 4*(36-4))
        self.assertEqual((weight_6, weight_2), (F(1, 12), F(3, 4)))
        p_sign = np.diag([(-1)**parity(p) for p in perms])
        np.testing.assert_array_equal(p_sign@gamma, -gamma@p_sign)
        # H_without_internal preserves logical Hamming parity; each internal swap flips it.
        l_sign = np.diag([(-1)**sum(b) for p in perms for b in bits])
        np.testing.assert_array_equal(l_sign@internal, -internal@l_sign)
        np.testing.assert_array_equal(l_sign@packet, packet@l_sign)
        sin3 = lambda x: x-x**3/6
        probability_lower = sin3(F(3, 4))**2/12+3*sin3(F(1, 4))**2/4
        self.assertEqual(probability_lower, F(8297, 98304))
        self.assertGreater(probability_lower, F(29, 100)**2)
        compressed_dyson_bound = F(1, 8)/(1-F(1, 48))
        self.assertEqual(compressed_dyson_bound, F(6, 47))
        final_bound = (F(29, 100)-compressed_dyson_bound)**2
        self.assertEqual(final_bound, F(582169, 22090000))
        self.assertGreater(final_bound, F(1, 40))
        OBS['uniform_certificate'] = dict(subjects=4, lambda_value=1, g=1, time='1/8',
            input='any specified matching basis state and logic |1111>',
            penalty_strength_scope='every Delta>0',
            free_odd_probability='sin^2(6*t)/12 + 3*sin^2(2*t)/4',
            sine_polynomial_probability_lower=old.fraction_data(probability_lower),
            compressed_even_Dyson_amplitude_bound=old.fraction_data(compressed_dyson_bound),
            nonreciprocal_and_all_one_probability_lower=old.fraction_data(final_bound),
            simple_strict_lower='1/40',
            arbitrary_unknown_matching_superposition_claimed=False,
            bound_inferred_from_large_Delta_scan=False)

    def test_04_complete_nonzero_internal_dynamics_counterexample(self):
        perms, bits, words, active, record, internal, packet, gamma, c, z, odd, all_one = model()
        initial_perm = old.matching_words(4)[0]
        col = perms.index(initial_perm)*16+15
        selected = odd & all_one
        matching = np.array([old.defects(p) == 0 for p in perms for b in bits])
        bound = F(582169, 22090000)
        reports = []
        for delta in (1, 64):
            h = np.diag(delta*active)+internal+packet
            u = old.unitary(h, 1/8)
            p_selected = float(np.sum(np.abs(u[selected, col])**2))
            p_nonmatching = float(np.sum(np.abs(u[~matching, col])**2))
            self.assertGreater(p_selected, float(bound))
            self.assertGreaterEqual(p_nonmatching, p_selected)
            # Direct compressed Dyson comparison with exact no-internal evolution.
            u0 = old.unitary(np.diag(delta*active)+packet, 1/8)
            compression_error = float(np.linalg.norm((u-u0)[all_one, col]))
            self.assertLess(compression_error, float(F(6, 47)))
            reports.append(dict(Delta=delta,
                odd_marker_and_all_one_probability=old.short(p_selected),
                total_nonreciprocal_marker_probability=old.short(p_nonmatching),
                compressed_amplitude_error=old.short(compression_error)))
        OBS['active_potential_counterexample'] = dict(numerical_cross_checks=reports,
            internal_and_packet_exchange_both_nonzero=True,
            unknown_input_uniform_matching_protection_not_inherited_from_445=True,
            no_information_discarded_or_postselected=True,
            temporary_nonreciprocal_relations_not_declared_forbidden_axiom=True)

    def test_05_complete_record_penalty_raw_extension(self):
        for n in (2, 3, 4):
            alphabet = tuple(range(n))+(-1,)
            local = tuple(it.product(alphabet, repeat=2))
            projectors = []
            for j in range(n):
                projectors.append(np.diag([int(pair in ((j, -1), (-1, j))) for pair in local]))
            for i, j in it.combinations(range(n), 2):
                np.testing.assert_array_equal(projectors[i]@projectors[j], np.zeros_like(projectors[i]))
            self.assertTrue(np.all(np.diag(sum(projectors)) <= 1))
        words = tuple(it.product((0, 1, -1), repeat=4))
        raw_values = [record_energy(w) for w in words]
        self.assertEqual(min(raw_values), 0)
        self.assertEqual(sorted(set(raw_values)), [0, 2])
        self.assertEqual(sum(e == 0 for e in raw_values), 4)
        _, _, _, active, record, internal, packet, gamma, c, z, odd, all_one = model()
        self.assertTrue(np.all(record >= 0))
        self.assertEqual(np.count_nonzero(record == 0), 48)
        self.assertEqual(np.count_nonzero(active == 0), 3)
        np.testing.assert_array_equal(record[:, None]*internal, internal*record[None, :])
        OBS['record_penalty'] = dict(local_projector='|j,blank><j,blank| + |blank,j><blank,j|',
            raw_two_subject_dimension=81, raw_zero_energy_dimension=4, raw_gap=2,
            single_marker_N4_zero_energy_dimension=48,
            active_penalty_N4_zero_energy_dimension=3,
            nonnegative_raw_extension_without_global_code_projector=True,
            complete_marker_readout_instead_of_A_only_is_a_changed_interaction=True)

    def test_06_record_replacement_exact_decoupling_and_protection(self):
        perms, bits, words, active, record, internal, packet, gamma, c, z, odd, all_one = model()
        delta, time = 64, 1/8
        x = np.array([[0, 1], [1, 0]], dtype=int)
        xsum = sum((transport.site_operator(x, i, 4) for i in range(4)),
                   np.zeros((16, 16), dtype=int))
        marker_h = np.diag([delta*old.defects(p) for p in perms])+gamma
        target_h = np.kron(marker_h, np.eye(16, dtype=int))+np.kron(np.eye(24, dtype=int), xsum)
        h = np.diag(delta*record)+internal+packet
        np.testing.assert_array_equal(c.T@h@c, target_h)
        u = old.unitary(h, time)
        factorized = c@np.kron(old.unitary(marker_h, time), old.unitary(xsum, time))@c.T
        error = old.opnorm(u-factorized)
        self.assertLess(error, 2e-12)
        outside = u@z-z@(z.T@u@z)
        maximum_leakage = old.opnorm(outside)**2
        bound = F(36, delta*delta)
        self.assertLess(maximum_leakage, float(bound))
        OBS['record_replacement'] = dict(Delta=delta, lambda_value=1, g=1,
            exact_marker_and_tagged_logic_tensor_sum=True,
            full_input_factorization_error=old.short(error),
            maximum_matching_code_leakage=old.short(maximum_leakage),
            inherited_all_time_all_input_probability_upper=old.fraction_data(bound),
            arbitrary_reference_covered_by_exact_operator_identities=True,
            conditional_tagged_logic_ZZ_of_447_retained=False,
            uniform_internal_field_required_for_this_factorization=True,
            physical_subject_transport_communication_ruled_out=False,
            adding_both_potentials_claimed_derived=False)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(round=448, baseline_round=447, date='2026-09-24',
        runtime=dict(python=platform.python_version(), numpy=np.__version__),
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        observations=OBS.copy(),
        scope=dict(common_marker_logic_carrier_exactly_closed=True,
            full_unknown_encoded_marker_logic_reference_preserved=True,
            physical_subject_readout_tracked_in_controlled_frame=True,
            active_potential_all_Delta_nonreciprocal_probability_certificate=True,
            complete_record_penalty_raw_positive_extension_supplied=True,
            record_penalty_uniform_field_exact_factorization_proved=True,
            roles_of_record_stability_and_logic_interaction_distinguished=True,
            fixed_reciprocal_code_exactly_preserved_by_packet_exchange=False,
            active_potential_protects_all_unknown_logic_as_Delta_grows=False,
            record_replacement_keeps_447_conditional_logic_ZZ=False,
            physical_transport_signal_ruled_out=False,
            either_potential_derived_from_429=False,
            all_possible_record_logic_mechanisms_refuted=False,
            spatial_dimension_or_full_GR_generated=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    report = run()
    if not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(report, ensure_ascii=False, indent=2))
