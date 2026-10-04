"""Round 446: whole-subject exchange, conserved tagged data, and actual readout.

Reuse the frozen 445 address model and its analytic certificate; no new controls
or separately conditioned data exchange are added.
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
import reciprocal_exchange_dynamics_audit as previous

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'whole_subject_exchange_audit_results.json'
OBS = {}


def data_perm(word, d):
    basis = tuple(it.product(range(d), repeat=len(word)))
    index = {x: i for i, x in enumerate(basis)}
    return np.array([index[tuple(x[j] for j in word)] for x in basis], dtype=int)


def permutation_matrix(rows):
    out = np.zeros((len(rows), len(rows)), dtype=np.int64)
    out[rows, np.arange(len(rows))] = 1
    return out


def site_operator(op, site, n):
    out = np.array([[1]])
    for j in range(n):
        out = np.kron(out, op if site == j else np.eye(op.shape[0], dtype=int))
    return out


@lru_cache(None)
def model():
    words, energy, gamma, z, _, _, k, _, _ = previous.four_model()
    n, d, dd = 4, 2, 16
    index = {w: i for i, w in enumerate(words)}
    transports = [permutation_matrix(data_perm(w, d)) for w in words]
    control_rows = np.concatenate([i*dd+data_perm(w, d) for i, w in enumerate(words)])
    c = permutation_matrix(control_rows)
    full_gamma = np.zeros((len(words)*dd, len(words)*dd), dtype=np.int64)
    for column, word in enumerate(words):
        for pair in it.combinations(range(n), 2):
            nxt = previous.swap(word, pair)
            row = index[nxt]
            identity = tuple(range(n))
            sp = permutation_matrix(data_perm(previous.swap(identity, pair), d))
            full_gamma[row*dd:(row+1)*dd, column*dd:(column+1)*dd] += sp
    full_energy = np.repeat(energy, dd)
    full_z = np.kron(z, np.eye(dd, dtype=int))
    code_rows = np.concatenate([i*dd+data_perm(w, d)
                                for i, w in enumerate(previous.matching_words(n))])
    code_c = permutation_matrix(code_rows)
    full_k = code_c@np.kron(np.array(k, dtype=float), np.eye(dd))@code_c.T
    return (words, energy, gamma, z, c, full_gamma, full_energy,
            full_z, code_c, full_k, transports)


def marginal_pure(psi, site, n=4, reference=2):
    shaped = psi.reshape((2,)*n+(reference,))
    other = [i for i in range(n) if i != site]
    a = shaped.transpose([site, n]+other).reshape(2*reference, -1)
    return a@a.conj().T


def probability_matrix(u_column, initial_word, words):
    q = np.zeros((4, 4))
    for prob, word in zip(np.abs(u_column)**2, words):
        for b in range(4):
            a = initial_word.index(word[b])
            q[b, a] += prob
    return q


class Audit(unittest.TestCase):
    def test_01_exact_whole_subject_conjugacy(self):
        counts = []
        for n, d in ((2, 3), (4, 2), (4, 3)):
            words = tuple(it.permutations(range(n)))
            basis = tuple(it.product(range(d), repeat=n))
            checked = 0
            for word in words:
                for pair in it.combinations(range(n), 2):
                    nxt = previous.swap(word, pair)
                    for x in basis:
                        lhs = previous.swap(tuple(x[j] for j in word), pair)
                        rhs = tuple(x[j] for j in nxt)
                        self.assertEqual(lhs, rhs)
                        checked += 1
            counts.append(dict(subjects=n, data_dimension=d, basis_action_checks=checked))
        _, _, gamma, _, c, whole, energy, _, _, _, _ = model()
        np.testing.assert_array_equal(c.T@whole@c, np.kron(gamma, np.eye(16, dtype=int)))
        np.testing.assert_array_equal(c.T@(energy[:, None]*c), np.diag(energy))
        np.testing.assert_array_equal(c.T@c, np.eye(384, dtype=int))
        OBS['conjugacy'] = dict(exact_basis_checks=counts, full_joint_dimension=384,
            relation_code_data_dimension=48,
            exact_joint_H_equals_C_address_H_C_dagger=True,
            arbitrary_joint_address_data_reference_preserved=True,
            control_frame_is_mathematics_not_a_preparation_protocol=True)

    def test_02_second_order_data_permutation_and_inherited_error(self):
        words, energy, gamma, z, c, whole, full_energy, full_z, cc, full_k, _ = model()
        # Four times the inverse on the one-hop support is integral (energies 2 or 4).
        r4 = np.array([4//e if e else 0 for e in full_energy], dtype=np.int64)
        exact_minus_four_k = full_z.T@whole@(r4[:, None]*(whole@full_z))
        np.testing.assert_array_equal(exact_minus_four_k, -4*full_k)
        matching_words = previous.matching_words(4)
        blocks = []
        for j, initial in enumerate(matching_words):
            for i, final in enumerate(matching_words):
                if i == j:
                    continue
                ti = permutation_matrix(data_perm(initial, 2))
                tf = permutation_matrix(data_perm(final, 2))
                block = full_k[16*i:16*(i+1), 16*j:16*(j+1)]
                np.testing.assert_array_equal(2*block, -tf@ti.T)
                mapping = [final.index(initial[a]) for a in range(4)]
                self.assertEqual(sum(a != b for a, b in enumerate(mapping)), 4)
                blocks.append(dict(initial=list(initial), final=list(final),
                                   physical_source_to_destination=mapping))
        epsilon, tau = F(1, 128), F(2)
        uaddr = previous.unitary(np.diag(energy)+float(epsilon)*gamma, float(tau/epsilon**2))
        framed = c@np.kron(uaddr, np.eye(16))@c.T
        direct = previous.unitary(np.diag(full_energy)+float(epsilon)*whole, float(tau/epsilon**2))
        direct_error = previous.opnorm((direct-framed)@full_z)
        self.assertLess(direct_error, 3e-9)
        effective = previous.unitary(full_k, float(tau))
        joint_error = previous.opnorm(framed@full_z-full_z@effective)
        addr_error = previous.opnorm(uaddr@z-z@previous.unitary(
            np.array(previous.four_model()[6], dtype=float), float(tau)))
        self.assertAlmostEqual(joint_error, addr_error, places=10)
        self.assertLess(joint_error, float(previous.bound(epsilon, tau)))
        OBS['effective_joint_model'] = dict(exact_second_order_blocks=blocks,
            direct_384_dimension_exponential_error=previous.short(direct_error),
            joint_all_input_error=previous.short(joint_error),
            inherited_445_error_bound=previous.short(previous.bound(epsilon, tau)),
            no_new_perturbation_fit_or_parameter_optimization=True,
            arbitrary_joint_code_data_reference_error_bound_inherited_exactly=True)

    def test_03_subject_readouts_and_conserved_tagged_data(self):
        words, _, _, _, c, whole, energy, _, _, _, _ = model()
        x = np.array([[0, 1], [1, 0]])
        z = np.diag([1, -1])
        counts = 0
        for op in (x, z):
            sites = [site_operator(op, b, 4) for b in range(4)]
            for b in range(4):
                physical = np.kron(np.eye(24, dtype=int), sites[b])
                expected = np.zeros((384, 384), dtype=np.int64)
                for j, word in enumerate(words):
                    expected[j*16:(j+1)*16, j*16:(j+1)*16] = sites[word[b]]
                np.testing.assert_array_equal(c.T@physical@c, expected)
                counts += 1
            for tag in range(4):
                tagged = np.zeros((384, 384), dtype=np.int64)
                for j, word in enumerate(words):
                    tagged[j*16:(j+1)*16, j*16:(j+1)*16] = sites[word.index(tag)]
                np.testing.assert_array_equal(whole@tagged, tagged@whole)
                np.testing.assert_array_equal(energy[:, None]*tagged, tagged*energy[None, :])
                np.testing.assert_array_equal(c.T@tagged@c, np.kron(np.eye(24, dtype=int), sites[tag]))
                counts += 1
        # Address-alone exchange has a different physical action on unknown data.
        initial, final = previous.matching_words(4)[:2]
        transport = permutation_matrix(data_perm(final, 2))@permutation_matrix(data_perm(initial, 2)).T
        self.assertFalse(np.array_equal(transport, np.eye(16, dtype=int)))
        OBS['readouts'] = dict(exact_X_Z_generator_checks=counts,
            physical_subject_readout_moves_to_current_address_label=True,
            complete_tagged_data_algebra_commutes_with_H=True,
            physical_subject_data_not_frozen=True,
            changing_representation_without_readout_would_be_wrong=True,
            independent_internal_tagged_data_processing_generated=False)

    def test_04_actual_channel_for_unknown_correlated_data_and_reference(self):
        words, energy, gamma, _, c, _, _, _, _, _, transports = model()
        initial = previous.matching_words(4)[0]
        initial_index = words.index(initial)
        epsilon, tau = F(1, 128), F(2)
        u = previous.unitary(np.diag(energy)+float(epsilon)*gamma, float(tau/epsilon**2))
        q = probability_matrix(u[:, initial_index], initial, words)
        np.testing.assert_allclose(q.sum(axis=0), np.ones(4), atol=1e-12)
        np.testing.assert_allclose(q.sum(axis=1), np.ones(4), atol=1e-12)
        rng = np.random.default_rng(446)
        psi = rng.normal(size=32)+1j*rng.normal(size=32)
        psi /= np.linalg.norm(psi)
        start_t = transports[initial_index]
        exact_output = []
        for j in range(24):
            exact_output.append(u[j, initial_index]*(
                np.kron(transports[j]@start_t.T, np.eye(2))@psi))
        physical_input = np.zeros((384, 2), dtype=complex)
        physical_input[initial_index*16:(initial_index+1)*16] = psi.reshape(16, 2)
        full_u = c@np.kron(u, np.eye(16))@c.T
        np.testing.assert_allclose(full_u@physical_input,
                                  np.array(exact_output).reshape(384, 2), atol=2e-12)
        residuals = []
        for b in range(4):
            actual = sum((marginal_pure(v, b) for v in exact_output),
                         np.zeros((4, 4), dtype=complex))
            predicted = sum((q[b, a]*marginal_pure(psi, a) for a in range(4)),
                            np.zeros((4, 4), dtype=complex))
            residual = previous.opnorm(actual-predicted)
            self.assertLess(residual, 2e-12)
            residuals.append(previous.short(residual))
        OBS['channel'] = dict(initial_matching=list(initial),
            known_matching_arbitrary_correlated_data_and_reference=True,
            exact_receiver_channel='rho_bR(t) = sum_a q_ba(t) rho_aR(0)',
            probability_matrix=[[previous.short(v) for v in row] for row in q],
            receiver_reference_residuals=residuals,
            classical_mixture_claim_for_unknown_initial_graph_coherence=False)

    def test_05_actual_signal_and_old_partner_selection_rule(self):
        words, energy, gamma, _, _, _, _, _, _, _, _ = model()
        initial = previous.matching_words(4)[0]
        a, old_partner, receiver = 0, initial[0], 3
        sine_lower, b, _ = previous.window_certificate()
        lower = F(4, 9)*sine_lower**2-b
        self.assertGreater(lower, F(7, 20))
        epsilon = F(1, 128)
        reports = []
        for tau in (F(19, 10), F(2), F(21, 10)):
            u = previous.unitary(np.diag(energy)+float(epsilon)*gamma, float(tau/epsilon**2))
            q = probability_matrix(u[:, words.index(initial)], initial, words)
            self.assertGreater(q[receiver, a], float(lower))
            leakage = sum(abs(u[j, words.index(initial)])**2
                          for j, w in enumerate(words) if previous.defects(w))
            self.assertLessEqual(q[old_partner, a], leakage+1e-12)
            self.assertLessEqual(leakage, float((6*epsilon)**2))
            for w in previous.matching_words(4):
                self.assertNotEqual(w[old_partner], initial[a])
            reports.append(dict(tau=str(tau),
                actual_cross_component_receiver_signal=previous.short(q[receiver, a]),
                actual_old_partner_signal=previous.short(q[old_partner, a]),
                actual_initial_matching_leakage=previous.short(leakage)))
        # For a single unknown message and blank other data the channel is
        # q*rho_bR + (1-q)*|0><0|_b tensor rho_R, including entangled input.
        q0 = q[receiver, a]
        bell = np.array([1, 0, 0, 1], dtype=complex)/np.sqrt(2)
        rho = np.outer(bell, bell.conj())
        channel = q0*rho+(1-q0)*np.kron(np.diag([1, 0]), np.eye(2)/2)
        self.assertAlmostEqual(float(np.trace(channel).real), 1, places=12)
        self.assertGreaterEqual(float(np.linalg.eigvalsh(channel).min()), -1e-12)
        # Two actual preparations |0000> and |1000> yield receiver trace distance q.
        difference = np.diag([-q0, q0])
        self.assertAlmostEqual(float(np.sum(np.abs(np.linalg.eigvalsh(difference)))/2), q0)
        OBS['signal'] = dict(sender=a, receiver=receiver, old_partner=old_partner,
            initial_edges=[[0, 1], [2, 3]],
            strict_cross_component_signal_window_lower=previous.fraction_data(lower),
            simple_signal_lower='7/20',
            all_time_old_partner_signal_upper=previous.fraction_data((6*epsilon)**2),
            numerical_cross_checks=reports,
            actual_physical_qubit_readout_used=True,
            actual_signal_distinguished_from_graph_change=True,
            universal_message_channel_is_identity=False,
            known_initial_matching_preparation_and_blank_data_explicit=True)


def run():
    OBS.clear()
    buffer = io.StringIO()
    result = unittest.TextTestRunner(stream=buffer, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(buffer.getvalue())
    return dict(round=446, baseline_round=445, date='2026-09-24',
        runtime=dict(python=platform.python_version(), numpy=np.__version__),
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        observations=OBS.copy(),
        scope=dict(whole_subject_swap_includes_address_and_data=True,
            exact_all_N_equal_data_dimension_controlled_permutation_conjugacy=True,
            unknown_joint_code_data_reference_error_inherited_from_445=True,
            fixed_subject_readout_transformed_explicitly=True,
            complete_address_tagged_data_algebra_conserved=True,
            actual_cross_component_data_signal_strictly_certified=True,
            old_partner_legal_code_transport_forbidden_in_this_model=True,
            separate_edge_conditioned_data_swap_added=False,
            actual_data_signal_identified_with_address_motion_only=False,
            arbitrary_unknown_graph_has_simple_classical_routing_channel=False,
            full_unknown_message_perfectly_transferred=False,
            consistency_energy_and_all_pair_access_derived_from_429=False,
            old_relation_locality_or_three_dimensional_space_generated=False,
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
