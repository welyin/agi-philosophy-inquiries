"""Round 441: operational propagation to a moving, identity-selected receiver.

Uses 440's exact path frame and Nachtergaele--Sims 1004.2086, Eq.2.30-2.32.
The new interface combines a unique-label displacement bound with the
complete receiver/reference channel. A path sector remains an input.
"""
import argparse
from fractions import Fraction
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import simultaneous_edge_exchange_audit as previous
import exchange_relation_audit as core

TARGET = Path(__file__).with_name('moving_subject_propagation_audit_results.json')
OBS = {}


def frame(n, j=1., kappa=1.):
    orders = [(0,)+p+(n-1,) for p in itertools.permutations(range(1, n-1))]
    m, d = len(orders), 2**n
    k = np.zeros((d*m, d*m), dtype=np.result_type(j, kappa, np.int64))
    moves = {}
    for r in range(n-1):
        swap = core.swap(n, r, r+1).real.astype(np.int64)
        k += j*np.kron(swap, np.eye(m, dtype=np.int64))
        if 1 <= r <= n-3:
            a = np.zeros((m, m), dtype=np.int64)
            for col, order in enumerate(orders):
                changed = list(order)
                changed[r], changed[r+1] = changed[r+1], changed[r]
                a[orders.index(tuple(changed)), col] = 1
            moves[r] = a
            k += kappa*np.kron(swap, a)
    v = np.zeros((d*m, d*m), dtype=np.int64)
    for g, order in enumerate(orders):
        for col in range(d):
            row = sum(((col >> (n-1-p)) & 1) << (n-1-label)
                      for p, label in enumerate(order))
            v[row*m+g, col*m+g] = 1
    return k, orders, moves, v


def label_positions(n, orders, label):
    return np.tile([p.index(label) for p in orders], 2**n)


def z_at(n, slot):
    return np.diag([1-2*((j >> (n-1-slot)) & 1) for j in range(2**n)])


def reader_z(n, orders, label):
    m = len(orders)
    out = np.zeros((2**n*m, 2**n*m), dtype=np.int64)
    for g, order in enumerate(orders):
        out += np.kron(z_at(n, order.index(label)), np.diag([int(i == g) for i in range(m)]))
    return out


def poisson_tail(z, ell):
    term = z**ell/math.factorial(ell)
    total = term
    for n in range(ell+1, ell+250):
        term *= z/n
        total += term
        if abs(term) < 1e-18*max(total, 1e-300):
            break
    return total


def receiver_bound(distance, t, j, kappa, radius=None, lam=1.):
    radius = distance//2 if radius is None else radius
    assert 0 <= radius < distance
    ell = distance-radius
    local = 2*poisson_tail(6*(abs(j)+abs(kappa))*abs(t), ell)
    motion = min(1., 2*math.exp(-2*lam*(radius+1)+4*abs(kappa)*math.sinh(lam)*abs(t)))
    return min(1., local+motion), local, motion


def kraus_to_label_receiver(isometry, n, graph_count, label):
    env = 2**(n-1)*graph_count
    out = np.zeros((env, 2, isometry.shape[1]), dtype=complex)
    for data in range(2**n):
        bits = [(data >> (n-1-i)) & 1 for i in range(n)]
        rest = 0
        for i, bit in enumerate(bits):
            if i != label:
                rest = 2*rest+bit
        for g in range(graph_count):
            out[rest*graph_count+g, bits[label], :] = isometry[data*graph_count+g, :]
    return out


def choi(kraus):
    dimension = kraus.shape[2]
    vec = kraus.reshape(len(kraus), 2*dimension)
    return vec.T@vec.conj()/dimension


class Audit(unittest.TestCase):
    def close(self, a, b, tol=1e-11):
        self.assertLess(float(np.max(np.abs(np.asarray(a)-np.asarray(b)))), tol)

    def test_01_exact_local_composite_exchange_completion(self):
        rows = []
        for n in (4, 5):
            k, orders, _, v = frame(n, 2, 1)
            graphs = [previous.path_graph(p) for p in orders]
            original, _ = previous.full_model(n, graphs)
            # full_model has J=1 and mu=1; compare its transformed form separately.
            k_one, _, _, _ = frame(n, 1, 1)
            self.assertTrue(np.array_equal(v.T@original@v,
                k_one+(n-1)*np.eye(len(k_one), dtype=np.int64)))
            # Full slot dimension (2n)^n is handled sparsely, not allocated.
            addresses = []
            for data in range(2**n):
                for order in orders:
                    addresses.append(tuple(2*label+((data >> (n-1-r)) & 1)
                                           for r, label in enumerate(order)))
            index = {a: i for i, a in enumerate(addresses)}
            lifted = np.zeros_like(k)
            for col, address in enumerate(addresses):
                for r in range(n-1):
                    # data-only SWAP, coefficient J=2
                    changed = list(address)
                    x, y = address[r], address[r+1]
                    changed[r], changed[r+1] = (x//2)*2+y % 2, (y//2)*2+x % 2
                    self.assertIn(tuple(changed), index)
                    lifted[index[tuple(changed)], col] += 2
                    # whole label+data SWAP, coefficient kappa=1
                    if 1 <= r <= n-3:
                        changed = list(address)
                        changed[r], changed[r+1] = changed[r+1], changed[r]
                        self.assertIn(tuple(changed), index)
                        lifted[index[tuple(changed)], col] += 1
            self.assertTrue(np.array_equal(lifted, k))
            rows.append(dict(subjects=n, slot_dimension=2*n,
                formal_full_slot_dimension=(2*n)**n, exact_code_dimension=len(k),
                sparse_intertwining_all_columns=len(addresses)))
        OBS['local_tensor_completion'] = dict(rows=rows,
            exact_all_unknown_data_graph_reference_embedding=True,
            composite_SWAP_identity='SWAP_(data+label)=SWAP_data*SWAP_label',
            fixed_slot_chain_and_unique_label_sector_are_inputs=True,
            controlled_encoding_physically_prepared_for_free=False,
            slot_local_dimension_uniform_in_N=False)

    def test_02_uniform_tag_displacement(self):
        n, b, beta = 5, 1, 1
        k, orders, moves, _ = frame(n, 3, 1)
        positions = label_positions(n, orders, b)
        right = np.zeros_like(k)
        for r, move in moves.items():
            op = np.kron(core.swap(n, r, r+1).real.astype(np.int64), move)
            right += op*((positions[:, None] == r+1) & (positions[None, :] == r))
        self.assertTrue(np.all(np.sum(np.abs(right), axis=0) <= 1))
        self.assertTrue(np.all(np.sum(np.abs(right), axis=1) <= 1))
        q = right.T@right
        self.assertTrue(np.array_equal(q@q, q))
        k0 = k-right-right.T
        self.assertFalse(np.any(k0*positions[None, :]-positions[:, None]*k0))
        weights = 4.**(positions-beta)
        tilted = weights[:, None]*k/weights[None, :]
        self.close(tilted, k0+4*right+right.T/4)
        anti = (tilted-tilted.T)/(2j)
        self.assertLessEqual(float(np.linalg.norm(anti, 2)), 15/4+1e-12)
        t, radius, lam = .1, 1, math.log(4)
        u = core.evolve(k, t)
        cols = np.flatnonzero(positions == beta)
        rows = np.flatnonzero(np.abs(positions-beta) > radius)
        actual = float(np.linalg.norm(u[np.ix_(rows, cols)], 2)**2)
        bound = 2*math.exp(-2*lam*(radius+1)+4*math.sinh(lam)*t)
        self.assertLessEqual(actual, bound)
        OBS['tag_motion'] = dict(
            unique_label_right_shift_partial_isometry_exact=True,
            tilted_H_identity_exact_binary_rational=True,
            tag_tail_bound='min(1,2 exp[-2 lambda(R+1)+4 abs(kappa) sinh(lambda) abs(t)])',
            proof_uniform_over_data_J_and_reference=True,
            norm_not_claimed_on_duplicate_label_sectors=True,
            finite_operator_witness=dict(subjects=n, J=3, kappa=1, time=t,
                radius=radius, lambda_value=lam,
                all_input_escape_probability_operator_norm=actual, upper_bound=bound))

    def test_03_bounded_local_receiver_extension(self):
        nlabels, b = 5, 2
        # Two slot data qubits and an untouched reference; label blocks explicit.
        ms = [core.swap(3, slot, 2) for slot in (0, 1)]
        ext = np.zeros((200, 200), complex)
        naive = np.zeros_like(ext)
        for a, c in itertools.product(range(nlabels), repeat=2):
            block = a*nlabels+c
            sl = slice(8*block, 8*(block+1))
            naive[sl, sl] = int(a == b)*ms[0]+int(c == b)*ms[1]
            ext[sl, sl] = int(a == b and c != b)*ms[0]+int(c == b and a != b)*ms[1]
            if (a == b)+(c == b) <= 1:
                self.close(ext[sl, sl], naive[sl, sl])
        self.assertAlmostEqual(float(np.linalg.norm(ext, 2)), 1.)
        self.assertAlmostEqual(float(np.linalg.norm(naive, 2)), 2.)
        OBS['bounded_receiver_extension'] = dict(
            local_two_slot_plus_reference_dimension=200,
            reference_dimension=2, extended_norm=1, naive_duplicate_label_norm=2,
            exact_agreement_on_zero_or_one_receiver_label=True,
            LR_constant_independent_of_receiver_window_and_reference_dimension=True,
            source='Nachtergaele-Sims 1004.2086 Eq.2.30-2.32; edge-chain counting',
            commutator_bound='2 ||A|| ||B|| sum_(n>=ell) (6g abs(t))^n/n!')

    def test_04_reference_complete_channel_witness(self):
        n, label_a, alpha, label_b, beta = 5, 0, 0, 3, 3
        j, kappa, t = .5, 1., .02
        k, orders, _, v = frame(n, j, kappa)
        m = len(orders)
        initial = [g for g, order in enumerate(orders)
                   if order[alpha] == label_a and order[beta] == label_b]
        self.assertEqual(len(initial), 2)
        cols = np.array([data*m+g for data in range(2**n) for g in initial])
        isometry = v@core.evolve(k, t)[:, cols]
        self.close(isometry.conj().T@isometry, np.eye(len(cols)))
        channel = kraus_to_label_receiver(isometry, n, m, label_b)
        m_in = len(cols)
        x_source = np.kron(np.kron(core.PAULI[0], np.eye(2**(n-1))), np.eye(len(initial)))
        gamma = 1/3
        damp0 = np.diag([1., math.sqrt(1-gamma)])
        damp1 = np.array([[0., math.sqrt(gamma)], [0., 0.]])
        damping = [np.kron(np.kron(op, np.eye(2**(n-1))), np.eye(len(initial)))
                   for op in (damp0, damp1)]
        base = choi(channel)
        outcomes = []
        bound, local, motion = receiver_bound(abs(alpha-beta), t, j, kappa)
        for name, enc in [('unitary_X', [x_source]), ('amplitude_damping', damping)]:
            altered = np.concatenate([channel@op for op in enc])
            difference = choi(altered)-base
            delta = float(np.sum(np.abs(np.linalg.eigvalsh((difference+difference.conj().T)/2)))/2)
            self.assertLessEqual(delta, bound)
            self.assertGreater(delta, 1e-12)
            outcomes.append(dict(encoding=name, entangled_reference_trace_distance=delta))
        OBS['receiver_reference_channel'] = dict(
            subjects=n, coherent_initial_graph_orders=len(initial),
            unknown_input_dimension=m_in, reference_dimension=m_in,
            initial_sender_slot=alpha, initial_receiver_slot=beta,
            J=j, kappa=kappa, time=t, witnesses=outcomes,
            proven_all_input_TP_encoding_upper_bound=bound,
            local_signal_contribution=local, receiver_motion_contribution=motion,
            witnesses_are_not_claimed_to_compute_the_diamond_norm=True)

    def test_05_fixed_slot_false_subject_signal(self):
        n, a, b = 4, 1, 2
        k, orders, _, _ = frame(n, 0, 1)
        m = len(orders)
        u = core.evolve(k, math.pi/2)
        rho_difference = np.zeros_like(k, dtype=complex)
        for sign, data in [(1, 0), (-1, 1 << (n-1-a))]:
            ket = np.zeros(len(k), complex)
            ket[data*m] = 1
            out = u@ket
            rho_difference += sign*np.outer(out, out.conj())
        moving_reader = reader_z(n, orders, b)
        fixed_slot = np.kron(z_at(n, b), np.eye(m))
        actual = float(np.trace(rho_difference@moving_reader).real)
        wrong = float(np.trace(rho_difference@fixed_slot).real)
        self.assertAlmostEqual(actual, 0.)
        self.assertAlmostEqual(wrong, 2.)
        OBS['identity_vs_slot'] = dict(J=0, kappa=1, time='pi/2',
            true_receiver_Z_gap=actual, old_slot_Z_gap=wrong,
            all_time_no_label_data_signal_at_J_zero_proved_from_original_frame=True,
            tagged_motion_is_not_itself_communication_between_subject_data=True)

    def test_06_all_size_rational_certificate(self):
        distance, radius, ell = 12, 6, 6
        t, g, kappa = Fraction(1, 100), Fraction(3, 2), Fraction(1)
        z = 6*g*t
        tail = z**ell/Fraction(math.factorial(ell))/(1-z/Fraction(ell+1))
        # lambda=log(4): 4 sinh(lambda)=15/2, exp(x)<=1/(1-x).
        exponent = Fraction(15, 2)*kappa*t
        motion = Fraction(2, 4**(2*(radius+1)))/(1-exponent)
        certificate = 2*tail+motion
        self.assertLess(certificate, Fraction(1, 10**8))
        OBS['all_size_certificate'] = dict(
            initial_distance=distance, receiver_radius=radius, time=str(t),
            g=str(g), kappa=str(kappa), lambda_choice='log(4)',
            local_series_tail_upper=str(tail), receiver_motion_upper=str(motion),
            final_trace_distance_upper=str(certificate),
            decimal_upper=float(certificate), strictly_below='1/100000000',
            arbitrary_finite_path_size_supporting_initial_separation=True,
            large_Hilbert_space_numerically_simulated=False,
            simple_bound='min(1,4 exp[6 e (abs(J)+abs(kappa)) abs(t)-d/2])',
            physical_speed_or_three_dimensional_space_derived=False)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=441, baseline_round=440, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(exact_local_composite_exchange_completion_verified=True,
            moving_unique_label_tail_uniform_over_unknown_data_proved=True,
            all_input_reference_TP_encoding_signal_bound_proved=True,
            dynamic_identity_readout_distinguished_from_fixed_slot=True,
            finite_rational_all_size_certificate_proved=True,
            existing_LR_and_prior_fixed_chain_results_explicitly_reused=True,
            new_dynamical_term_added_beyond_440=False,
            arbitrary_initial_label_delocalization_covered=False,
            initial_path_or_three_dimensional_geometry_generated=False,
            label_preparation_and_encoding_free=False,
            strict_relativistic_light_cone_proved=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert json.loads(TARGET.read_text(encoding='utf-8')) == result
    elif not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
