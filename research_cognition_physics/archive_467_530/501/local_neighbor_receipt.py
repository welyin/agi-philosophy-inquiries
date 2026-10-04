"""Round 501: local returned-label readout in a specified coloured variant.

The new local identity coupling and multilevel data are inputs. Neither a
qubit-only implementation nor automatic retries after root failure is claimed.
"""
import argparse
from fractions import Fraction as F
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import branching_tree_distance_audit as trees_model
import current_neighbor_detection as frozen500

HERE = Path(__file__).resolve().parent
TARGET = HERE/'local_neighbor_receipt_results.json'
OBS = {}


def model(trees, flip, n, root=0):
    """Unit J,kappa,g; port, graph, active colour ordering; query colour 0."""
    labels = [None]+[v for v in range(n) if v != root]
    ng = len(trees)
    base = np.kron(np.eye(n, dtype=np.int64), flip)
    for a, tree in enumerate(trees):
        adj = np.zeros((n, n), dtype=np.int64)
        for u, v in tree:
            adj[u, v] = adj[v, u] = 1
        ids = np.arange(n)*ng+a
        base[np.ix_(ids, ids)] += adj-np.diag(adj.sum(axis=0))
    h = np.kron(base, np.eye(n, dtype=np.int64))
    for v in range(n):
        if v == root:
            continue
        colour = labels.index(v)
        for a in range(ng):
            query = (v*ng+a)*n
            h[query, query+colour] += 1
            h[query+colour, query] += 1
    source = (root*ng+np.arange(ng))*n
    reply = np.zeros(len(h), bool)
    bad = reply.copy()
    for a, tree in enumerate(trees):
        for colour, label in enumerate(labels[1:], 1):
            row = (root*ng+a)*n+colour
            reply[row] = True
            bad[row] = tuple(sorted((root, label))) not in tree
    return h, labels, source, reply, bad


def six(root=0):
    trees, _, _, flip = frozen500.six_model()
    return trees, flip, model(trees, flip, 6, root)


def source_taylor(h, source, denominator=64, degree=12):
    """Integer Taylor columns: no dense object-valued matrix powers."""
    h = h.astype(object)
    power = np.eye(len(h), dtype=object)[:, source]
    real, imag = np.zeros_like(power), np.zeros_like(power)
    scale = denominator**degree*math.factorial(degree)
    for k in range(degree+1):
        coefficient = denominator**(degree-k)*math.factorial(degree)//math.factorial(k)
        (real if k % 2 == 0 else imag)[:] += ((-1)**((k+1)//2))*coefficient*power
        if k < degree:
            power = h@power
    bound = max(max(sum(abs(int(x)) for x in row) for row in h),
                max(sum(abs(int(x)) for x in row) for row in h.T))
    x = F(bound, denominator)
    assert x < degree+2
    tail = x**(degree+1)/(math.factorial(degree+1)*(1-x/(degree+2)))
    return real, imag, scale, tail, bound


def effect_bound(data, effect):
    real, imag, scale, tail, _ = data
    return frozen500.exact_effect_bounds(real, imag, scale,
                                         np.arange(real.shape[1]), effect, tail)


class Checks(unittest.TestCase):
    def test_01_full_tensor_qudit_restriction(self):
        n, local_d, ng = 3, 5, 3
        trees = sorted({trees_model.prufer_tree(n, (v,)) for v in range(n)}, key=lambda g: sorted(g))
        h, labels, _, _, _ = model(trees, np.zeros((ng, ng), dtype=np.int64), n)
        words = list(itertools.product(range(local_d), repeat=n))
        index = {w:k for k,w in enumerate(words)}
        full = np.zeros((len(words)*ng, len(words)*ng), dtype=np.int64)
        for column, word in enumerate(words):
            for a, tree in enumerate(trees):
                for u, v in tree:
                    target = list(word)
                    target[u], target[v] = target[v], target[u]
                    full[index[tuple(target)]*ng+a, column*ng+a] += 1
                for v in range(1, n):
                    if word[v] in (1, 2+v):
                        target = list(word)
                        target[v] = 2+v if word[v] == 1 else 1
                        full[index[tuple(target)]*ng+a, column*ng+a] += 1
        chosen = []
        for v in range(n):
            for a in range(ng):
                for label in labels:
                    word = [0]*n
                    word[v] = 1 if label is None else 2+label
                    chosen.append(index[tuple(word)]*ng+a)
        outside = sorted(set(range(len(full)))-set(chosen))
        self.assertFalse(np.any(full[np.ix_(outside, chosen)]))
        self.assertTrue(np.array_equal(full[np.ix_(chosen, chosen)]-2*np.eye(len(h), dtype=np.int64), h))
        OBS['full_tensor'] = dict(N=3, local_dimension=5, full_dimension=len(full),
            active_dimension=len(h), excluded_root_reply_exactly_dark=True,
            scalar_removed=2, restriction_exact=True, original_qubit_model_unchanged=False)

    def test_02_commutator_and_leading_map(self):
        cases = []
        tr, flip, _ = six()
        four = sorted({trees_model.prufer_tree(4, w)
                       for w in itertools.product(range(4), repeat=2)}, key=lambda g: sorted(g))
        ff = np.zeros((len(four), len(four)), dtype=np.int64)
        for a, tree in enumerate(four):
            for other, weight in trees_model.tree_flips(tree, 4).items():
                ff[four.index(other), a] = weight
        for n, trees, fmat, root in [(6, tr, flip, 0), (6, tr, flip, 1), (4, four, ff, 0)]:
            h, labels, source, reply, bad = model(trees, fmat, n, root)
            ng = len(trees)
            fg = np.kron(np.kron(np.eye(n, dtype=np.int64), fmat), np.eye(n, dtype=np.int64))
            comm = fg@h-h@fg
            degree = np.array([sum(root in e for e in tree) for tree in trees])
            charge = np.tile(np.repeat(degree, n), n)
            self.assertFalse(np.any(h*(charge[:, None]-charge[None, :])))
            root_max = max(np.sum(abs(comm), axis=1))
            self.assertLessEqual(root_max, 48)
            power = np.eye(len(h), dtype=np.int64)[:, source]
            for order in range(4):
                if order < 3:
                    self.assertFalse(np.any(power[reply]))
                else:
                    expected = np.zeros_like(power)
                    for a, tree in enumerate(trees):
                        for colour, label in enumerate(labels[1:], 1):
                            expected[(root*ng+a)*n+colour, a] = int(tuple(sorted((root, label))) in tree)
                    self.assertTrue(np.array_equal(power*reply[:, None], expected))
                    self.assertTrue(np.array_equal(expected.T@expected, np.diag(degree)))
                    self.assertFalse(np.any(power[bad]))
                power = h@power
            cases.append(dict(N=n, graph_states=ng, root=root,
                root_degrees=sorted(set(map(int, degree))), commutator_row_bound=int(root_max),
                third_order_map_exact=True, fourth_order_bad_nonzero=bool(np.any(power[bad]))))
        OBS['operator_crosschecks'] = cases

    def test_03_rational_all_input_root_receipt(self):
        _, _, (h, _, source, reply, bad) = six()
        data = source_taylor(h, source)
        lower, upper = effect_bound(data, reply)
        _, wrong = effect_bound(data, bad)
        self.assertGreater(lower, F(1, 10**13))
        self.assertLess(wrong/lower, F(1, 10000))
        OBS['rational_receipt_certificate'] = dict(time='1/64', Taylor_degree=12,
            h_norm_bound=data[-1], unitary_column_tail=str(data[-2]),
            success_lower=str(lower), success_upper=str(upper), bad_mass_upper=str(wrong),
            joint_bad_ratio_upper=str(wrong/lower),
            readable_success_lower=float(lower), readable_joint_bad_ratio_upper=float(wrong/lower),
            arithmetic='unbounded integer columns, Fraction, complex Gershgorin',
            all_six_graph_inputs_and_passive_reference=True,
            finite_model_time_outside_conservative_uniform_window=True)

    def test_04_complete_local_instrument_and_reference(self):
        trees, _, (h, labels, source, reply, bad) = six()
        u = frozen500.evolution(h, F(1, 64))
        up = u[:, source]
        psi = np.array([[1, 1j], [2j, -1], [1+1j, 2], [-1j, 1], [2, -2j], [1, 3j]], complex)
        psi /= np.linalg.norm(psi)
        evolved = up@psi
        gram = np.zeros((6, 6), complex)
        reference = np.zeros((2, 2), complex)
        outcomes = []
        for colour in range(1, 6):
            mask = np.zeros(len(h), bool)
            mask[np.arange(6)*6+colour] = True
            outcomes.append(mask)
        outcomes.append(~reply)
        self.assertTrue(np.all(np.sum(outcomes, axis=0) == 1))
        for mask in outcomes:
            gram += up[mask].conj().T@up[mask]
            branch = evolved[mask]
            reference += branch.T@branch.conj()
        s = float(np.linalg.norm(evolved[reply])**2)
        wrong = float(np.linalg.norm(evolved[bad])**2)
        self.assertLess(np.linalg.norm(gram-np.eye(6)), 2e-12)
        self.assertLess(np.linalg.norm(reference-psi.T@psi.conj()), 2e-12)
        self.assertGreater(s, 0)
        self.assertLess(wrong/s, 1/100)
        # Read record is retained: graph non-neighbor event after delay refers
        # to the FIXED returned label, not to the later packet's location.
        delayed = frozen500.evolution(h, F(1, 4096))
        wrong_delay = 0.0
        for colour, label in enumerate(labels[1:], 1):
            state = np.zeros_like(evolved)
            mask = outcomes[colour-1]
            state[mask] = evolved[mask]
            later = delayed@state
            missing = np.array([(0, label) not in tree for tree in trees])
            wrong_delay += float(np.linalg.norm(later[np.tile(np.repeat(missing, 6), 6)])**2)
        # Uses the certified finite-model bound rather than the general window.
        finite_eps = F(OBS['rational_receipt_certificate']['joint_bad_ratio_upper'])
        delayed_bound = (math.sqrt(float(finite_eps))+24/4096)**2
        self.assertLess(wrong_delay/s, delayed_bound)
        OBS['local_instrument_diagnostics'] = dict(success=s, abort=1-s,
            joint_wrong_given_success=wrong/s, delayed_joint_wrong_given_success=wrong_delay/s,
            root_only_outcomes=6, completeness_residual=float(np.linalg.norm(gram-np.eye(6))),
            nonselective_reference_residual=float(np.linalg.norm(reference-psi.T@psi.conj())),
            graph_measured=False, other_ports_read=False, selected_reference_unchanged_claim=False)

    def test_05_failure_does_not_reset(self):
        # N=2, J=g=1, kappa=0; ordering root-q, neighbor-q, neighbor-r, root-r.
        h = np.diag(np.ones(3, dtype=np.int64), 1)+np.diag(np.ones(3, dtype=np.int64), -1)
        data = source_taylor(h, np.array([0]), denominator=8, degree=12)
        outside = np.array([False, True, True, False])
        fail = np.array([True, True, True, False])
        out_lower, _ = effect_bound(data, outside)
        fail_lower, _ = effect_bound(data, fail)
        self.assertGreater(out_lower, F(1, 100))
        self.assertGreater(fail_lower, F(99, 100))
        self.assertTrue(np.all(~outside | fail))
        OBS['no_automatic_retry'] = dict(N=2, time='1/8',
            graph_H_zero=True, full_scalar_removed='-1', outside_lower=str(out_lower),
            fail_lower=str(fail_lower),
            conditional_outside_given_failure_strict_lower='1/100',
            same_graph_state_not_reset=True, general_retry_impossibility_claim=False)

    def test_06_uniform_budget_and_internal_resources(self):
        cmax, f, v, lip = 3, 24, 7, 48
        r = F(v*v*lip, 4)+F(v**4, 12)
        c = F(1, 6)
        t, delay = F(1, 131072), F(1, 1048576)
        self.assertLessEqual(v*t, 1)
        self.assertLessEqual(t, c/(2*r))
        s0 = c*c*t**6/4
        error_amplitude = 2*(r/c+f)*t
        eps = error_amplitude**2
        delayed = (error_amplitude+f*delay)**2
        relative_error = F(1, 1000)
        delta = s0*relative_error
        actual = (delayed+relative_error)/(1-relative_error)
        self.assertLess(actual, F(1, 100))
        self.assertEqual(s0, F(1, 730166745731460135262101046296576))
        # For finite-model source certificate, a less severe but still explicit budget.
        finite_s0 = F(1, 10**13)
        finite_delta = finite_s0/1000
        self.assertEqual(finite_delta, F(1, 10**16))
        finite_delay = F(1, 4096)
        finite_delayed = (F(1, 100)+f*finite_delay)**2
        finite_actual = (finite_delayed+relative_error)/(1-relative_error)
        self.assertLess(finite_actual, F(1, 100))
        n = 6
        data_bits = (n+1).bit_length()  # ceil(log2(N+2))
        id_bits = (n-1).bit_length()
        reader_bits = (n-1).bit_length()  # N-1 returns plus one failure
        listed = 2*n*data_bits+n*id_bits+n+reader_bits+n*(n-1)//2
        self.assertEqual(listed, 78)
        OBS['uniform_and_resource_budget'] = dict(C=cmax, f=f, speed_norm=v, lipschitz=lip,
            R=str(r), c=str(c), time=str(t), success_lower=str(s0),
            joint_error_upper=str(eps), delay=str(delay), delayed_joint_error_upper=str(delayed),
            total_complete_output_half_trace_budget=str(delta),
            actual_success_lower=str(s0-delta), actual_joint_error_upper=str(actual),
            actual_joint_error_strict_upper='1/100',
            finite_six_model_success_lower=str(finite_s0),
            finite_six_model_complete_error_budget=str(finite_delta),
            finite_six_model_delay=str(finite_delay),
            finite_six_model_actual_joint_error_upper=str(finite_actual),
            local_packet_dimension=8, packet_qubits=18, old_data_memory_qubits=18,
            trusted_identity_qubits=18, root_role_qubits=6, root_reader_qubits=3,
            relation_qubits=15, listed_qubits=listed,
            preparation_clock_phase_reference_and_coupling_hardware_not_generated=True,
            no_high_success_retry_or_experimental_sampling_claim=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Checks))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    dependency = HERE/'current_neighbor_detection_results.json'
    return dict(round=501, scientific_baseline_round=500,
        reused_frozen_rounds=[429, 440, 459, 472, 500], tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__,
        dependency_results_sha256={dependency.name: hashlib.sha256(dependency.read_bytes()).hexdigest()},
        scope=dict(local_root_readout_only=True, single_unknown_graph_and_passive_reference=True,
            root_degree_may_be_unknown=True, random_label_joint_error_bound=True,
            per_label_posterior_bound=False, new_multilevel_data_and_local_identity_interaction=True,
            original_qubit_H_implemented=False, autonomous_identity_source=False,
            automatic_failure_reset=False, high_success_retry_proved=False,
            physical_positions_or_dimension_generated=False,
            full_GR_goal_completed=False, phase_closure_triggered=False), observations=OBS.copy())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
