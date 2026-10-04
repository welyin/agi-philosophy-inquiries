"""Round 500: random current-neighbor output on one unknown dynamic graph.

The guarantee is a JOINT cq event, not a posterior bound for every label.
Exact identities and rational Taylor enclosures certify finite examples;
floating-point stopped-instrument checks are diagnostics only.
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
import dynamic_neighborhood_envelope as frozen499

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'current_neighbor_detection_results.json'
OBS = {}


def constants(c, j=1, kappa=1):
    assert c >= 1 and j != 0
    f = 2*c*(c-1)**2
    return f, abs(kappa)*f+abs(j)*c, 2*abs(kappa)*f+4*abs(j)*c


def six_model():
    trees, old_h, full, flip = frozen499.six_model()
    return trees, old_h-5*np.eye(36, dtype=np.int64), full, flip


def masks(trees, n, root):
    g = len(trees)
    p = np.repeat(np.arange(n) == root, g)
    near = np.array([v != root and tuple(sorted((root, v))) in tree
                     for v in range(n) for tree in trees], dtype=bool)
    bad = ~p & ~near
    return p, near, bad


def evolution(h, t):
    values, vectors = np.linalg.eigh(h.astype(float))
    return (vectors*np.exp(-1j*values*float(t)))@vectors.T.conj()


def exact_taylor(h, denominator=64, degree=12):
    """Integer real/imaginary numerators and a rational operator tail."""
    h = h.astype(object)
    order = len(h)
    scale = denominator**degree*math.factorial(degree)
    real = np.zeros((order, order), dtype=object)
    imag = real.copy()
    power = np.eye(order, dtype=object)
    for k in range(degree+1):
        coef = denominator**(degree-k)*math.factorial(degree)//math.factorial(k)
        part = real if k % 2 == 0 else imag
        part += ((-1)**((k+1)//2))*coef*power
        if k < degree:
            power = power@h
    norm_bound = max(max(sum(abs(int(x)) for x in row) for row in h),
                     max(sum(abs(int(x)) for x in row) for row in h.T))
    x = F(norm_bound, denominator)
    assert x < degree+2
    tail = x**(degree+1)/(math.factorial(degree+1)*(1-x/(degree+2)))
    return real, imag, scale, tail, norm_bound


def exact_effect_bounds(real, imag, scale, source, effect, tail):
    r = real[np.ix_(np.flatnonzero(effect), source)]
    z = imag[np.ix_(np.flatnonzero(effect), source)]
    ar, ai = r.T@r+z.T@z, r.T@z-z.T@r
    lows, highs = [], []
    for k in range(len(source)):
        radius = sum(abs(int(ar[k, l]))+abs(int(ai[k, l]))
                     for l in range(len(source)) if l != k)
        lows.append(F(int(ar[k, k])-radius, scale**2))
        highs.append(F(int(ar[k, k])+radius, scale**2))
    error = 2*tail+tail*tail
    return min(lows)-error, max(highs)+error


class Checks(unittest.TestCase):
    def test_01_original_restriction_and_blocks(self):
        trees, h, full, flip = six_model()
        self.assertEqual(full.shape, (384, 384))
        self.assertTrue(np.array_equal(h, h.T))
        rows = []
        for root in range(6):
            p, near, bad = masks(trees, 6, root)
            q = ~p
            t = h[np.ix_(q, p)]
            degrees = np.array([sum(root in e for e in tree) for tree in trees])
            self.assertTrue(np.array_equal(t.T@t, np.diag(degrees)))
            self.assertFalse(np.any(h[np.ix_(bad, p)]))
            self.assertTrue(np.all(p.astype(int)+near+bad == 1))
            degree_charge = np.tile(degrees, 6)
            self.assertFalse(np.any(h*(degree_charge[:, None]-degree_charge[None, :])))
            rows.append(dict(root=root, degree_values=sorted(set(map(int, degrees)))))
        OBS['exact_original_blocks'] = dict(full_dimension=384,
            one_excitation_dimension=36, common_scalar_removed=5,
            root_checks=rows, T_star_T_equals_root_degree=True, B_h_P_zero=True)

    def test_02_NNI_root_activity_and_degrees(self):
        graph_count = transition_count = 0
        maximum_touch = 0
        for n in range(3, 7):
            for word in itertools.product(range(n), repeat=n-2):
                tree = trees_model.prufer_tree(n, word)
                degrees = [sum(v in e for e in tree) for v in range(n)]
                c = max(degrees)
                f = constants(c)[0]
                touches = [0]*n
                graph_count += 1
                for other, weight in trees_model.tree_flips(tree, n).items():
                    self.assertEqual(degrees, [sum(v in e for e in other) for v in range(n)])
                    vertices = set().union(*map(set, tree.symmetric_difference(other)))
                    self.assertEqual(len(vertices), 4)
                    for v in vertices:
                        touches[v] += weight
                    transition_count += 1
                self.assertLessEqual(max(touches), f)
                maximum_touch = max(maximum_touch, max(touches))
        OBS['finite_NNI_crosscheck'] = dict(labelled_trees=graph_count,
            directed_transitions=transition_count, max_root_activity=maximum_touch,
            all_vertex_degrees_preserved=True, general_bound_proved_in_note=True)

    def test_03_rational_all_input_effect_certificate(self):
        trees, h, _, _ = six_model()
        real, imag, scale, tail, norm_bound = exact_taylor(h)
        t = F(1, 64)
        _, b, k = constants(3)
        self.assertLessEqual(t, F(1, k))
        rows = []
        for root in range(6):
            p, _, bad = masks(trees, 6, root)
            source = np.flatnonzero(p)
            lower, _ = exact_effect_bounds(real, imag, scale, source, ~p, tail)
            _, upper_bad = exact_effect_bounds(real, imag, scale, source, bad, tail)
            degree = sum(root in e for e in trees[0])
            self.assertGreater(lower, degree*t*t/4)
            self.assertLess(upper_bad/lower, b*b*t*t)
            rows.append(dict(root=root, click_lower=str(lower),
                bad_upper=str(upper_bad), joint_error_upper=str(upper_bad/lower),
                diagnostic_click_lower=float(lower),
                diagnostic_joint_error_upper=float(upper_bad/lower)))
        OBS['rational_full_input_certificate'] = dict(time=str(t), degree=12,
            operator_norm_bound=norm_bound, unitary_tail=str(tail),
            arithmetic='Python unbounded integers and Fraction; Gershgorin effects',
            root_certificates=rows, arbitrary_graph_and_passive_reference=True)

    def test_04_label_scope_and_unknown_degree(self):
        # All 16 labelled trees on four vertices include several root degrees.
        n = 4
        trees = sorted({trees_model.prufer_tree(n, w)
                        for w in itertools.product(range(n), repeat=n-2)}, key=lambda g: sorted(g))
        g = len(trees)
        flip = np.zeros((g, g), dtype=np.int64)
        for a, tree in enumerate(trees):
            for other, weight in trees_model.tree_flips(tree, n).items():
                flip[trees.index(other), a] = weight
        h = np.kron(np.eye(n, dtype=np.int64), flip)
        for a, tree in enumerate(trees):
            adj = np.zeros((n, n), dtype=np.int64)
            for u, v in tree:
                adj[u, v] = adj[v, u] = 1
            ix = np.arange(n)*g+a
            h[np.ix_(ix, ix)] += adj-np.diag(adj.sum(axis=0))
        p, _, bad = masks(trees, n, 0)
        degree = np.array([sum(0 in e for e in tree) for tree in trees])
        charge = np.tile(degree, n)
        self.assertFalse(np.any(h*(charge[:, None]-charge[None, :])))
        up = evolution(h, F(1, 64))[:, p]
        click = up[~p].conj().T@up[~p]
        wrong = up[bad].conj().T@up[bad]
        p0 = F(1, 4*64**2)
        eps = F(constants(3)[1]**2, 64**2)
        self.assertGreater(np.linalg.eigvalsh(click-float(p0)*np.diag(degree))[0], 0)
        self.assertGreater(np.linalg.eigvalsh(float(eps)*click-wrong)[0], 0)
        # z^3-3z+2=(z-1)^2(z+2); unit |z|, 0<t<2pi => nonzero.
        polynomial = np.polynomial.polynomial.polymul([1, -2, 1], [2, 1])
        self.assertEqual(polynomial.tolist(), [2, -3, 0, 1])
        path_h = np.array([[-1, 1, 0], [1, -2, 1], [0, 1, -1]])
        time = 1/64
        amp = (np.exp(1j*time)-1)**2*(np.exp(1j*time)+2)/6
        self.assertLess(abs(evolution(path_h, time)[2, 0]-amp), 1e-14)
        self.assertGreater(abs(amp)**2, 0)
        OBS['scope_and_degree'] = dict(root_degree_values=sorted(set(map(int, degree))),
            cross_degree_coherence_allowed=True, exact_degree_charge_commutes=True,
            four_vertex_effect_checks_are_diagnostics=True,
            path_counterexample_time='1/64', path_nonedge_click_probability=float(abs(amp)**2),
            fixed_label_2_error_exact=1, polynomial_identity_exact=True,
            per_label_posterior_bound_claimed=False)

    def test_05_same_graph_stopped_instrument_and_delay(self):
        trees, h, _, _ = six_model()
        root, m, time, delay = 0, 7, F(1, 64), F(1, 256)
        u, ud = evolution(h, time), evolution(h, delay)
        blocks = [u[v*6:(v+1)*6, root*6:(root+1)*6] for v in range(6)]
        psi = np.array([[1, 1j], [2j, -1], [1+1j, 2], [-1j, 1], [2, -2j], [1, 3j]], complex)
        psi /= np.linalg.norm(psi)
        initial_reference = psi.conj().T@psi
        ref = np.zeros((2, 2), complex)
        complete = np.zeros((6, 6), complex)
        prefix = np.eye(6, dtype=complex)
        success = wrong = delayed_wrong = 0.0
        for attempt in range(1, m+1):
            for label in range(1, 6):
                kraus = blocks[label]@prefix
                complete += kraus.conj().T@kraus
                branch = kraus@psi
                probability = float(np.linalg.norm(branch)**2)
                missing = np.array([tuple(sorted((root, label))) not in tree for tree in trees])
                success += probability
                wrong += float(np.linalg.norm(branch[missing])**2)
                ref += branch.conj().T@branch
                post = np.zeros((36, 2), complex)
                post[label*6:(label+1)*6] = branch
                delivered = ud@post
                delayed_wrong += float(np.linalg.norm(delivered[np.tile(missing, 6)])**2)
            prefix = blocks[root]@prefix
        abort_state = prefix@psi
        abort = float(np.linalg.norm(abort_state)**2)
        complete += prefix.conj().T@prefix
        ref += abort_state.conj().T@abort_state
        _, b, _ = constants(3)
        p0, eps = float(time*time/4), float(b*b*time*time)
        self.assertLess(np.linalg.norm(complete-np.eye(6)), 5e-13)
        self.assertAlmostEqual(success+abort, 1, places=12)
        self.assertLess(np.linalg.norm(ref-initial_reference), 5e-13)
        self.assertLessEqual(abort, (1-p0)**m)
        self.assertLessEqual(wrong/success, eps)
        self.assertLessEqual(delayed_wrong/success, float((b*time+24*delay)**2))
        OBS['stopped_instrument_diagnostics'] = dict(attempts=m, graph_copies=1,
            time_per_attempt=str(time), delay_after_each_actual_stop=str(delay),
            success=success, abort=abort, joint_wrong_given_success=wrong/success,
            delayed_fixed_record_wrong_given_success=delayed_wrong/success,
            completeness_residual=float(np.linalg.norm(complete-np.eye(6))),
            passive_reference_nonselective_residual=float(np.linalg.norm(ref-initial_reference)),
            graph_reset=False, common_terminal_time_graph_freeze=False,
            retained_first_click_and_abort_records=True)

    def test_06_finite_resource_and_implementation_certificate(self):
        f, b, k = constants(3)
        time, delay = F(1, 256), F(1, 4096)
        p0 = time*time/4
        m = int(3/p0)
        eps = b*b*time*time
        delayed = (b*time+f*delay)**2
        delta = F(1, 10000)
        self.assertEqual((f, b, k), (24, 27, 60))
        self.assertLessEqual(time, F(1, k))
        self.assertGreater(sum(F(3**v, math.factorial(v)) for v in range(9)), 20)
        self.assertEqual(m, 786432)
        self.assertLess(delayed, F(1, 80))
        actual_wrong = (F(1, 80)+delta)/(F(19, 20)-delta)
        self.assertLess(actual_wrong, F(1, 75))
        # A complete-instrument budget, never just a selected successful branch.
        source_error = delta/2
        per_attempt_error = delta/(2*m)
        self.assertEqual(source_error+m*per_attempt_error, delta)
        OBS['finite_resource_certificate'] = dict(C=3, J=1, kappa=1, f=f, b=b, K=k,
            time=str(time), p0=str(p0), attempts=m, ideal_abort_strict_upper='1/20',
            ideal_joint_error_upper=str(eps), delay=str(delay),
            delayed_joint_error_upper=str(delayed), total_complete_process_error=str(delta),
            source_half_trace_budget=str(source_error),
            each_complete_attempt_half_diamond_budget=str(per_attempt_error),
            actual_success_lower=str(F(19, 20)-delta), actual_joint_error_upper=str(actual_wrong),
            actual_joint_error_strict_upper='1/75',
            natural_waiting_sum_upper=str(m*time),
            raw_reader_bits_for_N_6=m*6, internal_data_memory_qubits_for_N_6=6,
            original_relation_qubits_for_N_6=15,
            controller_clock_preparation_isolation_not_generated=True,
            massive_protocol_actually_run=False)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Checks))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    deps = ['dynamic_neighborhood_envelope_results.json',
            'quantum_neighborhood_signal_audit_results.json']
    return dict(round=500, scientific_baseline_round=499,
        reused_frozen_rounds=[440, 443, 470, 472, 493, 499],
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__,
        dependency_results_sha256={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        scope=dict(random_reported_label_joint_current_relation_certificate=True,
            individual_label_posterior_guarantee=False, arbitrary_unknown_graph_and_passive_reference=True,
            root_degree_need_not_be_known=True, same_graph_finite_retries=True,
            leaf_only_radius_one_envelope=True, initial_probe_and_complete_readout_are_inputs=True,
            actual_spatial_positions_or_dimension_generated=False,
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
