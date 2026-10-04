"""Round 465: actual data reception under autonomous tree rewiring.

Scientific baseline 463, independent of round 464. The new propagation
contract allows unknown source--graph--reference correlations but prepares
all other data qubits in zero. No common fixed chain, threshold graph,
measurement of graph distance, or replacement by classical graph dynamics.
"""
import argparse
from fractions import Fraction as Q
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

import branching_tree_distance_audit as old

TARGET = Path(__file__).with_name('relational_propagation_bound_audit_results.json')
OBS = {}


def short(x):
    return float(f'{float(x):.12g}')


def unitary(h, t):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*t*values))@vectors.conj().T


def trace_distance(a, b):
    return float(np.abs(np.linalg.eigvalsh(a-b)).sum()/2)


def distances(tree, n):
    return np.array([[len(old.path_between(tree, n, a, b))-1
                      for b in range(n)] for a in range(n)])


def sector(j=1, kappa=1):
    trees, _, full, flip, _ = old.six_vertex_sector()
    n, ng = 6, len(trees)
    hopping = np.zeros((n*ng, n*ng))
    diagonal = np.zeros_like(hopping)
    for g, tree in enumerate(trees):
        adj = old.adjacency(tree, n)
        for v in range(n):
            diagonal[v*ng+g, v*ng+g] = -len(adj[v])
            for w in adj[v]:
                hopping[v*ng+g, w*ng+g] = 1
    graph = np.kron(np.eye(n), flip)
    h = kappa*graph+j*(diagonal+hopping)
    h_no_hop = kappa*graph+j*diagonal
    # Vacuum block shares the graph evolution, with the same constant removed.
    hv = np.zeros(((n+1)*ng, (n+1)*ng))
    hv[:ng, :ng] = kappa*flip
    hv[ng:, ng:] = h
    d = np.array([len(old.path_between(tree, n, v, 5))
                  for v in range(n) for tree in trees], dtype=float)
    return trees, full, flip, h, h_no_hop, hopping, hv, d


def anti_weighted(h, d):
    hw = h*d[None, :]/d[:, None]
    return (hw-hw.conj().T)/(2j)


def receiver_output(vector, reference_dim):
    # vector shape is (vacuum + six excitation positions, graph, reference).
    arr = vector.reshape(7, 6, reference_dim)
    # Environment distinguishes vacuum from each nonreceiver excitation.
    q = np.zeros((2, 6, 6, reference_dim), complex)
    q[0, 0] = arr[0]
    q[1, 0] = arr[6]
    for v in range(5):
        q[0, v+1] = arr[v+1]
    mat = q.transpose(0, 3, 1, 2).reshape(2*reference_dim, -1)
    return mat@mat.conj().T


def source_embedding(source_graph_reference, far):
    # Unknown qubit at label 0; all other data fixed zero.
    nr = source_graph_reference.shape[-1]
    out = np.zeros((7, 6, nr), complex)
    out[0, far] = source_graph_reference[0]
    out[1, far] = source_graph_reference[1]
    return out.reshape(42, nr)


class Audit(unittest.TestCase):
    def close(self, a, b, tol=2e-11):
        self.assertLess(np.linalg.norm(a-b), tol)

    def test_01_arbitrary_vertex_distance_flip_count(self):
        total_trees, total_pairs, total_transitions = 0, 0, 0
        max_weighted = Q(0)
        for n in range(2, 7):
            for word in itertools.product(range(n), repeat=n-2):
                tree = old.prufer_tree(n, word)
                adj = old.adjacency(tree, n)
                c = max(map(len, adj))
                initial = distances(tree, n)
                changed = np.zeros((n, n), dtype=int)
                row = [[Q(0) for _ in range(n)] for _ in range(n)]
                flips = old.tree_flips(tree, n)
                for target, multiplicity in flips.items():
                    self.assertEqual(multiplicity, 1)
                    after = distances(target, n)
                    self.assertLessEqual(int(np.abs(after-initial).max()), 1)
                    changed += after != initial
                    removed = tree-target
                    for a, b in itertools.product(range(n), repeat=2):
                        if after[a, b] == initial[a, b]:
                            continue
                        path = old.path_between(tree, n, a, b)
                        path_edges = {old.edge(u, v) for u, v in zip(path[:-1], path[1:])}
                        self.assertTrue(removed & path_edges)
                        x, y = int(initial[a, b]+1), int(after[a, b]+1)
                        row[a][b] += abs(Q(x, y)-Q(y, x))/2
                    total_transitions += 1
                self.assertTrue(np.all(changed <= 2*initial*(c-1)**2))
                local_max = max(x for r in row for x in r)
                self.assertLessEqual(local_max, 2*(c-1)**2)
                max_weighted = max(max_weighted, local_max)
                total_trees += 1
                total_pairs += n*n
        OBS['all_vertex_flip_count'] = dict(
            vertices_range=[2, 6], labeled_trees=total_trees,
            ordered_vertex_pairs=total_pairs, directed_flips=total_transitions,
            no_leaf_condition=True, change_per_flip_at_most_one=True,
            changed_distance_requires_removed_path_edge=True,
            flip_count_bound='2*d_G(v,b)*(C-1)^2',
            graph_weighted_row_bound='2*(C-1)^2',
            maximum_observed_weighted_row=str(max_weighted),
            arbitrary_size_claim_from_analytic_proof=True)

    def test_02_exact_sector_and_weighted_generators(self):
        trees, full, flip, h, h0, hop, hv, d = sector()
        ng = 6
        ids = list(range(ng))+[(1 << (5-v))*ng+g for v in range(6) for g in range(ng)]
        embed = np.eye(384)[:, ids]
        # Old full H removed only mu*5; remove the additional J*5 scalar.
        self.close((full-5*np.eye(384))@embed, embed@hv)
        self.close(h-h0, hop)
        constants = []
        for j, kappa in ((1, 1), (Q(1, 2), -1), (-2, Q(1, 3))):
            _, _, _, hh, hh0, hopping, _, dd = sector(float(j), float(kappa))
            bg = 2*abs(kappa)*4
            bound = bg+Q(9, 4)*abs(j)
            actual = np.linalg.norm(anti_weighted(hh, dd), 2)
            self.assertLessEqual(actual, float(bound)+1e-12)
            self.assertLessEqual(np.linalg.norm(anti_weighted(hh0, dd), 2), float(bg)+1e-12)
            weighted_hop = float(j)*hopping*dd[None, :]/dd[:, None]
            self.assertLessEqual(np.linalg.norm(weighted_hop, 2), 6*abs(float(j))+1e-12)
            constants.append(dict(J=str(j), kappa=str(kappa), B=str(bound),
                                  antihermitian_norm=short(actual)))
        OBS['exact_sector_and_generators'] = dict(
            full_dimension=384, vacuum_plus_excitation_dimension=42,
            full_intertwining_error=short(np.linalg.norm((full-5*np.eye(384))@embed-embed@hv)),
            data_positions_can_be_internal_vertices=True,
            D='d_G(v,b)+1', generator_bound='B=2*|kappa|*(C-1)^2+3*C*|J|/4',
            weighted_hopping_norm_bound='2*C*|J|', cases=constants)

    def test_03_actual_reception_operator_bound(self):
        trees, _, _, h, h0, hop, hv, d = sector()
        t, c, bound = .02, 3, 41/4
        u = unitary(h, t)
        far = np.array([g for g, tree in enumerate(trees)
                        if len(old.path_between(tree, 6, 0, 5))-1 >= 3])
        received = u[30:36, far]
        eta = math.exp(bound*t)*min(1, 2*c*t)/4
        maximum = np.linalg.norm(received, 2)**2
        self.assertLessEqual(maximum, eta**2)
        # Whole graph input, with its initial inverse distance retained.
        k = u[30:36, :6]
        alpha = math.exp(bound*t)*min(1, 2*c*t)
        positive_difference = alpha**2*np.diag(1/d[:6]**2)-k.conj().T@k
        self.assertGreater(np.linalg.eigvalsh(positive_difference).min(), -1e-12)
        _, _, _, zero_h, _, _, _, _ = sector(0, 1)
        self.close(unitary(zero_h, t)[30:36, :6], np.zeros((6, 6)))
        OBS['actual_reception'] = dict(t=str(Q(1, 50)), C=3, J=1, kappa=1,
            initial_distance_minimum=3, unknown_far_graph_dimension=len(far),
            largest_reception_probability=short(maximum), probability_bound=short(eta**2),
            arbitrary_graph_input_weighted_inequality_min_eigenvalue=short(np.linalg.eigvalsh(positive_difference).min()),
            bound_formula='p_b <= exp(2*B*abs(t))*min(1,2*C*abs(J*t))^2*Tr[rho_G/(d_G(a,b)+1)^2]',
            J_zero_exact_no_reception=True, graph_not_measured=True)

    def test_04_unknown_source_graph_reference_and_channels(self):
        trees, _, _, h, h0, hop, hv, d = sector()
        far = np.array([g for g, tree in enumerate(trees)
                        if len(old.path_between(tree, 6, 0, 5))-1 >= 3])
        rng = np.random.default_rng(465)
        nr = 3
        psi = rng.normal(size=(2, len(far), nr))+1j*rng.normal(size=(2, len(far), nr))
        psi /= np.linalg.norm(psi)
        rho_r = np.einsum('sgr,sgt->rt', psi, psi.conj())
        gram_s = np.einsum('sgr,tgr->st', psi, psi.conj())
        self.assertLess(np.trace(gram_s@gram_s).real, .99)
        t = .05
        u = unitary(hv, t)
        replacer = np.kron(np.diag([1, 0]), rho_r)
        effects = {
            'identity': [np.eye(2)],
            'phase': [np.diag([1, -1])],
            'amplitude_damping_half': [np.diag([1, 1/math.sqrt(2)]),
                                       np.array([[0, 1/math.sqrt(2)], [0, 0]])]}
        outputs, rows = {}, []
        for name, kraus in effects.items():
            output = np.zeros((2*nr, 2*nr), complex)
            for e in kraus:
                prepared = np.einsum('st,tgr->sgr', e, psi)
                evolved = u@source_embedding(prepared, far)
                output += receiver_output(evolved, nr)
            outputs[name] = output
            self.close(output.reshape(2, nr, 2, nr).trace(axis1=0, axis2=2), rho_r)
            p = np.trace(output[nr:, nr:]).real
            distance = trace_distance(output, replacer)
            self.assertLessEqual(distance, p+math.sqrt(p*(1-p))+1e-12)
            rows.append(dict(channel=name, probability=short(p),
                             receiver_reference_replacer_distance=short(distance),
                             actual_p_block_bound=short(p+math.sqrt(p*(1-p)))))
        eta = math.exp((41/4)*t)*min(1, 6*t)/4
        for name in ('phase', 'amplitude_damping_half'):
            value = trace_distance(outputs['identity'], outputs[name])
            self.assertLessEqual(value, min(1, 2*(eta+eta**2))+1e-12)
            rows.append(dict(comparison='identity versus '+name, trace_distance=short(value)))
        OBS['unknown_source_graph_reference'] = dict(
            initial_source_graph_reference_dimensions=[2, 4, nr],
            source_graph_reference_product_not_assumed=True,
            all_other_data_initially_zero=True,
            receiver_reads_its_own_fixed_subject_factor=True,
            no_graph_localization_channel_or_postselection=True,
            encoding_channels=['identity', 'phase', 'amplitude_damping_half'],
            channel_bound=short(min(1, 2*(eta+eta**2))), values=rows)

    def test_05_coherent_message_requires_square_root_scale(self):
        # Same model's two-vertex tree: F=0, no graph dynamics. This only
        # audits the message quantifier, not a new communication mechanism.
        swap = old.core.swap(2, 0, 1)
        t = .2
        u = unitary(swap, t)
        plus = np.array([1, 0, 1, 0], complex)/math.sqrt(2)
        minus = np.array([1, 0, -1, 0], complex)/math.sqrt(2)
        outputs = []
        for vec in (plus, minus):
            out = (u@vec).reshape(2, 2)
            outputs.append(out.T@out.conj())
        distance = trace_distance(*outputs)
        p_one = math.sin(t)**2
        self.assertAlmostEqual(distance, abs(math.sin(t)), places=13)
        self.assertGreater(distance, p_one)
        OBS['coherent_message_boundary'] = dict(
            model='two-subject tree with same data SWAP, F=0',
            t='1/5', excited_input_reception_probability=short(p_one),
            plus_minus_message_trace_distance=short(distance),
            exact_formula='D_plus_minus=|sin(J*t)|=sqrt(p_one)',
            reception_probability_not_general_trace_distance=True)

    def test_06_size_independent_rational_certificate(self):
        c, j, kappa, d0, t = 3, Q(1), Q(1), 100, Q(1, 100)
        bound = 2*kappa*(c-1)**2+Q(3, 4)*c*j
        eta_upper = min(Q(1), 2*c*j*t)/(d0+1)/(1-bound*t)
        self.assertEqual(bound, Q(41, 4))
        self.assertEqual(eta_upper, Q(24, 36259))
        probability = eta_upper**2
        two_encoding = 2*(eta_upper+probability)
        self.assertLess(probability, Q(1, 2000000))
        self.assertLess(two_encoding, Q(7, 5000))
        OBS['arbitrary_size_certificate'] = dict(
            C=c, J=str(j), kappa=str(kappa), minimum_initial_distance=d0,
            t=str(t), B=str(bound), amplitude_upper=str(eta_upper),
            probability_upper=str(probability), two_encoding_trace_distance_upper=str(two_encoding),
            probability_less_than='1/2000000', two_encoding_distance_less_than='7/5000',
            no_large_Hilbert_matrix_built=True,
            continuous_time_exact_zero_cone_claimed=False,
            logarithmic_distance_time_bound_not_a_linear_light_cone=True)


def run():
    OBS.clear()
    out = io.StringIO()
    result = unittest.TextTestRunner(stream=out).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(out.getvalue())
    return dict(round=465, baseline_round=463, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(same_autonomous_H_as_rounds_440_to_443=True,
            initial_bounded_degree_tree_family_is_input=True,
            all_non_source_data_zero_is_additional_contract=True,
            arbitrary_source_graph_reference_state_allowed=True,
            arbitrary_local_source_CPTP_encodings_allowed=True,
            arbitrary_unknown_background_data_bound_proved=False,
            predetermined_common_chain_used=False, threshold_used=False,
            graph_distance_measured_or_postselected=False,
            internal_preparation_or_readout_device_derived=False,
            spatial_dimension_derived=False, full_GR_goal_completed=False,
            independent_from_round_464=True, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert result == json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
