"""Round 268: paired edge refinement amplifies a cut at bounded degree.

Fresh channels and a marked matching are inputs. This is an atomic graph
rewrite, not a derivation of a physical channel or a distributed protocol.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np
from ancestral_quotient_audit import copy_state, fingerprint, graph_data
from vertex_split_cycle_audit import subcubic_graphs, edges, statistics
from fisher_hanoi_scaling_audit import bfs


def port_state(graph):
    return {(v,): {p: ((w,), graph[w].index(v)) for p, w in enumerate(row)}
            for v, row in enumerate(graph)}


def check_ports(state):
    for u, row in state.items():
        assert len(row) <= 3
        assert len({v for v, _ in row.values()}) == len(row)
        for p, (v, q) in row.items():
            assert u != v and state[v][q] == (u, p)


def refine(state, u, v, flips=(0, 0)):
    """Replace an existing marked pair; residual old ports use distinct children."""
    if u == v or u not in state or v not in state:
        raise ValueError('Distinct live endpoints are required.')
    if v not in {w for w, _ in state[u].values()}:
        raise ValueError('Endpoints must already share a primitive edge.')
    old = {u: dict(state[u]), v: dict(state[v])}
    residual = {a: [(p, w, q) for p, (w, q) in sorted(old[a].items()) if w not in (u, v)]
                for a in (u, v)}
    if any(len(row) > 2 for row in residual.values()):
        raise ValueError('Degree cap exceeded.')
    children = {a: (a+(0,), a+(1,)) for a in (u, v)}
    assert not any(w in state for row in children.values() for w in row)
    del state[u], state[v]
    for row in children.values():
        a, b = row
        state[a] = {0: (b, 0)}
        state[b] = {0: (a, 0)}
    for i in (0, 1):
        a, b = children[u][i], children[v][i]
        state[a][1] = (b, 1)
        state[b][1] = (a, 1)
    for a, flip in zip((u, v), flips):
        for i, (_, w, q) in enumerate(residual[a]):
            child = children[a][i ^ flip]
            state[child][2] = (w, q)
            state[w][q] = (child, 2)
    return state


def cut_size(state, roots):
    return sum((u[0] in roots) != (v[0] in roots)
               for u, row in state.items() for v, _ in row.values() if u < v)


@lru_cache(maxsize=None)
def exhaustive():
    graph_count = rewrite_count = cut_checks = 0
    for graph in subcubic_graphs():
        graph_count += 1
        before = statistics(graph)
        initial = port_state(graph)
        for u, v in edges(graph):
            for flips in itertools.product((0, 1), repeat=2):
                state = refine(copy_state(initial), (u,), (v,), flips)
                check_ports(state)
                after = statistics(graph_data(state)[0])
                assert after['vertices'] == before['vertices']+2
                assert after['edges'] == before['edges']+3
                assert after['cycle_rank'] == before['cycle_rank']+1
                assert after['components'] == before['components']
                for mask in range(1, 16):
                    roots = {i for i in range(5) if mask >> i & 1}
                    assert cut_size(state, roots) == cut_size(initial, roots)+int((u in roots) != (v in roots))
                    cut_checks += 1
                rewrite_count += 1
    return {'five_vertex_graphs': graph_count, 'paired_rewrites': rewrite_count,
            'ancestral_cut_increment_checks': cut_checks}


@lru_cache(maxsize=None)
def ladder_state(depth):
    state = {(0,): {1: ((1,), 1)}, (1,): {1: ((0,), 1)}}
    for _ in range(depth):
        pairs = [(u, row[1][0]) for u, row in sorted(state.items()) if u < row[1][0]]
        for u, v in pairs:
            refine(state, u, v)
    return state


def rail_order(state):
    """Recover a path order from adjacency, without input spatial coordinates."""
    side = {w for w in state if w[0] == 0}
    neighbors = {w: sorted(v for v, _ in state[w].values() if v in side) for w in side}
    if len(side) == 1:
        return sorted(side)
    start = min(w for w, row in neighbors.items() if len(row) == 1)
    order, previous = [start], None
    while True:
        following = [w for w in neighbors[order[-1]] if w != previous]
        if not following:
            return order
        previous = order[-1]
        order.append(following[0])
        assert len(order) <= len(side)


def ladder_certificate(depth):
    state = ladder_state(depth)
    first = rail_order(state)
    second = [state[w][1][0] for w in first]
    assert len(set(first+second)) == len(state)
    coordinate = {w: (side, i) for side, row in enumerate((first, second)) for i, w in enumerate(row)}
    actual = {tuple(sorted((coordinate[u], coordinate[v]))) for u, row in state.items()
              for v, _ in row.values() if u < v}
    length = len(first)
    expected = {((0, i), (1, i)) for i in range(length)}
    expected |= {((side, i), (side, i+1)) for side in (0, 1) for i in range(length-1)}
    assert actual == expected
    prefix = set(first[:length//2]+second[:length//2])
    prefix_cut = sum((u in prefix) != (v in prefix)
                     for u, row in state.items() for v, _ in row.values() if u < v)
    return {'depth': depth, **statistics(graph_data(state)[0]),
            'vertices_per_root': length, 'root_interface_edges': cut_size(state, {0}),
            'longitudinal_half_vertices': len(prefix), 'longitudinal_half_cut': prefix_cut,
            'diameter': length}


def data_embedding():
    # Old endpoint data remain at child zero. Fresh child data are initialized.
    j = np.zeros((16, 4), complex)
    for a, b in itertools.product((0, 1), repeat=2):
        j[8*a+2*b, 2*a+b] = 1
    return j


def report():
    return {'round': 268, 'exhaustive': exhaustive(),
            'paired_growth_ladder': [ladder_certificate(n) for n in range(1, 10)],
            'new_primitive_channels_per_rewrite_net': 3,
            'new_cross_ancestry_channels_per_selected_cross_edge': 1,
            'scope': 'A supplied two-endpoint channel-creation primitive increases selected ancestral cuts while preserving degree at most three. Persistent ports make disjoint marked pair rewrites commute abstractly. Repeated aligned refinement is exactly a ladder, hence does not select three-dimensional geometry; data isometry does not certify new channel capacity.'}


class Audit(unittest.TestCase):
    def test_01_graph_counts_and_degree_capacity(self):
        self.assertEqual(exhaustive()['five_vertex_graphs'], 768)
        self.assertGreater(exhaustive()['paired_rewrites'], 1000)

    def test_02_crossing_cut_amplification(self):
        self.assertEqual(exhaustive()['ancestral_cut_increment_checks'],
                         15*exhaustive()['paired_rewrites'])

    def test_03_disjoint_pairs_commute_even_with_external_adjacency(self):
        checked = 0
        for graph in subcubic_graphs():
            initial = port_state(graph)
            for first, second in itertools.combinations(edges(graph), 2):
                if set(first) & set(second):
                    continue
                a, b = copy_state(initial), copy_state(initial)
                for u, v in (first, second):
                    refine(a, (u,), (v,))
                for u, v in (second, first):
                    refine(b, (u,), (v,))
                self.assertEqual(fingerprint(a), fingerprint(b))
                checked += 1
        self.assertGreater(checked, 1000)

    def test_04_complete_marked_cut_doubles_interface(self):
        for depth in range(1, 10):
            row = ladder_certificate(depth)
            self.assertEqual(row['root_interface_edges'], 2**depth)
            self.assertEqual(row['vertices'], 2**(depth+1))
            self.assertLessEqual(row['max_degree'], 3)
            self.assertEqual(row['edges'], 3*2**depth-2)

    def test_05_independent_ladder_distances(self):
        for depth in range(1, 6):
            state = ladder_state(depth)
            first = rail_order(state)
            second = [state[w][1][0] for w in first]
            graph, words = graph_data(state)
            coord = {w: (a, i) for a, row in enumerate((first, second)) for i, w in enumerate(row)}
            for u, word in enumerate(words):
                distance = bfs(graph, u)
                a, i = coord[word]
                self.assertEqual(distance, [abs(i-coord[w][1])+abs(a-coord[w][0]) for w in words])

    def test_06_ladder_root_balls_and_transverse_cut(self):
        state = ladder_state(9)
        graph, words = graph_data(state)
        root = rail_order(state)[0]
        distance = np.array(bfs(graph, words.index(root)))
        for radius in (1, 2, 4, 16, 64, 128, 511):
            self.assertEqual(int(np.sum(distance <= radius)), 2*radius+1)
        for depth in range(1, 10):
            row = ladder_certificate(depth)
            self.assertEqual(row['longitudinal_half_cut'], 2)
            self.assertEqual(row['longitudinal_half_vertices'], 2**depth)

    def test_07_arbitrary_reference_data_preserved(self):
        j = data_embedding()
        np.testing.assert_allclose(j.conj().T@j, np.eye(4))
        rng = np.random.default_rng(268)
        psi = rng.normal(size=20)+1j*rng.normal(size=20)
        psi /= np.linalg.norm(psi)
        encoded = np.kron(np.eye(5), j)@psi
        decoded = np.kron(np.eye(5), j.conj().T)@encoded
        np.testing.assert_allclose(decoded, psi, atol=1e-14)

    def test_08_old_root_record_is_copied_only_as_a_basis_record(self):
        record = np.zeros((4, 2))
        record[0, 0] = record[3, 1] = 1
        np.testing.assert_allclose(record.T@record, np.eye(2))
        plus = np.ones(2)/np.sqrt(2)
        encoded = (record@plus).reshape(2, 2)
        np.testing.assert_allclose(encoded@encoded.T, np.eye(2)/2)
        self.assertFalse(np.allclose(encoded@encoded.T, np.outer(plus, plus)))

    def test_09_rewrite_requires_existing_pair_and_live_endpoints(self):
        state = port_state(((1,), (0,), ()))
        with self.assertRaises(ValueError):
            refine(copy_state(state), (0,), (2,))
        refine(state, (0,), (1,))
        with self.assertRaises(ValueError):
            refine(state, (0,), (1,))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result['runtime'] = {'python': platform.python_version(), 'numpy': np.__version__}
    result['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    payload = json.dumps(result, ensure_ascii=False, indent=2)+'\n'
    target = Path(__file__).with_name('paired_interface_growth_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
