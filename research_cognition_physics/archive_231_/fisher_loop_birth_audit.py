"""Round 260: Fisher-type loop birth with port and quantum resource accounting.

Triangle replacement and the ability to create its channels are additional
inputs. Quantum instruments select a supplied graph rewrite, not physical laws.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np
import vertex_split_cycle_audit as graphlib
from contact_route_closure_audit import adjacency, distances
from mutual_proposal_interaction_audit import I, X, Z, PLUS, tensor
from record_preserving_feedback_audit import weak


def embedding(gadget, flags=False):
    k = len(gadget)
    d = 4 if flags else 2
    result = np.zeros((d**k, d))
    for initial in range(d):
        r, logical = divmod(initial, 2) if flags else (0, initial)
        for rest in itertools.product((0, 1), repeat=k-1):
            bits = (logical,)+rest
            sign = (-1)**sum(bits[a]*bits[b] for a, b in graphlib.edges(gadget))
            index = 0
            for bit in bits:
                index = d*index+(2*r+bit if flags else bit)
            result[index, initial] = sign/np.sqrt(2**(k-1))
    return result


def instrument(sign, eta=.6, q=X, flags=False):
    gadget = graphlib.TRIANGLE if sign == 1 else graphlib.EDGE
    k = weak(q, sign, eta)
    if flags:
        k = np.kron(I, k)
    return embedding(gadget, flags)@k


def triangle_all(graph):
    # Neighbor order represents transported local port slots, not physical coordinates.
    slots = [{neighbor: i for i, neighbor in enumerate(row)} for row in graph]
    es = [(3*v+a, 3*v+b) for v in range(len(graph)) for a, b in graphlib.edges(graphlib.TRIANGLE)]
    es += [(3*u+slots[u][v], 3*v+slots[v][u]) for u, v in graphlib.edges(graph)]
    return adjacency(3*len(graph), es)


@lru_cache(maxsize=None)
def topology_audit():
    local_cases = graph_count = distance_pairs = 0
    for graph in graphlib.subcubic_graphs():
        graph_count += 1
        before = graphlib.statistics(graph)
        for v in range(len(graph)):
            for slots in itertools.permutations(range(3), len(graph[v])):
                assignment = dict(zip(graph[v], slots))
                output = graphlib.replace_vertex(graph, v, graphlib.TRIANGLE, assignment)
                after = graphlib.statistics(output)
                assert after['vertices'] == before['vertices']+2
                assert after['edges'] == before['edges']+3
                assert after['components'] == before['components']
                assert after['cycle_rank'] == before['cycle_rank']+1
                assert after['max_degree'] <= 3
                local_cases += 1
        decorated = triangle_all(graph)
        old_dist = distances(graph)
        new_dist = distances(decorated)
        collapse = np.arange(len(decorated))//3
        lifted = old_dist[np.ix_(collapse, collapse)]
        assert np.all(lifted <= new_dist)
        assert np.all(new_dist <= 2*lifted+1)
        distance_pairs += len(decorated)**2
    return {'five_vertex_graphs': graph_count, 'local_triangle_rewrites': local_cases,
            'single_decoration_distance_pairs': distance_pairs,
            'distance_bound': 'd_G(collapse(x),collapse(y)) <= d_F(x,y) <= 2*d_G(collapse(x),collapse(y))+1'}


def growth_summary():
    graph = ((),)
    loops = 0
    for step in range(24):
        v = step % len(graph)
        triangle = step % 3 != 0
        if triangle:
            assignment = dict(zip(graph[v], range(len(graph[v]))))
            graph = graphlib.replace_vertex(graph, v, graphlib.TRIANGLE, assignment)
            loops += 1
        else:
            graph = graphlib.replace_vertex(graph, v, graphlib.EDGE,
                                            next(graphlib.binary_assignments(graph, v)))
    return {'events': 24, 'triangle_events': loops, 'binary_events': 24-loops,
            'fresh_logical_qubits': 24+loops, 'new_internal_channels': 24+2*loops,
            'final': graphlib.statistics(graph)}


def report():
    eta = .6
    plus = instrument(1, eta)
    probabilities = {}
    for name, psi in (('X_plus', PLUS), ('X_minus', np.array([1., -1.])/np.sqrt(2)),
                      ('Z_plus', np.array([1., 0.]))):
        value = plus@psi
        probabilities[name] = float(np.vdot(value, value).real)
    maps = [instrument(sign, eta) for sign in (1, -1)]
    return {'round': 260, 'topology': topology_audit(), 'mixed_growth': growth_summary(),
            'eta': eta, 'one_event_cycle_birth_probabilities': probabilities,
            'instrument_completeness_max_error': float(np.max(np.abs(sum(b.T@b for b in maps)-I))),
            'expected_cycle_increment_X_plus': probabilities['X_plus'],
            'expected_new_logical_qubits_X_plus': 1+probabilities['X_plus'],
            'expected_new_internal_channels_X_plus': 1+2*probabilities['X_plus'],
            'scope': 'Known Fisher-type triangle gadgets supply local cycle creation in bounded-degree channel graphs. Supplied weak instruments control binary versus triangle birth, with explicit ancillas and channel-creation capability. Finite decoration obeys a distance bound; no dimension, continuum spacetime, unitarity on all graph sectors, or gravity is derived.'}


class Audit(unittest.TestCase):
    def test_01_exhaustive_triangle_capacity_and_cycles(self):
        self.assertGreater(topology_audit()['local_triangle_rewrites'], 15000)

    def test_02_single_decoration_metric_bound(self):
        self.assertEqual(topology_audit()['single_decoration_distance_pairs'], 172800)

    def test_03_mixed_growth_resource_equalities(self):
        result = growth_summary()
        self.assertEqual(result['triangle_events'], 16)
        self.assertEqual(result['final'], {'vertices': 41, 'edges': 56, 'components': 1,
                                          'cycle_rank': 16, 'max_degree': 3})

    def test_04_embeddings_preserve_entire_unknown_input(self):
        for gadget in (graphlib.EDGE, graphlib.TRIANGLE):
            for flags in (False, True):
                j = embedding(gadget, flags)
                np.testing.assert_allclose(j.T@j, np.eye(j.shape[1]), atol=3e-14)
                d = j.shape[1]
                maximally_entangled = np.eye(d).reshape(-1)/np.sqrt(d)
                encoded = np.kron(j, np.eye(d))@maximally_entangled
                decoded = np.kron(j.T, np.eye(d))@encoded
                np.testing.assert_allclose(decoded, maximally_entangled, atol=3e-14)

    def test_05_full_instrument_is_complete(self):
        for eta in (0., .6, 1.):
            for flags in (False, True):
                bs = [instrument(sign, eta, flags=flags) for sign in (1, -1)]
                np.testing.assert_allclose(sum(b.T@b for b in bs), np.eye(bs[0].shape[1]), atol=3e-14)

    def test_06_old_record_on_each_child(self):
        old_flag = np.kron(Z, I)
        for gadget in (graphlib.EDGE, graphlib.TRIANGLE):
            j = embedding(gadget, flags=True)
            for at in range(len(gadget)):
                new_flag = tensor([old_flag if pos == at else np.eye(4) for pos in range(len(gadget))])
                np.testing.assert_allclose(new_flag@j, j@old_flag, atol=3e-14)

    def test_07_state_dependent_loop_birth(self):
        actual = report()['one_event_cycle_birth_probabilities']
        for name, expected in [('X_plus', .8), ('X_minus', .2), ('Z_plus', .5)]:
            self.assertAlmostEqual(actual[name], expected)

    def test_08_disjoint_quantum_rewrites_commute(self):
        for first, second in itertools.product((1, -1), repeat=2):
            a, b = instrument(first), instrument(second)
            ab = np.kron(np.eye(a.shape[0]), b)@np.kron(a, I)
            ba = np.kron(a, np.eye(b.shape[0]))@np.kron(I, b)
            np.testing.assert_allclose(ab, ba, atol=3e-14)

    def test_09_readout_disturbance_remains_accounted(self):
        # Decode the known gadget after each branch, then discard its outcome.
        eta = .6
        decoded = [embedding(g).T@instrument(s, eta)
                   for s, g in ((1, graphlib.TRIANGLE), (-1, graphlib.EDGE))]
        output = sum(k@Z@k.T for k in decoded)
        np.testing.assert_allclose(output, .8*Z, atol=3e-14)
        self.assertAlmostEqual(float(sum(abs(np.trace(k))**2 for k in decoded)/4), .9)

    def test_10_no_triangle_with_fewer_than_three_simple_vertices(self):
        for n in (1, 2, 3):
            possible = tuple(itertools.combinations(range(n), 2))
            cyclic = []
            for mask in range(1 << len(possible)):
                graph = adjacency(n, [edge for k, edge in enumerate(possible) if mask >> k & 1])
                if graphlib.components(graph) == 1 and graphlib.statistics(graph)['cycle_rank'] > 0:
                    cyclic.append(graph)
            self.assertEqual(len(cyclic), int(n == 3))


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
    target = Path(__file__).with_name('fisher_loop_birth_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
