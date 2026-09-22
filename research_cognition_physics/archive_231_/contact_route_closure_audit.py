"""Round 257: certified overlay contacts do not create primitive channels.

Introduction depth is a bookkeeping depth, not elapsed time. Shortest path
certificates retain their primitive hops; no new physical edge is inferred.
"""
import argparse
from collections import deque
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np


def adjacency(n, edges):
    result = [set() for _ in range(n)]
    for u, v in edges:
        result[u].add(v); result[v].add(u)
    return tuple(tuple(sorted(row)) for row in result)


def distances(graph):
    n = len(graph)
    result = np.full((n, n), np.inf)
    for start in range(n):
        result[start, start] = 0
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v in graph[u]:
                if np.isinf(result[start, v]):
                    result[start, v] = result[start, u]+1
                    queue.append(v)
    return result


def initial_routes(graph):
    return {(u, v): (u, v) for u, row in enumerate(graph) for v in row}


def introduce(routes):
    """Take one snapshot composition of certified routes, retaining shortest witnesses."""
    result = dict(routes)
    for (u, w), first in routes.items():
        for (w2, v), second in routes.items():
            if w2 != w or u == v:
                continue
            path = first+second[1:]
            # Loop erasure is sound: all remaining consecutive pairs were primitive edges.
            simple = []
            for node in path:
                if node in simple:
                    simple = simple[:simple.index(node)+1]
                else:
                    simple.append(node)
            path = tuple(simple)
            old = result.get((u, v))
            if old is None or (len(path), path) < (len(old), old):
                result[u, v] = path
    return result


def overlay(n, routes):
    return adjacency(n, routes)


def route_cost_distances(n, routes):
    d = np.full((n, n), np.inf)
    np.fill_diagonal(d, 0)
    for (u, v), path in routes.items():
        d[u, v] = len(path)-1
    for w in range(n):
        d = np.minimum(d, d[:, w, None]+d[None, w, :])
    return d


@lru_cache(maxsize=None)
def exhaustive():
    n = 5
    possible = tuple(itertools.combinations(range(n), 2))
    stages = 0
    for mask in range(1 << len(possible)):
        graph = adjacency(n, [e for k, e in enumerate(possible) if mask >> k & 1])
        base = distances(graph)
        routes = initial_routes(graph)
        for depth in range(4):
            observed = set(routes)
            expected = {(u, v) for u in range(n) for v in range(n)
                        if 0 < base[u, v] <= 2**depth}
            assert observed == expected, (mask, depth)
            for (u, v), path in routes.items():
                assert path[0] == u and path[-1] == v
                assert len(path)-1 == base[u, v]
                assert all(b in graph[a] for a, b in zip(path, path[1:]))
            np.testing.assert_array_equal(route_cost_distances(n, routes), base)
            np.testing.assert_array_equal(np.isfinite(distances(overlay(n, routes))), np.isfinite(base))
            stages += 1
            routes = introduce(routes)
    return {'graphs': 1 << len(possible), 'depths_per_graph': 4,
            'graph_depth_cases': stages, 'all_support_cost_and_component_checks_passed': True}


def chain_report():
    graph = adjacency(9, [(u, u+1) for u in range(8)])
    routes = initial_routes(graph)
    rows = []
    for depth in range(4):
        d = distances(overlay(9, routes))
        rows.append({'introduction_depth': depth, 'overlay_edges': len(routes)//2,
                     'unweighted_overlay_diameter': int(d.max()),
                     'primitive_cost_endpoint_distance': int(route_cost_distances(9, routes)[0, 8]),
                     'endpoint_contact_known': (0, 8) in routes,
                     'endpoint_route': list(routes.get((0, 8), ()))})
        routes = introduce(routes)
    return rows


def reduced_b(rho):
    return np.trace(rho.reshape(2, 2, 2, 2), axis1=0, axis2=2)


def report():
    return {'round': 257, 'exhaustive': exhaustive(), 'nine_vertex_path': chain_report(),
            'primitive_communication_edges_created': 0,
            'primitive_quantum_gates_created': 0,
            'scope': 'Classical introductions over supplied reliable primitive channels create certified overlay routes only. Connected components and primitive-hop distances are unchanged. Shared history or entanglement alone supplies neither signaling channels nor entangling gates.'}


class Audit(unittest.TestCase):
    def test_01_all_five_vertex_graphs(self):
        self.assertEqual(exhaustive()['graph_depth_cases'], 4096)

    def test_02_line_closed_forms(self):
        rows = chain_report()
        self.assertEqual([r['overlay_edges'] for r in rows], [8, 15, 26, 36])
        self.assertEqual([r['unweighted_overlay_diameter'] for r in rows], [8, 4, 2, 1])
        self.assertEqual([r['primitive_cost_endpoint_distance'] for r in rows], [8]*4)
        self.assertEqual(rows[-1]['endpoint_route'], list(range(9)))

    def test_03_introductions_do_not_join_disconnected_components(self):
        graph = adjacency(6, [(0, 1), (1, 2), (3, 4), (4, 5)])
        routes = initial_routes(graph)
        for _ in range(6):
            routes = introduce(routes)
        self.assertEqual(len(routes), 12)
        self.assertNotIn((0, 5), routes)

    def test_04_same_common_record_does_not_determine_channels(self):
        common_records = ('root_record',)*4
        first = adjacency(4, [(0, 1), (1, 2), (2, 3)])
        second = adjacency(4, [(0, 1), (2, 3)])
        self.assertEqual(len(set(common_records)), 1)
        self.assertTrue(np.isfinite(distances(first)[0, 3]))
        self.assertTrue(np.isinf(distances(second)[0, 3]))

    def test_05_overlay_support_is_relabeling_covariant(self):
        graph = adjacency(5, [(0, 1), (1, 2), (1, 3), (3, 4)])
        permutation = (3, 0, 4, 1, 2)
        renamed = adjacency(5, [(permutation[u], permutation[v])
                               for u, row in enumerate(graph) for v in row])
        a, b = initial_routes(graph), initial_routes(renamed)
        for _ in range(4):
            self.assertEqual({(permutation[u], permutation[v]) for u, v in a}, set(b))
            a, b = introduce(a), introduce(b)
        # Chosen shortest witness may depend on label tie-break; existence and costs do not.

    def test_06_entanglement_without_messages_cannot_signal(self):
        bell = np.array([1., 0., 0., 1.])/np.sqrt(2)
        rho = np.outer(bell, bell)
        identity = np.eye(2)
        for strength in (0., .2, .7, 1.):
            ks = (np.array([[1., 0.], [0., np.sqrt(1-strength)]]),
                  np.array([[0., np.sqrt(strength)], [0., 0.]]))
            output = sum((a := np.kron(k, identity))@rho@a.T for k in ks)
            np.testing.assert_allclose(reduced_b(output), identity/2, atol=2e-14)

    def test_07_conditioning_on_unreceived_result_is_not_signaling(self):
        bell = np.array([1., 0., 0., 1.])/np.sqrt(2)
        rho = np.outer(bell, bell)
        branches = []
        for value in (0, 1):
            p = np.diag([int(value == 0), int(value == 1)])
            k = np.kron(p, np.eye(2))
            branches.append(reduced_b(k@rho@k))
        np.testing.assert_allclose(sum(branches), np.eye(2)/2, atol=2e-14)
        self.assertGreater(float(np.max(np.abs(2*branches[0]-2*branches[1]))), .99)

    def test_08_overlay_does_not_grant_an_entangling_gate(self):
        # Target ZZ unitary maps a product input to an NPT state; a classical route
        # with separable ancillary resources allows only mixtures of product branches.
        psi = np.exp(-1j*np.pi/4*np.array([1, -1, -1, 1]))/2
        rho = np.outer(psi, psi.conj())
        partial_transpose = rho.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
        self.assertAlmostEqual(float(np.linalg.eigvalsh(partial_transpose).min()), -.5)


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
    target = Path(__file__).with_name('contact_route_closure_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
