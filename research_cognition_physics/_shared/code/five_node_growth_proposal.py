"""Unnumbered review of the user's +5 node-count proposal; not a growth law.

Small counterexamples separate node counts from adjacency and recursive closure.
No physical dimension is inferred from these finite graphs.
"""
import argparse
from collections import Counter, deque
import hashlib
import json
from pathlib import Path
import platform
import numpy as np


def adjacency(n, edges):
    a = np.zeros((n, n), dtype=int)
    for u, v in edges:
        assert u != v
        a[u, v] = a[v, u] = 1
    return a


def seed():
    return adjacency(5, [(i, (i+1) % 5) for i in range(5)])


def expand_five(inner_step):
    # Both examples preserve exactly the same original five-cycle and add
    # five vertices, five spokes and five inner edges.
    edges = [(i, (i+1) % 5) for i in range(5)]
    edges += [(i, i+5) for i in range(5)]
    edges += [(i+5, (i+inner_step) % 5+5) for i in range(5)]
    return adjacency(10, edges)


def metrics(a):
    n = len(a)
    pair_distances = []
    for root in range(n):
        distance = [-1]*n
        distance[root] = 0
        todo = deque([root])
        while todo:
            u = todo.popleft()
            for v in np.flatnonzero(a[u]):
                if distance[v] < 0:
                    distance[v] = distance[u]+1
                    todo.append(int(v))
        assert min(distance) >= 0
        pair_distances.extend(distance[root+1:])
    return {'nodes': n, 'edges': int(a.sum()//2),
            'degree_counts': dict(sorted(Counter(map(int, a.sum(axis=1))).items())),
            'diameter': max(pair_distances),
            'unordered_pair_distance_counts': dict(sorted(Counter(pair_distances).items())),
            'mean_distance': sum(pair_distances)/len(pair_distances)}


def copy_and_link(a):
    # A separate possible completion: retain G, copy G, link paired copies.
    # It doubles population but also increases every inherited degree by one.
    n = len(a)
    b = np.zeros((2*n, 2*n), dtype=int)
    b[:n, :n] = b[n:, n:] = a
    b[:n, n:] = b[n:, :n] = np.eye(n, dtype=int)
    return b


def report():
    a, b = expand_five(1), expand_five(2)
    assert np.array_equal(a[:5, :5], seed())
    assert np.array_equal(b[:5, :5], seed())
    ma, mb = metrics(a), metrics(b)
    assert ma['nodes'] == mb['nodes'] == 10
    assert ma['edges'] == mb['edges'] == 15
    assert ma['degree_counts'] == mb['degree_counts'] == {3: 10}
    assert ma['diameter'] == 3 and mb['diameter'] == 2
    assert ma['unordered_pair_distance_counts'] == {1: 15, 2: 20, 3: 10}
    assert mb['unordered_pair_distance_counts'] == {1: 15, 2: 30}
    g = seed()
    copies = []
    for n in range(5):
        assert len(g) == 5*2**n
        assert np.all(g.sum(axis=1) == 2+n)
        copies.append({'generation': n, 'nodes': len(g), 'degree': 2+n})
        g = copy_and_link(g)
    return {
        'date': '2026-09-22',
        'kind': 'Unnumbered user-proposal review; main research remains through round 273.',
        'user_clarification': '+5 means additional nodes, not numeric relabeling.',
        'fixed_increment_extension_if_repeated': [5+5*n for n in range(6)],
        'independent_groups_if_all_five_are_retained': [
            {'generation': n, 'group_sizes': [i+5*n for i in range(1, 6)],
             'total': 15+25*n} for n in range(4)],
        'seed': metrics(seed()),
        'same_budget_extensions': [
            {'inner_step': 1, **ma},
            {'inner_step': 2, **mb}],
        'copy_and_link_is_a_different_rule': copies,
        'checks': {
            'original_seed_preserved_in_both': True,
            'equal_new_nodes_and_edges': True,
            'equal_cubic_degree': True,
            'different_distance_distributions': True,
            'copy_and_link_degree_growth_verified': True},
        'scope': 'Arithmetic recursion allows population growth but does not specify adjacency. Two explicit five-to-ten extensions have equal node, edge and degree budgets and different graph distances. A supplied copy-and-link completion doubles population but is not closed under a fixed degree-three cap. No unique spatial geometry, node identities, quantum channels or physical resource source is derived.'
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = report()
    result['runtime'] = {'python': platform.python_version(), 'numpy': np.__version__}
    result['code_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    payload = json.dumps(result, ensure_ascii=False, indent=2)+'\n'
    path = Path(__file__).with_name('five_node_growth_proposal_results.json')
    if args.write_results:
        if path.exists() and path.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing proposal result.')
        path.write_text(payload, encoding='utf-8')
    print(payload)
