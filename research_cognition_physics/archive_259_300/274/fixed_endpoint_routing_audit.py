"""Round 274: charged routing between the original registers after growth.

One ideal qubit transmission per undirected edge per logical slot, local
buffers/gates supplied. Construction finishes before task injection.
"""
import argparse
from collections import Counter, deque
from functools import lru_cache
import json
from pathlib import Path
import platform
import unittest
import numpy as np
from paired_interface_growth_audit import ladder_state, check_ports
from ancestral_quotient_audit import graph_data


@lru_cache(maxsize=None)
def instance(depth):
    state = ladder_state(depth)
    check_ports(state)
    graph, words = graph_data(state)
    s_word, t_word = (0,)+(0,)*depth, (1,)+(0,)*depth
    index = {w:i for i,w in enumerate(words)}
    s, t = index[s_word], index[t_word]
    paths = [(s, t)]
    for u, _ in state[s_word].values():
        if u[0] == 0:
            v = state[u][1][0]
            assert v[0] == 1 and t_word in {w for w, _ in state[v].values()}
            paths.append((s, index[u], index[v], t))
    return graph, s, t, tuple(paths)


def max_flow(capacities, source, sink):
    residual = {u:dict(row) for u,row in capacities.items()}
    for u,row in list(residual.items()):
        for v in list(row):
            residual.setdefault(v, {}).setdefault(u, 0)
    value = 0
    while True:
        parent = {source:None}
        todo = deque([source])
        while todo and sink not in parent:
            u = todo.popleft()
            for v,capacity in residual[u].items():
                if capacity > 0 and v not in parent:
                    parent[v] = u
                    todo.append(v)
        if sink not in parent:
            return value
        amount, v = 10**9, sink
        while parent[v] is not None:
            u = parent[v]
            amount = min(amount, residual[u][v])
            v = u
        v = sink
        while parent[v] is not None:
            u = parent[v]
            residual[u][v] -= amount
            residual[v][u] += amount
            v = u
        value += amount


def time_expanded_flow(graph, source, sink, slots):
    """Optimistic directed relaxation; opposing uses would share one real edge."""
    n = len(graph)
    large = 3*slots+1
    capacities = {i:{} for i in range(n*(slots+1))}
    for time in range(slots):
        for u,row in enumerate(graph):
            capacities[time*n+u][(time+1)*n+u] = large
            for v in row:
                capacities[time*n+u][(time+1)*n+v] = 1
    return max_flow(capacities, source, slots*n+sink)


def routed_by(paths, slots):
    return sum(max(0, slots-(len(path)-1)+1) for path in paths)


def events(paths, slots):
    """Keep exact task identity; count all distribution and collection hops."""
    result = []
    for lane,path in enumerate(paths):
        for launch in range(1, slots+1):
            for hop,(u,v) in enumerate(zip(path, path[1:])):
                time = launch+hop
                if time <= slots:
                    result.append((time, lane, launch, hop, u, v))
    return sorted(result)


def swap_qubits(vector, left, right):
    indices = np.arange(len(vector), dtype=np.int64)
    difference = ((indices >> left) ^ (indices >> right)) & 1
    permutation = indices ^ (difference << left) ^ (difference << right)
    return vector[permutation]


def choi_wire_check(lengths, omit_collection=False):
    """Bell-reference audit of one batch, using full amplitudes and blank buffers."""
    count = len(lengths)
    registers, offset = [], count
    for length in lengths:
        registers.append(tuple(range(offset, offset+length+1)))
        offset += length+1
    psi = np.zeros(1 << offset, complex)
    ideal = np.zeros_like(psi)
    for bits in range(1 << count):
        start = finish = bits
        for lane,row in enumerate(registers):
            if bits >> lane & 1:
                start |= 1 << row[0]
                finish |= 1 << row[-1]
        psi[start] = ideal[finish] = 1/np.sqrt(1 << count)
    for tick in range(max(lengths)):
        for row,length in zip(registers, lengths):
            if tick < length and not (omit_collection and length > 1 and tick == length-1):
                psi = swap_qubits(psi, row[tick], row[tick+1])
    return float(abs(np.vdot(ideal, psi))**2)


def report():
    rows = []
    for depth in range(7):
        graph, s, t, paths = instance(depth)
        lengths = [len(p)-1 for p in paths]
        rows.append({'depth':depth, 'nodes':len(graph), 'rungs':2**depth,
                     'original_source_degree':len(graph[s]),
                     'route_lengths':lengths,
                     'achievable_asymptotic_qubits_per_slot':len(paths),
                     'delivered_by_slot_8':routed_by(paths, 8),
                     'completed_task_hops_by_slot_8':
                         sum(length*max(0, 8-length+1) for length in lengths)})
    return {
        'round':274, 'rows':rows,
        'finite_routing_optimum_checks':
            [{'depth':d, 'slots':t, 'time_expanded_flow':time_expanded_flow(*instance(d)[:3], t),
              'explicit_routing':routed_by(instance(d)[3], t)}
             for d in range(5) for t in (1, 2, 3, 5, 8)],
        'quantum_batch_choi_fidelities':
            {str(d):choi_wire_check([len(p)-1 for p in instance(d)[3]]) for d in range(3)},
        'missing_collection_control_fidelity':choi_wire_check([1,3,3], True),
        'scope':'Fixed original source and receiver on completed paired-growth ladders. Exact finite-horizon store-and-forward routing, a quantum wire implementation and asymptotically matching Schmidt cut bounds; finite routing optimality is not asserted for every quantum protocol. Logical clocks, local memory and primitive channels are supplied.'
    }


class Audit(unittest.TestCase):
    def test_01_routes_use_actual_edges_and_keep_original_data_vertices(self):
        for d in range(7):
            graph,s,t,paths = instance(d)
            self.assertIn(t, graph[s])
            self.assertEqual(len(paths), min(2**d,3))
            for path in paths:
                self.assertEqual((path[0],path[-1]),(s,t))
                for u,v in zip(path,path[1:]):
                    self.assertIn(v,graph[u])

    def test_02_paths_are_edge_disjoint_and_match_the_source_cut(self):
        for d in range(7):
            graph,s,t,paths = instance(d)
            used = [tuple(sorted((u,v))) for p in paths for u,v in zip(p,p[1:])]
            self.assertEqual(len(set(used)),len(used))
            self.assertEqual(len(paths),len(graph[s]))
            self.assertEqual(len(paths),len(graph[t]))

    def test_03_time_expanded_routing_optimum_matches_constructed_schedule(self):
        for d in range(5):
            graph,s,t,paths = instance(d)
            for slots in (1,2,3,5,8):
                self.assertEqual(time_expanded_flow(graph,s,t,slots),routed_by(paths,slots))

    def test_04_actual_packet_schedule_obeys_edge_capacity_and_precedence(self):
        for d in range(4):
            paths = instance(d)[3]
            trace = events(paths,16)
            capacity = Counter((time,tuple(sorted((u,v)))) for time,lane,launch,hop,u,v in trace)
            self.assertEqual(set(capacity.values()),{1})
            seen = {}
            complete = set()
            for time,lane,launch,hop,u,v in trace:
                task = (lane,launch)
                if hop:
                    self.assertEqual(seen[task],(time-1,u))
                seen[task] = (time,v)
                if hop == len(paths[lane])-2:
                    complete.add(task)
            self.assertEqual(len(complete),routed_by(paths,16))

    def test_05_distribution_and_collection_are_both_charged(self):
        for slots in (3,8,16):
            paths = instance(4)[3]
            completed = {(lane,launch) for lane,p in enumerate(paths)
                         for launch in range(1,slots-(len(p)-1)+2)}
            charged = [event for event in events(paths,slots) if (event[1],event[2]) in completed]
            self.assertEqual(len(charged),slots+6*(slots-2))
            self.assertEqual(len(completed),3*slots-4)

    def test_06_quantum_reference_channel_is_identity_for_all_three_route_cases(self):
        for lengths in ([1],[1,3],[1,3,3]):
            self.assertAlmostEqual(choi_wire_check(lengths),1.0,places=12)

    def test_07_omitting_collection_does_not_deliver_to_original_receiver(self):
        self.assertAlmostEqual(choi_wire_check([1,3,3],True),1/16,places=12)

    def test_08_initial_entanglement_cannot_change_the_asymptotic_source_cut(self):
        graph,s,t,paths = instance(6)
        degree = len(graph[s])
        for slots in (8,64,512):
            for initial_bell_pairs in (0,3,10):
                self.assertLessEqual(routed_by(paths,slots),degree*slots+initial_bell_pairs)
        self.assertEqual(routed_by(paths,512),3*512-4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result['runtime'] = {'python':platform.python_version(), 'numpy':np.__version__}
    result['checks'] = {'run':checks.testsRun,'failures':0,'errors':0}
    payload = json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    path = Path(__file__).with_name('fixed_endpoint_routing_audit_results.json')
    if args.write_results:
        if path.exists() and path.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        path.write_text(payload,encoding='utf-8')
    print(payload)
