"""Round 285: anonymous finite-ball mapping with one fixed anchor.

The mapper sees only local degree, exit/arrival ports and the anchor detector.
Static reversible channels and stable readable memory are explicit inputs.
Finite-radius adaptation of marked-homebase BFS (Chalopin, Das, Kosowski,
OPODIS 2010, Proposition 9); not a claim of a new general mapping algorithm.
"""
import unittest
import numpy as np
from growing_stream_audit import main
from observable_cover_audit import complete_four, cube
from operational_dimension_audit import decorated_torus, distances
from fisher_hanoi_scaling_audit import hanoi


class PortProbe:
    """Simulator internals are not part of the mapper's observation interface."""
    def __init__(self, rows, root=0):
        self._rows = [tuple(row) for row in rows]
        self._position = self._root = root
        self.steps = self.reads = 0
        for v, row in enumerate(self._rows):
            assert v not in row and len(row) == len(set(row))
            assert all(v in self._rows[w] for w in row)

    def degree(self):
        return len(self._rows[self._position])

    def step(self, port):
        old = self._position
        self._position = self._rows[old][port]
        self.steps += 1
        return self._rows[self._position].index(old)

    def at_anchor(self):
        self.reads += 1
        return self._position == self._root


def reverse_tape(tape):
    return tuple((q, p) for p, q in reversed(tape))


def follow(probe, tape):
    """Replay a known valid root path or a recorded actual return path."""
    for p, q in tape:
        if probe.step(p) != q:
            raise RuntimeError('A stored route changed: the static-port model failed.')


def same_endpoint(probe, alpha, beta, read_anchor=None):
    """Start and finish at root, even when the attempted inverse route fails."""
    follow(probe, alpha)
    actual = []
    valid = True
    for p, expected in reverse_tape(beta):
        if p >= probe.degree():
            valid = False
            break
        q = probe.step(p)
        actual.append((p, q))
        if q != expected:
            valid = False
            break
    answer = bool((read_anchor or probe.at_anchor)()) if valid else False
    follow(probe, reverse_tape(actual))
    follow(probe, reverse_tape(alpha))
    return answer


def extend(probe, tape, port):
    follow(probe, tape)
    incoming = probe.step(port)
    target_degree = probe.degree()
    extension = tape + ((port, incoming),)
    follow(probe, reverse_tape(extension))
    return extension, target_degree


def map_ball(probe, radius, read_anchor=None):
    """Breadth-first discovery; no vertex IDs, coordinates or graph size oracle."""
    if radius < 0:
        raise ValueError('radius must be nonnegative')
    representatives, degrees = [()], [probe.degree()]
    adjacency = []
    queries = probes = 0
    cursor = 0
    while cursor < len(representatives):
        tape = representatives[cursor]
        local = []
        for port in range(degrees[cursor]):
            extended, degree = extend(probe, tape, port)
            probes += 1
            target = None
            for i, candidate in enumerate(representatives):
                queries += 1
                if same_endpoint(probe, extended, candidate, read_anchor):
                    target = i
                    break
            if target is None and len(tape) < radius:
                target = len(representatives)
                representatives.append(extended)
                degrees.append(degree)
            local.append(target)  # None is an edge leaving the requested ball.
        adjacency.append(local)
        cursor += 1
    return {'representatives':representatives, 'adjacency':adjacency,
            'root_distances':[len(t) for t in representatives],
            'degrees':degrees, 'comparisons':queries, 'edge_probes':probes,
            'steps':probe.steps, 'marker_reads':probe.reads,
            'stored_port_pairs':sum(map(len,representatives))}


def endpoint(rows, root, tape):
    """Ground-truth oracle: used only outside the mapping algorithm."""
    v = root
    for p, q in tape:
        w = rows[v][p]
        assert rows[w][q] == v
        v = w
    return v


def certify_map(rows, root, radius, value):
    ids = [endpoint(rows,root,t) for t in value['representatives']]
    metric = distances(rows,root)
    expected = {i for i,d in enumerate(metric) if 0 <= d <= radius}
    if set(ids) != expected or len(ids) != len(expected):
        return False
    for i, row in enumerate(value['adjacency']):
        actual = [None if j is None else ids[j] for j in row]
        want = [w if w in expected else None for w in rows[ids[i]]]
        if actual != want or value['root_distances'][i] != int(metric[ids[i]]):
            return False
    return True


def bounds(radius, count, degree=3):
    comparisons = degree*count*count
    return {'comparisons_upper':comparisons,
            'steps_upper':2*(radius+1)*degree*count+(4*radius+2)*comparisons,
            'representative_pairs_upper':count*radius,
            'adjacency_entries_upper':degree*count}


def all_tapes(rows, root, radius):
    tapes = [()]
    layer = [(root,())]
    for _ in range(radius):
        new = []
        for v,tape in layer:
            for p,w in enumerate(rows[v]):
                pair = (p,rows[w].index(v))
                new.append((w,tape+(pair,)))
        layer = new
        tapes.extend(t for _,t in layer)
    return tapes


def shuffled(rows, seed):
    rng = np.random.default_rng(seed)
    return [tuple(int(x) for x in rng.permutation(row)) for row in rows]


def case(name, rows, radius, root=0):
    probe = PortProbe(rows,root)
    value = map_ball(probe,radius)
    assert probe._position == root and certify_map(rows,root,radius,value)
    n = len(value['representatives'])
    return {'name':name,'radius':radius,'nodes_in_ball':n,
            'induced_edges':sum(j is not None for row in value['adjacency'] for j in row)//2,
            'edge_probes':value['edge_probes'],'comparisons':value['comparisons'],
            'steps':value['steps'],'marker_reads':value['marker_reads'],
            'stored_port_pairs':value['stored_port_pairs'],
            'physical_anchor_tokens':1,'correct_against_hidden_oracle':True,
            **bounds(radius,n)}


def report():
    torus,_ = decorated_torus(2,3)
    hanoi_rows,_ = hanoi(3)
    return {'round':285,
            'scope':'Conditional classical mapping of a static reversible port graph; no dimension selection or quantum/spacetime derivation.',
            'cases':[case('K4',complete_four(),1),
                     case('cube',cube(),1),case('cube',cube(),2),case('cube',cube(),3),
                     case('Hanoi_27',hanoi_rows,3),
                     case('decorated_2D_torus_36',torus,3)],
            'exit_only_false_positive':{
                'graph':[[1],[0,2],[1]],'root':1,'alpha':[[0,0]],'beta':[[1,0]],
                'distinct_endpoints':[0,2],'root_detector_after_exit_only_inverse':True,
                'expected_arrival_port':1,'actual_arrival_port':0,
                'paired_port_protocol_correctly_rejects':True},
            'conditions':['Static undirected simple graph, degree at most 3',
                          'Locally numbered exit and readable arrival ports, reversible traversal',
                          'One unique persistent root token with local contact detection',
                          'Reliable route memory, chosen finite radius, no unknown-data disturbance claim'],
            'resource_claim':'Upper bounds for this explicit BFS protocol; no optimality or physical energy bound.'}


class Audit(unittest.TestCase):
    def test_01_pair_inverse(self):
        rows = shuffled(cube(),12)
        for root in range(8):
            for tape in all_tapes(rows,root,2):
                probe = PortProbe(rows,root)
                follow(probe,tape)
                follow(probe,reverse_tape(tape))
                self.assertEqual(probe._position,root)

    def test_02_endpoint_relation_all_pairs(self):
        for rows in (complete_four(),cube(),[[1],[0,2],[1]]):
            for root in range(len(rows)):
                tapes = all_tapes(rows,root,2)
                probe = PortProbe(rows,root)
                for alpha in tapes:
                    for beta in tapes:
                        before = probe.steps
                        self.assertEqual(same_endpoint(probe,alpha,beta),
                                         endpoint(rows,root,alpha)==endpoint(rows,root,beta))
                        self.assertEqual(probe._position,root)
                        self.assertLessEqual(probe.steps-before,2*(len(alpha)+len(beta)))

    def test_03_exit_only_counterexample(self):
        rows = [[1],[0,2],[1]]
        alpha,beta = ((0,0),),((1,0),)
        probe = PortProbe(rows,1)
        follow(probe,alpha)
        arrived = probe.step(0)
        self.assertTrue(probe.at_anchor())
        self.assertNotEqual(arrived,1)
        self.assertFalse(same_endpoint(PortProbe(rows,1),alpha,beta))

    def test_04_failure_restoration(self):
        # The inverse requires a missing exit at a leaf.
        rows = [[1,2],[0],[0,3],[2]]
        probe = PortProbe(rows)
        self.assertFalse(same_endpoint(probe,((0,0),),((1,0),(1,0),(0,1))))
        self.assertEqual(probe._position,0)
        # The first inverse step is traversable, but has a wrong arrival port.
        probe = PortProbe([[1],[0,2],[1]],1)
        self.assertFalse(same_endpoint(probe,((0,0),),((1,0),)))
        self.assertEqual(probe._position,1)

    def test_05_maps_and_radius_zero(self):
        torus,_ = decorated_torus(2,3)
        hr,_ = hanoi(2)
        for rows in (complete_four(),cube(),torus,hr,[[]]):
            for radius in range(4):
                probe = PortProbe(rows)
                value = map_ball(probe,radius)
                self.assertTrue(certify_map(rows,0,radius,value))
                self.assertEqual(probe._position,0)

    def test_06_arbitrary_ports_and_roots(self):
        for seed in range(4):
            rows = shuffled(cube(),seed)
            for root in (0,3,7):
                value = map_ball(PortProbe(rows,root),2)
                self.assertTrue(certify_map(rows,root,2,value))

    def test_07_resource_bounds(self):
        for row in report()['cases']:
            self.assertLessEqual(row['comparisons'],row['comparisons_upper'])
            self.assertLessEqual(row['steps'],row['steps_upper'])
            self.assertLessEqual(row['stored_port_pairs'],row['representative_pairs_upper'])
            self.assertLessEqual(row['marker_reads'],row['comparisons'])

    def test_08_boundary_edges_and_private_interface(self):
        # A wrapper deliberately exposes no simulator coordinates or vertex IDs.
        simulator = PortProbe(complete_four())
        class LocalOnly:
            degree = simulator.degree
            step = simulator.step
            at_anchor = simulator.at_anchor
            @property
            def steps(self): return simulator.steps
            @property
            def reads(self): return simulator.reads
        value = map_ball(LocalOnly(),1)
        self.assertTrue(certify_map(complete_four(),0,1,value))
        self.assertEqual(sum(j is not None for row in value['adjacency'] for j in row),12)
        self.assertEqual(sorted(value['root_distances']),[0,1,1,1])


if __name__ == '__main__':
    main(__name__,'anchored_position_audit',report)
