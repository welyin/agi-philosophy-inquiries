"""Round 287: physically grounded overlap matching between anonymous maps.

Reuses marked-homebase mapping from round 285 / Chalopin-Das-Kosowski (2010).
Unique local markers, transferable port names and stable reversible channels
are explicit inputs; no graph coordinates are exposed to the matcher.
"""
import math
import unittest
from growing_stream_audit import main
from anchored_position_audit import (PortProbe, map_ball, follow, reverse_tape,
                                    same_endpoint, endpoint, shuffled)
from observable_cover_audit import cube


class TaggedProbe(PortProbe):
    def __init__(self, rows, tags, own):
        self._tags = dict(tags)
        self._own = own
        super().__init__(rows,self._tags[own])

    def detect(self, tag):
        self.reads += 1
        return self._position == self._tags[tag]

    def at_anchor(self):
        return self.detect(self._own)


def cycle(n):
    assert n >= 3
    return [((v+1)%n,(v-1)%n) for v in range(n)]


def discover_bridge(probe, survey, target):
    """A bounded survey only; None means not found within that survey."""
    for tape in survey['representatives']:
        follow(probe,tape)
        found = probe.detect(target)
        follow(probe,reverse_tape(tape))
        if found:
            return tape
    return None


def record_bridge(probe, exits, target):
    """Walk a supplied local control sequence, check target, and return home."""
    actual = []
    for port in exits:
        if not 0 <= port < probe.degree():
            follow(probe,reverse_tape(actual))
            return None
        incoming = probe.step(port)
        actual.append((port,incoming))
    valid = probe.detect(target)
    follow(probe,reverse_tape(actual))
    return tuple(actual) if valid else None


def match_overlap(probe, chart_a, chart_b, bridge, read_anchor=None):
    if bridge is None:
        raise ValueError('A verified bridge is required; map shape is insufficient.')
    start_steps,start_reads = probe.steps,probe.reads
    matches,queries = {},0
    for i,alpha in enumerate(chart_a['representatives']):
        for j,beta in enumerate(chart_b['representatives']):
            queries += 1
            if same_endpoint(probe,alpha,bridge+beta,read_anchor):
                matches[i] = j
                break
    return {'matches':matches,'comparisons':queries,
            'steps':probe.steps-start_steps,'marker_reads':probe.reads-start_reads}


def overlap_oracle(rows,roots,charts):
    left = [endpoint(rows,roots[0],t) for t in charts[0]['representatives']]
    right = {endpoint(rows,roots[1],t):j for j,t in enumerate(charts[1]['representatives'])}
    return {i:right[v] for i,v in enumerate(left) if v in right}


def packet_bits(chart,radius):
    """One fixed-width encoding; no claim of optimal compression."""
    n = len(chart['representatives'])
    length_width = max(1,math.ceil(math.log2(radius+1)))
    index_width = max(1,math.ceil(math.log2(n+1)))
    pairs = sum(map(len,chart['representatives']))
    # Two 2-bit ports per pair, a length and degree per vertex, and adjacency.
    return {'route_bits':4*pairs+n*length_width,
            'map_bits':4*pairs+n*(length_width+2)+sum(chart['degrees'])*index_width,
            'header_included':False}


def pair_case(name,rows,roots,radii,survey_radius):
    tags = {'A':roots[0],'B':roots[1]}
    probes = [TaggedProbe(rows,tags,t) for t in ('A','B')]
    charts = [map_ball(p,r) for p,r in zip(probes,radii)]
    setup_steps = sum(p.steps for p in probes)
    setup_reads = sum(p.reads for p in probes)
    a = probes[0]
    before = a.steps
    survey = charts[0] if survey_radius==radii[0] else map_ball(a,survey_radius)
    survey_steps = a.steps-before
    before,reads_before = a.steps,a.reads
    bridge = discover_bridge(a,survey,'B')
    bridge_steps,bridge_reads = a.steps-before,a.reads-reads_before
    if bridge is None:
        return {'name':name,'bridge_found':False,'survey_radius':survey_radius}
    value = match_overlap(a,*charts,bridge)
    assert value['matches'] == overlap_oracle(rows,roots,charts)
    assert a._position == roots[0]
    na,nb = [len(c['representatives']) for c in charts]
    packet = packet_bits(charts[1],radii[1])
    return {'name':name,'radii':list(radii),'chart_sizes':[na,nb],
            'overlap_size':len(value['matches']),'bridge_length':len(bridge),
            'map_setup_steps':setup_steps,'map_setup_marker_reads':setup_reads,
            'additional_survey_steps':survey_steps,
            'bridge_search_steps':bridge_steps,'bridge_search_reads':bridge_reads,
            'bridge_search_steps_upper':2*survey_radius*len(survey['representatives']),
            'match_comparisons':value['comparisons'],'match_steps':value['steps'],
            'match_marker_reads':value['marker_reads'],
            'match_comparisons_upper':na*nb,
            'match_steps_upper':2*(sum(radii)+len(bridge))*na*nb,
            'unique_marker_tokens':2,'common_local_port_dialect':True,
            'map_packet_payload_bits':packet['map_bits'],
            'payload_bit_edge_uses':packet['map_bits']*len(bridge),
            'packet_header_and_control_traffic_included':False,
            'truth_check_passed':True}


def identical_map_counterexample():
    rows = cycle(12)
    charts = [map_ball(PortProbe(rows,root),1) for root in (0,2,5)]
    comparable = lambda c:(c['representatives'],c['adjacency'],c['degrees'])
    assert comparable(charts[0]) == comparable(charts[1]) == comparable(charts[2])
    overlaps = [len(overlap_oracle(rows,(0,r),(charts[0],charts[j])))
                for j,r in ((1,2),(2,5))]
    return {'graph':'C12','radius':1,'roots':[0,2,5],
            'rooted_port_maps_identical':True,'overlaps_with_first':overlaps,
            'naive_root_to_root_identity_correct':False}


def duplicated_marker_counterexample():
    rows = cycle(6)
    class DuplicateProbe(PortProbe):
        def at_anchor(self):
            self.reads += 1
            return self._position in (0,3)
    result = map_ball(DuplicateProbe(rows),3)
    return {'graph':'C6','indistinguishable_markers_at':[0,3],
            'physical_nodes_in_requested_ball':6,
            'reported_nodes':len(result['representatives']),
            'reported_adjacency':result['adjacency'],
            'reason':'A duplicated marker violates the uniqueness input and yields a quotient map.'}


def report():
    return {'round':287,
            'scope':'Conditional finite classical chart overlap from physically verified bridges; no coordinates, dimension selection or spacetime generation.',
            'identical_map_counterexample':identical_map_counterexample(),
            'cases':[pair_case('C12_overlap',cycle(12),(0,2),(1,1),2),
                     pair_case('C12_disjoint',cycle(12),(0,5),(1,1),5),
                     pair_case('cube_partial',cube(),(0,1),(1,1),1),
                     pair_case('cube_larger_overlap',cube(),(0,3),(2,2),2),
                     pair_case('cube_shuffled_ports',shuffled(cube(),17),(0,3),(2,2),2)],
            'duplicate_marker_counterexample':duplicated_marker_counterexample(),
            'resource_scope':'Packet payload counted separately from marker identity, framing, control traffic, bandwidth and energy; no free global broadcast assumed.'}


class Audit(unittest.TestCase):
    def test_01_same_shape_not_same_overlap(self):
        self.assertEqual(identical_map_counterexample()['overlaps_with_first'],[1,0])

    def test_02_all_cycle_root_pairs(self):
        rows = cycle(8)
        for a in range(8):
            for b in range(8):
                tags = {'A':a,'B':b}
                pa,pb = TaggedProbe(rows,tags,'A'),TaggedProbe(rows,tags,'B')
                ca,cb = map_ball(pa,2),map_ball(pb,2)
                w = record_bridge(pa,[0]*((b-a)%8),'B')
                out = match_overlap(pa,ca,cb,w)
                self.assertEqual(out['matches'],overlap_oracle(rows,(a,b),(ca,cb)))
                self.assertEqual(pa._position,a)

    def test_03_port_renaming(self):
        for seed in range(6):
            row = pair_case('renamed',shuffled(cube(),seed),(0,3),(2,2),2)
            self.assertEqual(row['overlap_size'],6)

    def test_04_bounded_discovery_is_not_absence(self):
        p = TaggedProbe(cycle(12),{'A':0,'B':5},'A')
        chart = map_ball(p,1)
        self.assertIsNone(discover_bridge(p,chart,'B'))
        with self.assertRaises(ValueError):
            match_overlap(p,chart,chart,None)
        self.assertIsNotNone(record_bridge(p,[0]*5,'B'))
        self.assertEqual(p._position,0)

    def test_05_different_bridges_same_overlap(self):
        rows = cycle(8)
        p = TaggedProbe(rows,{'A':0,'B':2},'A')
        ca,cb = map_ball(p,2),map_ball(PortProbe(rows,2),2)
        short = record_bridge(p,[0]*2,'B')
        long = record_bridge(p,[1]*6,'B')
        self.assertNotEqual(short,long)
        self.assertEqual(match_overlap(p,ca,cb,short)['matches'],
                         match_overlap(p,ca,cb,long)['matches'])
        self.assertIsNone(record_bridge(p,[0],'B'))

    def test_06_duplicate_marker_is_not_safe(self):
        row = duplicated_marker_counterexample()
        self.assertEqual((row['physical_nodes_in_requested_ball'],row['reported_nodes']),(6,3))

    def test_07_resources(self):
        for row in report()['cases']:
            self.assertLessEqual(row['match_comparisons'],row['match_comparisons_upper'])
            self.assertLessEqual(row['match_steps'],row['match_steps_upper'])
            self.assertLessEqual(row['bridge_search_steps'],row['bridge_search_steps_upper'])
            self.assertEqual(row['payload_bit_edge_uses'],
                             row['map_packet_payload_bits']*row['bridge_length'])

    def test_08_overlap_edges_and_local_interface(self):
        rows = cube()
        p = TaggedProbe(rows,{'A':0,'B':1},'A')
        ca,cb = map_ball(p,2),map_ball(PortProbe(rows,1),2)
        w = discover_bridge(p,ca,'B')
        class LocalOnly:
            degree = p.degree
            step = p.step
            at_anchor = p.at_anchor
            @property
            def steps(self): return p.steps
            @property
            def reads(self): return p.reads
        matches = match_overlap(LocalOnly(),ca,cb,w)['matches']
        self.assertEqual(matches,overlap_oracle(rows,(0,1),(ca,cb)))
        for i,j in matches.items():
            for port,neighbor in enumerate(ca['adjacency'][i]):
                if neighbor in matches:
                    self.assertEqual(cb['adjacency'][j][port],matches[neighbor])


if __name__ == '__main__':
    main(__name__,'chart_overlap_audit',report)

