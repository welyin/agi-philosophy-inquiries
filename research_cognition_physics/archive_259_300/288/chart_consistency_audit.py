"""Round 288: chart synchronization, cycle tests and physical grounding.

Permutation synchronization follows Pachauri-Kondor-Singh (NIPS 2013).
Consistency of labels is distinguished from correct physical correspondence.
"""
from collections import deque
from itertools import permutations, product
import math
import unittest
import numpy as np
from growing_stream_audit import main
from chart_overlap_audit import (TaggedProbe, cycle, discover_bridge,
                                match_overlap, overlap_oracle)
from anchored_position_audit import map_ball, certify_map
from observable_cover_audit import cube
from noisy_anchor_audit import ball_upper, majority_error, exponential_bound


def compose(after,before):
    return tuple(after[x] for x in before)


def inverse(p):
    out = [0]*len(p)
    for i,j in enumerate(p):
        out[j] = i
    return tuple(out)


def edge_pair(edges,i,j,p):
    edges[i,j] = tuple(p)
    edges[j,i] = inverse(p)


def synchronize(count,edges):
    """Full bijections on a common finite set; not partial-overlap completion."""
    size = len(next(iter(edges.values())))
    identity = tuple(range(size))
    for (i,j),p in edges.items():
        assert sorted(p)==list(identity) and edges[j,i]==inverse(p)
    frames = {0:identity}
    queue = deque([0])
    tree_edges = []
    while queue:
        i = queue.popleft()
        for j in sorted(v for u,v in edges if u==i):
            if j not in frames:
                frames[j] = compose(edges[i,j],frames[i])
                tree_edges.append((i,j))
                queue.append(j)
    if len(frames)!=count:
        raise ValueError('The comparison graph must be connected.')
    disagreements = []
    for (i,j),p in sorted(edges.items()):
        if i<j:
            want = compose(frames[j],inverse(frames[i]))
            if p!=want:
                disagreements.append([i,j])
    return {'frames':frames,'tree_edges':tree_edges,
            'inconsistent_edges':disagreements,'consistent':not disagreements,
            'independent_cycles':len(edges)//2-count+1}


def permutation_matrix(p):
    out = np.zeros((len(p),len(p)))
    out[p,np.arange(len(p))] = 1.
    return out


def block_matrix(count,edges):
    n = len(next(iter(edges.values())))
    return np.block([[np.eye(n) if i==j else permutation_matrix(edges[j,i])
                      for j in range(count)] for i in range(count)])


def partial_triple_failures(ab,bc,ac):
    """Check only points for which all three chart representations exist."""
    checked = failed = 0
    for x,y in ab.items():
        if y in bc and x in ac:
            checked += 1
            failed += bc[y]!=ac[x]
    return checked,failed


class NoisyTags(TaggedProbe):
    def __init__(self,rows,tags,own,p,m,seed):
        majority_error(p,m)
        super().__init__(rows,tags,own)
        self.p,self.m = p,m
        self.rng = np.random.default_rng(seed)
        self.blocks = self.errors = 0

    def detect(self,tag):
        truth = super().detect(tag)
        observations = [truth]+[super(NoisyTags,self).detect(tag) for _ in range(self.m-1)]
        flips = self.rng.random(self.m)<self.p
        answer = bool(np.count_nonzero(np.logical_xor(observations,flips))>self.m//2)
        self.blocks += 1
        self.errors += answer!=truth
        return answer


def triple_run(p=0.,m=1,seed=0):
    rows = cube()
    roots = (0,1,3)
    tags = dict(enumerate(roots))
    probes = [NoisyTags(rows,tags,i,p,m,seed+i*1000) for i in range(3)]
    charts = [map_ball(probe,2) for probe in probes]
    map_errors = sum(not certify_map(rows,roots[i],2,chart)
                     for i,chart in enumerate(charts))
    pair_maps, bridge_lengths = {},[]
    physical_errors = 0
    for i,j in ((0,1),(1,2),(0,2)):
        bridge = discover_bridge(probes[i],charts[i],j)
        if bridge is None:
            return {'correct':False,'bridge_not_found':True}
        value = match_overlap(probes[i],charts[i],charts[j],bridge)
        pair_maps[i,j] = value['matches']
        physical_errors += value['matches']!=overlap_oracle(rows,(roots[i],roots[j]),
                                                           (charts[i],charts[j]))
        bridge_lengths.append(len(bridge))
    checked,failed = partial_triple_failures(pair_maps[0,1],pair_maps[1,2],pair_maps[0,2])
    return {'correct':physical_errors==0 and map_errors==0,'bridge_not_found':False,
            'chart_sizes':[len(c['representatives']) for c in charts],
            'pair_overlap_sizes':[len(pair_maps[e]) for e in ((0,1),(1,2),(0,2))],
            'bridge_lengths':bridge_lengths,'triple_positions_checked':checked,
            'triple_failures':failed,'physical_pair_errors':physical_errors,
            'physical_map_errors':map_errors,
            'edge_traversals':sum(probe.steps for probe in probes),
            'marker_read_blocks':sum(probe.blocks for probe in probes),
            'raw_marker_reads':sum(probe.reads for probe in probes),
            'wrong_read_blocks':sum(probe.errors for probe in probes),
            'markers':3,'p':p,'repetitions':m,
            'scope':'Map construction, bridge search and pair matching; communication overhead is separate.'}


def all_consistent_but_wrong():
    rows = cycle(6)
    roots = (0,1,2)
    charts = [map_ball(TaggedProbe(rows,dict(enumerate(roots)),i),3) for i in range(3)]
    assert charts[0]['adjacency']==charts[1]['adjacency']==charts[2]['adjacency']
    edges = {}
    for i,j in ((0,1),(1,2),(0,2)):
        edge_pair(edges,i,j,range(6))  # Same rooted map index; no physical bridge.
    sync = synchronize(3,edges)
    wrong = []
    for i,j in ((0,1),(1,2),(0,2)):
        expected = overlap_oracle(rows,(roots[i],roots[j]),(charts[i],charts[j]))
        wrong.append(sum(expected[x]!=edges[i,j][x] for x in range(6)))
    vals = np.linalg.eigvalsh(block_matrix(3,edges))
    return {'graph':'C6','origins':list(roots),'all_rooted_maps_identical':True,
            'all_cycle_tests_pass':sync['consistent'],
            'wrong_physical_matches_per_pair':wrong,
            'block_eigenvalues_rounded':np.round(vals,12).tolist(),
            'block_psd':bool(vals.min()>-1e-10),'block_rank':int(np.count_nonzero(vals>1e-10))}


def sparse_cycle_counterexample():
    edges = {}
    for i,j in ((0,1),(1,2),(2,3)):
        edge_pair(edges,i,j,range(3))
    edge_pair(edges,3,0,(1,2,0))
    value = synchronize(4,edges)
    around = tuple(range(3))
    for e in ((0,1),(1,2),(2,3),(3,0)):
        around = compose(edges[e],around)
    return {'comparison_graph':'C4','triangle_checks_available':0,
            'all_edge_inverse_checks_pass':True,'whole_cycle_permutation':list(around),
            'independent_cycle_count':value['independent_cycles'],
            'global_consistency':value['consistent']}


def confidence_case(p=.1,delta=.01):
    size = ball_upper(2)
    # Three map builds + three bounded bridge searches + three pair matchings.
    q = 3*(3*size**2)+3*size+3*size**2
    m = math.ceil(math.log(q/delta)/(-math.log(2*math.sqrt(p*(1-p)))))
    m += m%2==0
    return {'radius':2,'degree_upper':3,'ball_size_upper':size,
            'all_phases_queries_before_first_error_upper':q,
            'p':p,'target_delta':delta,'repetitions':m,
            'failure_upper_exact_union':q*majority_error(p,m),
            'failure_upper_exponential':q*exponential_bound(p,m),
            'good_event_raw_reads_upper':q*m,
            'conditions':'All three true origins lie in each relevant search ball; fresh independent marker read flips, unique tags, exact ports, stable graph and memories.'}


def report():
    risk = confidence_case()
    runs = [triple_run(.1,risk['repetitions'],seed) for seed in range(60)]
    return {'round':288,
            'scope':'Classical map identity and permutation consistency; no Lorentzian geometry, curvature or gauge dynamics derived.',
            'grounded_partial_triple':triple_run(),
            'coherent_wrong_matching':all_consistent_but_wrong(),
            'sparse_cycle_counterexample':sparse_cycle_counterexample(),
            'whole_protocol_confidence':risk,
            'independent_noise_runs':{'trials':len(runs),
                                     'correct_maps_and_matches':sum(r['correct'] for r in runs),
                                     'seed_start':0,'repetitions':risk['repetitions'],
                                     'purpose':'Implementation regression; confidence comes from the analytic bound.'},
            'scope_limit':'Full permutation synchronization assumes equal identified domains. Partial charts use only actual common points; missing overlap is not an identity or zero map.'}


class Audit(unittest.TestCase):
    def test_01_partial_grounded_cocycle(self):
        result = triple_run()
        self.assertTrue(result['correct'])
        self.assertEqual(result['pair_overlap_sizes'],[6,6,6])
        self.assertEqual((result['triple_positions_checked'],result['triple_failures']),(5,0))

    def test_02_exhaustive_three_point_permutations(self):
        perms = list(permutations(range(3)))
        consistent = 0
        for ab,bc,ac in product(perms,repeat=3):
            edges = {}
            for i,j,p in ((0,1,ab),(1,2,bc),(0,2,ac)):
                edge_pair(edges,i,j,p)
            expected = compose(bc,ab)==ac
            self.assertEqual(synchronize(3,edges)['consistent'],expected)
            self.assertEqual(np.linalg.eigvalsh(block_matrix(3,edges)).min()>-1e-9,expected)
            consistent += expected
        self.assertEqual(consistent,36)

    def test_03_tree_reconstruction_and_label_changes(self):
        frames = [(0,1,2,3),(2,0,3,1),(3,2,0,1),(1,3,2,0)]
        edges = {}
        for i,j in ((0,1),(1,2),(2,3),(3,0),(0,2)):
            edge_pair(edges,i,j,compose(frames[j],inverse(frames[i])))
        self.assertTrue(synchronize(4,edges)['consistent'])
        labels = [tuple(np.random.default_rng(i).permutation(4)) for i in range(4)]
        changed = {(i,j):compose(labels[j],compose(p,inverse(labels[i])))
                   for (i,j),p in edges.items()}
        self.assertTrue(synchronize(4,changed)['consistent'])
        for i,j in edges:
            self.assertEqual(compose(frames[j],inverse(frames[i])),edges[i,j])

    def test_04_triangle_checks_can_miss_a_cycle(self):
        row = sparse_cycle_counterexample()
        self.assertFalse(row['global_consistency'])
        self.assertEqual(row['triangle_checks_available'],0)
        self.assertEqual(row['whole_cycle_permutation'],[1,2,0])

    def test_05_consistency_is_not_physical_truth(self):
        row = all_consistent_but_wrong()
        self.assertTrue(row['all_cycle_tests_pass'] and row['block_psd'])
        self.assertEqual(row['wrong_physical_matches_per_pair'],[6,6,6])
        self.assertEqual(row['block_rank'],6)

    def test_06_correct_block_spectrum(self):
        frames = [(0,1,2),(1,2,0),(2,0,1),(1,0,2)]
        edges = {(i,j):compose(frames[j],inverse(frames[i]))
                 for i in range(4) for j in range(4) if i!=j}
        vals = np.linalg.eigvalsh(block_matrix(4,edges))
        self.assertLess(np.max(abs(vals-np.array([0]*9+[4]*3))),1e-12)

    def test_07_noise_budget_all_phases(self):
        row = confidence_case()
        self.assertEqual(row['all_phases_queries_before_first_error_upper'],1230)
        self.assertLessEqual(row['failure_upper_exponential'],.01)
        exact = triple_run()
        self.assertLessEqual(exact['marker_read_blocks'],1230)
        for seed in range(30):
            sample = triple_run(.1,row['repetitions'],seed)
            self.assertTrue(sample['correct'])
            self.assertEqual(sample['raw_marker_reads'],
                             sample['marker_read_blocks']*row['repetitions'])

    def test_08_partial_domains_and_disconnection(self):
        self.assertEqual(partial_triple_failures({0:1,2:3},{1:4},{0:4,2:8}),(1,0))
        self.assertEqual(partial_triple_failures({0:1},{1:4},{0:5}),(1,1))
        edges = {}
        edge_pair(edges,0,1,(0,1))
        with self.assertRaises(ValueError):
            synchronize(3,edges)


if __name__ == '__main__':
    main(__name__,'chart_consistency_audit',report)
