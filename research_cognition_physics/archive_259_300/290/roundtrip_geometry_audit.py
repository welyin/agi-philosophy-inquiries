"""Round 290: round-trip metrics, directed loops and synchronization choices.

Finite graph specialization of standard synchronization/round-trip methods
(compare Minguzzi, arXiv:1009.3005). No Lorentz metric or gravity is assumed
to follow from these classical delay-network calculations.
"""
from collections import deque
import unittest
import numpy as np
from growing_stream_audit import main
from affine_clock_audit import AffineNetwork, timestamp_data, calibrate, fixture
from observable_cover_audit import cube


def split_intercepts(h):
    symmetric = {(i,j):(h[i,j]+h[j,i])/2 for i,j in h}
    antisymmetric = {(i,j):(h[i,j]-h[j,i])/2 for i,j in h}
    if any(s<=0 for s in symmetric.values()):
        raise ValueError('positive corrected pair round trips required')
    return symmetric,antisymmetric


def symmetric_zero_points(count,h):
    symmetric,anti = split_intercepts(h)
    potentials = {0:0.}
    queue = deque([0])
    while queue:
        i = queue.popleft()
        for u,j in h:
            if u==i and j not in potentials:
                potentials[j] = potentials[i]+anti[i,j]
                queue.append(j)
    if len(potentials)!=count:
        raise ValueError('connected graph required')
    residual = {(i,j):anti[i,j]-(potentials[j]-potentials[i])
                for i,j in h if i<j}
    return potentials,residual


def shortest_paths(count,weights):
    matrix = np.full((count,count),np.inf)
    np.fill_diagonal(matrix,0.)
    for (i,j),weight in weights.items():
        matrix[i,j] = min(matrix[i,j],weight)
    for k in range(count):
        matrix = np.minimum(matrix,matrix[:,k,None]+matrix[None,k,:])
    if np.min(np.diag(matrix)) < -1e-9:
        raise ValueError('negative cycle: inconsistent with positive stationary delays')
    if not np.all(np.isfinite(matrix)):
        raise ValueError('strongly connected graph required')
    return matrix


def metrics(count,h):
    symmetric,_ = split_intercepts(h)
    reverse_path = shortest_paths(count,symmetric)
    directed_reading = shortest_paths(count,h)
    independent_routes = (directed_reading+directed_reading.T)/2
    return reverse_path,independent_routes


def loop_sum(weights,path):
    assert path[0]==path[-1]
    return sum(weights[i,j] for i,j in zip(path,path[1:]))


def loop_duration(net,path,hold_ticks=0.):
    """Single start/end clock plus measured relay waits in calibrated units."""
    assert path[0]==path[-1]
    relative,_,_ = calibrate(timestamp_data(net),len(net.rates))
    root = path[0]
    start = stamp = 10.
    hold_correction = 0.
    for k,(i,j) in enumerate(zip(path,path[1:])):
        stamp = net.receive(i,j,stamp)
        if k<len(path)-2:
            stamp += hold_ticks
            hold_correction += hold_ticks/relative[j]
    return float((stamp-start)/relative[root]-hold_correction)


def triangle(clockwise=1.,counterclockwise=3.):
    delays = {}
    for i,j in ((0,1),(1,2),(2,0)):
        delays[i,j],delays[j,i] = clockwise,counterclockwise
    return AffineNetwork([1.,.8,1.3],[4.,-3.,9.],delays)


def summarize(name,net):
    n = len(net.rates)
    _,h,_ = calibrate(timestamp_data(net),n)
    symmetric,anti = split_intercepts(h)
    potential,residual = symmetric_zero_points(n,h)
    reverse_path,independent = metrics(n,h)
    return {'name':name,'nodes':n,'undirected_edges':len(h)//2,
            'independent_cycles':len(h)//2-n+1,
            'maximum_symmetric_gauge_residual':float(max(abs(x) for x in residual.values())),
            'all_edges_can_be_symmetric_in_one_zero_point_choice':bool(max(abs(x) for x in residual.values())<1e-10),
            'same_reversed_path_metric':reverse_path.tolist(),
            'independently_optimized_roundtrip_metric':independent.tolist()}


def report():
    biased,reciprocal,weak = triangle(),triangle(2.,2.),triangle(1.8,2.2)
    _,h,_ = calibrate(timestamp_data(biased),3)
    symmetric,anti = split_intercepts(h)
    forward,backward = [0,1,2,0],[0,2,1,0]
    cases = [summarize('biased_triangle',biased),summarize('reciprocal_triangle',reciprocal),
             summarize('weakly_biased_triangle',weak),summarize('five_edge_network',fixture())]
    return {'round':290,
            'scope':'Two operational time metrics and loop-delay invariants on a supplied stationary classical network; no derived meter scale, light speed, Lorentz symmetry or gravity.',
            'triangle_comparison':{
                'both_models_pair_roundtrips':[4.,4.,4.],
                'biased_clockwise_loop':loop_duration(biased,forward),
                'biased_counterclockwise_loop':loop_duration(biased,backward),
                'reciprocal_clockwise_loop':loop_duration(reciprocal,forward),
                'reciprocal_counterclockwise_loop':loop_duration(reciprocal,backward),
                'biased_antisymmetric_loop_sum':float(loop_sum(anti,forward)),
                'processing_corrected_biased_loop':loop_duration(biased,forward,7.),
                'extra_closed_loop_probe_transmissions':2*(len(forward)-1),
                'interpretation':'Identical immediate edge echoes do not determine multi-edge directed circulation.'},
            'cases':cases,
            'negative_local_intercept_example':float(min(h.values())),
            'resource_scope':'Local endpoint timestamps and stable additive relay rules are assumed. Loop probes cost their actual hops and wait corrections; message encoding, clock energy and bandwidth remain unmodeled.',
            'logical_limits':['Zero circulation suffices for the two metrics to agree; equality of the metrics alone does not imply zero circulation.',
                              'A nonzero loop-delay asymmetry is compatible with a classical asymmetric network; it is not a proof of rotation, curvature or Sagnac physics.',
                              'Choosing the reference clock rescales both metrics; multiplying by a speed to obtain lengths requires extra physical input.']}


class Audit(unittest.TestCase):
    def test_01_parameter_gauge_and_reference_units(self):
        net = fixture()
        _,h,_ = calibrate(timestamp_data(net),4)
        _,other,_ = calibrate(timestamp_data(net.shifted([0,.2,-.3,.4])),4)
        left,right = metrics(4,h),metrics(4,other)
        for a,b in zip(left,right):
            self.assertLess(np.max(abs(a-b)),1e-12)
        _,h2,_ = calibrate(timestamp_data(net),4,reference=2)
        for a,b in zip(left,metrics(4,h2)):
            self.assertLess(np.max(abs(b-a*(net.rates[2]/net.rates[0]))),1e-12)

    def test_02_tree_can_be_symmetrized(self):
        net = AffineNetwork([1.,.8,1.2,1.4],[2.,-4.,1.,7.],
                            {(0,1):1.,(1,0):3.,(1,2):2.,(2,1):1.,(1,3):4.,(3,1):2.})
        relative,h,_ = calibrate(timestamp_data(net),4)
        symmetric,_ = split_intercepts(h)
        potential,residual = symmetric_zero_points(4,h)
        self.assertLess(max(abs(x) for x in residual.values()),1e-12)
        model = AffineNetwork([relative[i] for i in range(4)],
                              [relative[i]*potential[i] for i in range(4)],symmetric)
        _,h2,_ = calibrate(timestamp_data(model),4)
        self.assertLess(max(abs(h[e]-h2[e]) for e in h),1e-12)
        a,b = metrics(4,h)
        self.assertLess(np.max(abs(a-b)),1e-12)

    def test_03_same_echo_different_loops(self):
        a,b = triangle(),triangle(2.,2.)
        for i,j in a.delays:
            self.assertEqual(a.delays[i,j]+a.delays[j,i],b.delays[i,j]+b.delays[j,i])
        self.assertAlmostEqual(loop_duration(a,[0,1,2,0]),3.)
        self.assertAlmostEqual(loop_duration(a,[0,2,1,0]),9.)
        self.assertAlmostEqual(loop_duration(b,[0,1,2,0]),6.)

    def test_04_metric_equality_not_zero_circulation(self):
        _,h,_ = calibrate(timestamp_data(triangle(1.8,2.2)),3)
        a,b = metrics(3,h)
        self.assertLess(np.max(abs(a-b)),1e-12)
        _,anti = split_intercepts(h)
        self.assertAlmostEqual(loop_sum(anti,[0,1,2,0]),-.6)

    def test_05_two_route_policies_differ(self):
        _,h,_ = calibrate(timestamp_data(triangle()),3)
        a,b = metrics(3,h)
        self.assertAlmostEqual(a[0,1],2.)
        self.assertAlmostEqual(b[0,1],1.5)
        self.assertTrue(np.all(b<=a+1e-12))

    def test_06_metric_axioms_and_ellipticity_bound(self):
        rows = cube()
        for seed in range(6):
            rng = np.random.default_rng(seed)
            delays = {(i,j):float(rng.uniform(.4,3.)) for i,row in enumerate(rows) for j in row}
            net = AffineNetwork(rng.uniform(.7,1.4,8),rng.uniform(-20.,20.,8),delays)
            _,h,_ = calibrate(timestamp_data(net),8)
            a,b = metrics(8,h)
            epsilon = min(2*delays[i,j]/(delays[i,j]+delays[j,i]) for i,j in delays)
            self.assertTrue(np.all(b>=epsilon*a-1e-11))
            self.assertTrue(np.all(b<=a+1e-11))
            for d in (a,b):
                self.assertLess(np.max(abs(d-d.T)),1e-11)
                self.assertGreater(np.min(d[~np.eye(8,dtype=bool)]),0)
                for k in range(8):
                    self.assertTrue(np.all(d<=d[:,k,None]+d[None,k,:]+1e-11))

    def test_07_loops_and_relay_waits(self):
        net = triangle()
        _,h,_ = calibrate(timestamp_data(net),3)
        for path in ([0,1,2,0],[0,2,1,0],[1,2,0,1]):
            duration = loop_duration(net,path)
            self.assertAlmostEqual(duration,loop_sum(h,path))
            self.assertAlmostEqual(duration,loop_duration(net,path,11.))

    def test_08_clock_zero_changes_and_negative_cycles(self):
        _,h,_ = calibrate(timestamp_data(triangle()),3)
        self.assertLess(min(h.values()),0.)
        f = [3.,-2.,1.]
        changed = {(i,j):value+f[j]-f[i] for (i,j),value in h.items()}
        for a,b in zip(metrics(3,h),metrics(3,changed)):
            self.assertLess(np.max(abs(a-b)),1e-12)
        self.assertAlmostEqual(loop_sum(h,[0,1,2,0]),loop_sum(changed,[0,1,2,0]))
        with self.assertRaises(ValueError):
            shortest_paths(3,{(0,1):-2.,(1,2):-2.,(2,0):1.})


if __name__ == '__main__':
    main(__name__,'roundtrip_geometry_audit',report)

