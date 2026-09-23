"""Round 292: FIFO arrival maps, real echoes and a quasi-static certificate.
Reuses the time-dependent shortest path framework (Brian C. Dean).
Known exogenous arrival functions are inputs, not inferred spacetime dynamics.
"""
import heapq
import itertools
import unittest
import numpy as np
from growing_stream_audit import main
from roundtrip_geometry_audit import shortest_paths


def path_arrival(arcs,path,time):
    for i,j in zip(path,path[1:]):
        time=arcs[i,j](time)
    return float(time)


def earliest(count,arcs,source,time):
    """Time-dependent Dijkstra: requires positive transit times and FIFO maps."""
    labels=np.full(count,np.inf)
    labels[source]=time
    queue=[(time,source)]
    while queue:
        t,i=heapq.heappop(queue)
        if t>labels[i]:
            continue
        for u,j in arcs:
            if u!=i:
                continue
            arrival=float(arcs[u,j](t))
            if arrival<=t:
                raise ValueError('strictly positive transit times required')
            if arrival<labels[j]:
                labels[j]=arrival
                heapq.heappush(queue,(arrival,j))
    if not np.all(np.isfinite(labels)):
        raise ValueError('strongly connected graph required')
    return labels


def simple_paths(count,arcs,source,target):
    if source==target:
        yield [source]
        return
    def walk(path):
        for u,j in arcs:
            if u==path[-1] and j not in path:
                nxt=path+[j]
                if j==target:
                    yield nxt
                else:
                    yield from walk(nxt)
    yield from walk([source])


def echo(count,arcs,source,target,time):
    at_target=earliest(count,arcs,source,time)[target]
    at_source=earliest(count,arcs,target,at_target)[source]
    return float((at_source-time)/2)


def echo_matrix(count,arcs,time):
    return np.array([[echo(count,arcs,i,j,time) for j in range(count)]
                     for i in range(count)])


def reference_metric(count,weights):
    directed=shortest_paths(count,weights)
    return (directed+directed.T)/2


def linear_chain(kappa=.1):
    return {(i,j):(lambda t,k=kappa:(1+k)*t+1.)
            for i,j in ((0,1),(1,0),(1,2),(2,1))}


def oscillating_network(seed=0,amplitude=.2):
    rng=np.random.default_rng(seed)
    edges=((0,1),(1,2),(2,3),(3,0),(0,2))
    weights={}
    arcs={}
    for i,j in edges:
        for u,v in ((i,j),(j,i)):
            base=float(rng.uniform(1.,3.))
            phase=float(rng.uniform(0,2*np.pi))
            weights[u,v]=base
            arcs[u,v]=lambda t,b=base,p=phase:t+b*(1+amplitude*np.sin(.1*t+p))
    # d/dt of every arrival map >= 1 - 3 * amplitude * .1 > 0.
    return weights,arcs


def triangle_excess(matrix):
    n=len(matrix)
    return max(float(matrix[i,k]-matrix[i,j]-matrix[j,k])
               for i,j,k in itertools.product(range(n),repeat=3))


def window_sufficient(count,max_weight,epsilon,start,end):
    """Conservative horizon for static optimal outward and return paths."""
    duration=2*(1+epsilon)*(count-1)*max_weight
    return duration, bool(start+duration<=end)


def report():
    arcs=linear_chain()
    one_way=np.array([earliest(3,arcs,i,0.) for i in range(3)])
    actual=echo_matrix(3,arcs,0.)
    weights,varying=oscillating_network()
    reference=reference_metric(4,weights)
    measured=echo_matrix(4,varying,5.)
    ratios=measured[np.triu_indices(4,1)]/reference[np.triu_indices(4,1)]
    duration,covered=window_sufficient(4,3.,.2,0.,30.)
    return {'round':292,
            'scope':'Conditional FIFO time-dependent propagation and approximate static time metrics; known arrival functions and calibrated rate units are extra inputs, not derived Lorentz geometry or gravity.',
            'reciprocal_fifo_chain':{
                'edge_delay':'1 + 0.1 * departure_time, for time >= 0',
                'one_way_durations_at_zero':one_way.tolist(),
                'actual_half_roundtrips_at_zero':actual.tolist(),
                'one_way_triangle_excess':triangle_excess(one_way),
                'echo_triangle_excess':triangle_excess(actual),
                'snapshot_metric':[[0.,1.,2.],[1.,0.,1.],[2.,1.,0.]]},
            'bounded_variation_example':{
                'relative_delay_bound':.2, 'departure_time':5.,
                'reference_static_metric':reference.tolist(),
                'actual_half_roundtrips':measured.tolist(),
                'min_upper_triangle_ratio':float(min(ratios)),
                'max_upper_triangle_ratio':float(max(ratios)),
                'maximum_origin_asymmetry':float(np.max(abs(measured-measured.T))),
                'approximate_triangle_factor':1.2/.8},
            'finite_window_certificate':{
                'nodes':4, 'maximum_reference_edge_delay':3.,
                'relative_bound':.2, 'conservative_roundtrip_horizon':duration,
                'window':[0.,30.], 'static_optimal_route_stays_in_window':covered},
            'counterexample_to_waiting_never_helping_without_fifo':{
                'immediate_departure_arrival':10.,
                'wait_one_then_depart_arrival':2.},
            'zero_point_covariance':'Local arrival maps transform by endpoint conjugation; closed-loop elapsed time is unaffected by constant display shifts.',
            'limits':['Snapshot shortest paths still form a mathematical metric, but need not equal actual journey durations.',
                      'The FIFO algorithm assumes arrival functions are known; endpoint samples do not establish their entire future behavior.',
                      'Relative variation bounds must hold at every actual traversal, not only at the initial departure.',
                      'A rate-normalized common time coordinate remains a modeling input; clock drift and endogenous routing load require a joint model.']}


class Audit(unittest.TestCase):
    def test_01_static_echo_reduces_to_round290_metric(self):
        weights,_=oscillating_network()
        arcs={e:(lambda t,d=d:t+d) for e,d in weights.items()}
        expected=reference_metric(4,weights)
        for time in (0.,3.,100.):
            np.testing.assert_allclose(echo_matrix(4,arcs,time),expected,atol=1e-12)

    def test_02_fifo_dijkstra_matches_exhaustive_paths(self):
        for seed in range(4):
            _,arcs=oscillating_network(seed)
            for time,source in itertools.product((0.,4.,30.),range(4)):
                result=earliest(4,arcs,source,time)
                for target in range(4):
                    exhaustive=min(path_arrival(arcs,p,time)
                                   for p in simple_paths(4,arcs,source,target))
                    self.assertAlmostEqual(result[target],exhaustive)

    def test_03_reciprocal_fifo_echo_fails_triangle(self):
        arcs=linear_chain(.1)
        actual=echo_matrix(3,arcs,0.)
        self.assertAlmostEqual(actual[0,1],1.05)
        self.assertAlmostEqual(actual[0,2],2.3205)
        self.assertAlmostEqual(triangle_excess(actual),.2205)
        np.testing.assert_allclose(actual,actual.T)
        self.assertGreater(actual[0,2],actual[0,1]+actual[1,2])

    def test_04_time_shifted_composition_triangle(self):
        for seed in range(3):
            _,arcs=oscillating_network(seed)
            for time in (0.,4.,20.):
                for x,y in itertools.product(range(4),repeat=2):
                    from_x=earliest(4,arcs,x,time)
                    via_y=earliest(4,arcs,y,from_x[y])
                    self.assertTrue(np.all(from_x<=via_y+1e-12))

    def test_05_uniform_variation_bounds_for_real_roundtrips(self):
        for seed in range(8):
            weights,arcs=oscillating_network(seed)
            ref=reference_metric(4,weights)
            for time in (0.,1.,10.,50.):
                measured=echo_matrix(4,arcs,time)
                self.assertTrue(np.all(measured>=.8*ref-1e-12))
                self.assertTrue(np.all(measured<=1.2*ref+1e-12))

    def test_06_launch_snapshot_does_not_bound_later_hops(self):
        # All edges have delay exactly 1 at t=0, but an echo samples later times.
        arcs=linear_chain(.1)
        self.assertEqual(arcs[0,1](0.),1.)
        self.assertGreater(echo(3,arcs,0,1,0.),1.)
        duration,okay=window_sufficient(4,3.,.2,0.,30.)
        self.assertAlmostEqual(duration,21.6)
        self.assertTrue(okay)
        self.assertFalse(window_sufficient(4,3.,.2,10.,30.)[1])

    def test_07_non_fifo_wait_counterexample(self):
        arrival=lambda t:t+(10. if t<1. else 1.)
        self.assertLess(arrival(1.),arrival(0.))
        self.assertEqual(arrival(1.),2.)
        self.assertEqual(arrival(0.),10.)

    def test_08_local_zero_point_conjugation(self):
        _,arcs=oscillating_network()
        beta=[4.,-3.,9.,2.]
        local={e:(lambda x,e=e:arcs[e](x-beta[e[0]])+beta[e[1]]) for e in arcs}
        for path in ([0,1,2,0],[1,0,3,2,1],[0,2,3]):
            for time in (0.,5.,10.):
                observed=path_arrival(local,path,time+beta[path[0]])
                self.assertAlmostEqual(observed,
                                       path_arrival(arcs,path,time)+beta[path[-1]])

    def test_09_approximate_triangle_and_origin_symmetry(self):
        for seed in range(8):
            weights,arcs=oscillating_network(seed)
            ref=reference_metric(4,weights)
            matrix=echo_matrix(4,arcs,5.)
            self.assertTrue(np.all(abs(matrix-matrix.T)<=.4*ref+1e-12))
            for x,y,z in itertools.product(range(4),repeat=3):
                self.assertLessEqual(matrix[x,z],1.5*(matrix[x,y]+matrix[y,z])+1e-12)


if __name__ == '__main__':
    main(__name__,'temporal_distance_audit',report)

