"""Round 293: local clock reparametrization and closed-loop obstructions.
Increasing clock changes preserve event order, not interval arithmetic.
"""
import unittest
import numpy as np
from growing_stream_audit import main
from temporal_distance_audit import path_arrival, linear_chain, echo_matrix


def clocks():
    scales=[1.2,.8,1.5]
    rates=[.3,.4,.25]
    offsets=[2.,-3.,4.]
    forward=[lambda t,a=a,s=s,b=b:a*np.sinh(s*t)+b
             for a,s,b in zip(scales,rates,offsets)]
    inverse=[lambda x,a=a,s=s,b=b:np.arcsinh((x-b)/a)/s
             for a,s,b in zip(scales,rates,offsets)]
    return forward,inverse


def rechart(arcs,forward,inverse):
    return {(i,j):(lambda x,i=i,j=j:forward[j](arcs[i,j](inverse[i](x))))
            for i,j in arcs}


def loop_f(x):
    return x+1.


def loop_g(x):
    return x+2.+.2*np.tanh(x)


def inverse_g(y):
    low,high=y-2.2,y-1.8
    for _ in range(70):
        mid=(low+high)/2
        if loop_g(mid)<y:
            low=mid
        else:
            high=mid
    return (low+high)/2


def single_loop_clock(x):
    """Fundamental-domain conjugacy: u(G(x))=u(x)+1, domain R."""
    value=float(x)
    count=0
    while value<0.:
        value=loop_g(value)
        count-=1
    while value>=2.:
        value=inverse_g(value)
        count+=1
    return count+value/2.


def star_network():
    return {(0,1):lambda x:x+.4, (1,0):lambda x:x+.6,
            (0,2):lambda x:x+.4,
            (2,0):lambda x:loop_g(x-.4)}


def chain_clock(t,kappa=.1):
    return np.log1p(kappa*np.asarray(t))/np.log1p(kappa)


def order_certificate(first,second,error):
    if error<0:
        raise ValueError('nonnegative output uncertainty required')
    if first+error<second-error:
        return -1
    if second+error<first-error:
        return 1
    return 0


def report():
    arcs=star_network()
    fg=[0,1,0,2,0]
    gf=[0,2,0,1,0]
    phi,inv=clocks()
    changed=rechart(arcs,phi,inv)
    original=[path_arrival(arcs,p,0.) for p in (fg,gf)]
    transformed=[path_arrival(changed,p,phi[0](0.)) for p in (fg,gf)]
    chain=linear_chain()
    raw=echo_matrix(3,chain,0.)
    after=chain_clock(2*raw)/2
    starts=np.array([-2.,0.,2.])
    single_raw=np.array([loop_g(x)-x for x in starts])
    single_fixed=[single_loop_clock(loop_g(x))-single_loop_clock(x) for x in starts]
    return {'round':293,
            'scope':'Order and conjugacy invariants of supplied deterministic local arrival maps; no physical clock law, absolute time, Lorentz symmetry or quantum noncommutativity derived.',
            'round292_chain_reclocking':{
                'new_coordinate':'u(t)=log(1+0.1*t)/log(1.1), t > -10',
                'original_half_roundtrips':raw.tolist(),
                'reclocked_half_roundtrips':after.tolist(),
                'each_edge_is_unit_translation_in_new_coordinate':True},
            'closed_loop_order_witness':{
                'start_reading':0., 'first_F_then_G':original[0],
                'first_G_then_F':original[1], 'difference':original[0]-original[1],
                'after_nonlinear_clock_change':transformed,
                'comparison_preserved':bool(transformed[0]>transformed[1]),
                'edge_transmissions_for_two_words':len(fg)+len(gf)-2,
                'comparison_with_per_output_error_0_01':order_certificate(*original,.01),
                'comparison_with_per_output_error_0_1':order_certificate(*original,.1)},
            'single_loop_straightening':{
                'input_readings':starts.tolist(), 'raw_elapsed':single_raw.tolist(),
                'new_elapsed':single_fixed,
                'fundamental_interval':[0.,2.],
                'choice_inside_fundamental_interval_is_not_unique':True},
            'boundaries':['Only order-preserving relabeling is treated as gauge; a calibrated physical clock can carry extra interval structure.',
                          'Changing the actual hardware, triggering rule or probe load is not merely changing its displayed coordinate.',
                          'Noncommuting loop maps rule out simultaneous constant-delay coordinates, not time coordinates in general.',
                          'The ideal two-word comparison assumes the same initial event and noninteracting reproducible channels; precision and launch error require separate budgets.']}


class Audit(unittest.TestCase):
    def test_01_path_conjugacy_telescopes(self):
        arcs=star_network()
        phi,inv=clocks()
        changed=rechart(arcs,phi,inv)
        for path in ([0,1,0],[0,2,0],[1,0,2],[0,1,0,2,0],[2,0,1,0,2]):
            for t in np.linspace(-3,3,13):
                got=path_arrival(changed,path,phi[path[0]](t))
                want=phi[path[-1]](path_arrival(arcs,path,t))
                self.assertAlmostEqual(got,want,places=11)

    def test_02_order_but_not_interval_ratio_is_preserved(self):
        phi,_=clocks()
        times=np.array([-2.,-1.,0.,1.,2.])
        for f in phi:
            self.assertTrue(np.all(np.diff(f(times))>0))
        displayed=phi[0](np.array([0.,1.,2.]))
        self.assertGreater((displayed[2]-displayed[1])/(displayed[1]-displayed[0]),1.)

    def test_03_constant_loop_can_look_time_dependent(self):
        phi,inv=clocks()
        transformed=lambda x:phi[0](inv[0](x)+1.)
        starts=phi[0](np.array([0.,2.,4.]))
        elapsed=np.array([transformed(x)-x for x in starts])
        self.assertGreater(np.ptp(elapsed),.3)
        np.testing.assert_allclose([inv[0](transformed(x))-inv[0](x) for x in starts],1.)

    def test_04_round292_chain_becomes_static_after_reclocking(self):
        for t in (-5.,0.,10.,30.):
            self.assertAlmostEqual(chain_clock(1.1*t+1)-chain_clock(t),1.)
        raw=echo_matrix(3,linear_chain(),0.)
        np.testing.assert_allclose(chain_clock(2*raw)/2,
                                   [[0.,1.,2.],[1.,0.,1.],[2.,1.,0.]])

    def test_05_loop_noncommutation_survives_recharting(self):
        arcs=star_network()
        phi,inv=clocks()
        changed=rechart(arcs,phi,inv)
        for x in (-2.,0.,2.):
            a=path_arrival(arcs,[0,1,0,2,0],x)
            b=path_arrival(arcs,[0,2,0,1,0],x)
            self.assertGreater(a,b)
            a2=path_arrival(changed,[0,1,0,2,0],phi[0](x))
            b2=path_arrival(changed,[0,2,0,1,0],phi[0](x))
            self.assertGreater(a2,b2)

    def test_06_constant_unequal_loop_delays_still_commute(self):
        phi,inv=clocks()
        f=lambda x:phi[0](inv[0](x)+3.)
        g=lambda x:phi[0](inv[0](x)+9.)
        for x in np.linspace(-1,5,11):
            self.assertAlmostEqual(f(g(x)),g(f(x)),places=9)
            self.assertGreater(g(x),f(x))

    def test_07_finite_output_resolution_only_certifies_large_gaps(self):
        a,b=loop_g(loop_f(0.)),loop_f(loop_g(0.))
        self.assertEqual(order_certificate(a,b,.01),1)
        self.assertEqual(order_certificate(a,b,.1),0)
        for ea in (-.01,.01):
            for eb in (-.01,.01):
                self.assertGreater(a+ea,b+eb)

    def test_08_single_loop_fundamental_domain_conjugacy(self):
        samples=np.linspace(-15.,15.,151)
        values=np.array([single_loop_clock(x) for x in samples])
        self.assertTrue(np.all(np.diff(values)>0))
        for x in samples:
            self.assertAlmostEqual(single_loop_clock(loop_g(x))-single_loop_clock(x),1.,places=11)


if __name__ == '__main__':
    main(__name__,'local_clock_gauge_audit',report)
