"""Round 335: local cone convergence, accumulated ray lag, and late protocols."""
import unittest
import numpy as np
from growing_stream_audit import main
from disformal_probe_geometry_audit import distance as old_distance, arrival as old_arrival

B=.4
NODES,WEIGHTS=np.polynomial.legendre.leggauss(80)


def deficit(t,b=B):
    q=b/np.asarray(t)**2
    if b<0 or np.any(q>=1):
        raise ValueError('Require nonnegative coupling and positive probe lapse.')
    return q/(1+np.sqrt(1-q))


def lag(start,end=np.inf,b=B,nodes=NODES,weights=WEIGHTS):
    deficit(start,b)
    if end<start:
        raise ValueError('End precedes emission.')
    lower=0. if np.isinf(end) else (start/end)**(2/3)
    u=(1+lower)/2+(1-lower)*nodes/2
    integrand=u/(1+np.sqrt(1-b/start**2*u**3))
    return 1.5*b*start**(-4/3)*(1-lower)/2*float(weights@integrand)


def lag_bounds(start,end=np.inf,b=B):
    deficit(start,b)
    power=start**(-4/3)-(0. if np.isinf(end) else end**(-4/3))
    return 3*b*power/8,3*b*power/(4*(1+np.sqrt(1-b/start**2)))


def arrival_delay(start,arrival_a,b=B):
    if arrival_a<start:
        raise ValueError('Arrival precedes emission.')
    if b==0 or arrival_a==start:
        return 0.
    def residual(delay):
        # Stable conformal distance between arrival_a and arrival_a+delay.
        advance=1.5*arrival_a**(2/3)*np.expm1((2/3)*np.log1p(delay/arrival_a))
        return advance-lag(start,arrival_a+delay,b)
    left=0.; right=2*arrival_a**(1/3)*lag(start,b=b)
    while residual(right)<0:
        right*=2
    for _ in range(70):
        mid=(left+right)/2
        if residual(mid)<0:
            left=mid
        else:
            right=mid
    return (left+right)/2


def late_protocol(start,physical_length=.5,b=B):
    travel_a=start*np.expm1(1.5*np.log1p(2*physical_length/(3*start)))
    delay=arrival_delay(start,start+travel_a,b)
    return {'emission_time':start,'fixed_initial_physical_length':physical_length,
            'initial_speed_deficit':float(deficit(start,b)),
            'A_travel_time':float(travel_a),'B_minus_A_arrival_time':float(delay),
            'delay_over_leading_late_estimate':float(delay/(physical_length*b/(2*start**2)))}


def report():
    memory=lag(1.)
    return {'round':335,
            'scope':'The constant-coupling classical weak-probe branch of round 334, a=t^(1/3), D_A=0, D_B=0.6 and b=2D_B/3=0.4, analytically extended over t>=1. Exact characteristic rays in the stipulated second-order model; ideal A clocks; no finite probe backreaction or wave-packet/UV completion claim. Late background dilution and lower probe frequency are distinct limits.',
            'integrated_comoving_lag_to_infinity':memory,
            'analytic_infinite_lag_bounds':list(lag_bounds(1.)),
            'fixed_early_emission_increasing_receiver_distance':[
                {'A_arrival_time':t,'instantaneous_speed_deficit_at_A_arrival':float(deficit(t)),
                 'B_minus_A_arrival_time':float(arrival_delay(1.,t)),
                 'delay_divided_by_a_at_A_arrival':float(arrival_delay(1.,t)/t**(1/3)),
                 'delay_over_A_arrival_time':float(arrival_delay(1.,t)/t)}
                for t in (2.,10.,100.,1000.)],
            'late_emission_fixed_physical_length':[late_protocol(t) for t in (1.,2.,4.,16.,64.)],
            'next_interface':'Determine whether the background dilution survives field-dependent couplings and finite matter occupation; do not equate vanishing instantaneous mismatch with erasure of measurement history.'}


class Checks(unittest.TestCase):
    def test_01_rationalized_deficit_matches_speed(self):
        t=np.array([1.,2.,8.,64.])
        np.testing.assert_allclose(deficit(t),1-np.sqrt(1-B/t**2),atol=1e-16)

    def test_02_pointwise_bounds(self):
        for start in (1.,2.,8.):
            t=np.geomspace(start,start*1000,100)
            value=deficit(t)
            self.assertTrue(np.all(value>=B/(2*t*t)*(1-1e-15)))
            self.assertTrue(np.all(value<=B/(t*t*(1+np.sqrt(1-B/start**2)))*(1+1e-15)))

    def test_03_lag_agrees_with_independent_legacy_simpson_integral(self):
        for end in (1.5,2.,4.):
            exact_a=1.5*(end**(2/3)-1)
            self.assertLess(abs(lag(1.,end)-(exact_a-old_distance(end,.6))),2e-14)

    def test_04_finite_and_infinite_lag_bounds(self):
        for start in (1.,2.,8.):
            for end in (2*start,20*start,np.inf):
                low,high=lag_bounds(start,end)
                self.assertLessEqual(low,lag(start,end)+1e-15)
                self.assertLessEqual(lag(start,end),high+1e-15)

    def test_05_infinite_tail_is_controlled(self):
        for end in (4.,16.,64.):
            remainder=lag(1.)-lag(1.,end)
            self.assertAlmostEqual(remainder,lag(end),places=14)
            low,high=lag_bounds(end)
            self.assertLessEqual(low,remainder+1e-15)
            self.assertLessEqual(remainder,high+1e-15)

    def test_06_arrival_matches_the_original_protocol(self):
        ta=old_arrival(.5,0.); tb=old_arrival(.5,.6)
        self.assertAlmostEqual(arrival_delay(1.,ta),tb-ta,places=12)

    def test_07_instantaneous_convergence_does_not_erase_early_ray_history(self):
        early=arrival_delay(1.,10.); far=arrival_delay(1.,1000.)
        self.assertLess(deficit(1000.),deficit(10.)/9000)
        self.assertGreater(far,early)
        self.assertLess(abs(far/1000.**(1/3)-lag(1.))/lag(1.),.002)
        self.assertLess(far/1000.,early/10.)

    def test_08_late_local_protocol_converges(self):
        rows=[late_protocol(t) for t in (4.,16.,64.)]
        self.assertGreater(rows[0]['B_minus_A_arrival_time']/rows[1]['B_minus_A_arrival_time'],14.)
        self.assertGreater(rows[1]['B_minus_A_arrival_time']/rows[2]['B_minus_A_arrival_time'],15.)
        self.assertLess(abs(rows[-1]['delay_over_leading_late_estimate']-1),.02)

    def test_09_probe_frequency_does_not_select_a_common_principal_cone(self):
        # Local principal symbol: omega^2=(1-b/t^2) k^2/a^2.
        t=2.; a=t**(1/3)
        speeds=[a*np.sqrt((1-B/t**2)*k*k/a**2)/k for k in (1.,10.,100.)]
        np.testing.assert_allclose(speeds,np.sqrt(1-B/t**2),atol=2e-16)
        self.assertGreater(1-speeds[0],.05)

    def test_10_quadrature_is_subdominant(self):
        nodes,weights=np.polynomial.legendre.leggauss(40)
        self.assertLess(abs(lag(1.)-lag(1.,nodes=nodes,weights=weights)),2e-15)

    def test_11_domain_and_identical_coupling_limit(self):
        self.assertEqual(lag(1.,4.,b=0.),0.)
        self.assertEqual(arrival_delay(1.,4.,b=0.),0.)
        with self.assertRaises(ValueError):
            lag(.5,b=B)


if __name__=='__main__':
    main(__name__,'late_probe_memory_audit',report)
