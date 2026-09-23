"""Round 334: identical Einstein background with inequivalent matter characteristic cones."""
import unittest
import numpy as np
from growing_stream_audit import main


def background(t):
    return t**(1/3),1/(3*t),np.sqrt(2/3)/t


def metric(t,coupling=0.):
    a,_,velocity=background(t)
    c=1-coupling*velocity**2
    if c<=0:
        raise ValueError('Probe metric is not on the positive lapse branch.')
    return np.diag([-c,a*a,a*a,a*a])


def cone_speed(t,coupling=0.):
    _,_,velocity=background(t)
    return np.sqrt(1-coupling*velocity**2)


def kinetic_density(t,coupling=0.):
    g=metric(t,coupling)
    return np.sqrt(-np.linalg.det(g))*np.linalg.inv(g)


def probe_lagrangian(lapse,scale,phi_velocity,probe_velocity,gradient_squared,coupling):
    c=lapse*lapse-coupling*phi_velocity**2
    return .5*scale**3*probe_velocity**2/np.sqrt(c)-.5*scale*np.sqrt(c)*gradient_squared


def distance(end,coupling=0.,samples=4096):
    # Composite Simpson integration of radial null propagation, from t=1.
    t=np.linspace(1.,end,samples+1); a,_,_=background(t)
    integrand=cone_speed(t,coupling)/a
    weights=np.ones(samples+1); weights[1:-1:2]=4; weights[2:-1:2]=2
    return (end-1)/(3*samples)*float(weights@integrand)


def arrival(comoving_distance,coupling=0.):
    left,right=1.,4.
    assert distance(right,coupling)>comoving_distance
    for _ in range(60):
        mid=(left+right)/2
        if distance(mid,coupling)<comoving_distance:
            left=mid
        else:
            right=mid
    return (left+right)/2


def report():
    arrival_a=arrival(.5,0.); arrival_b=arrival(.5,.6)
    rows=[]
    for t in (1.,2.,4.):
        a,h,v=background(t)
        rows.append({'time':t,'a':float(a),'H':float(h),'scalar_density':float(v*v/2),
                     'einstein_constraint':float(3*h*h-v*v/2),
                     'speed_A_in_g_orthonormal_frame':1.,
                     'speed_B_in_g_orthonormal_frame':float(cone_speed(t,.6)),
                     'generalized_covariant_metric_eigenvalues_B_over_A':[float(x) for x in np.linalg.eigvals(np.linalg.solve(metric(t),metric(t,.6)))]})
    return {'round':334,
            'scope':'An explicitly stipulated generally covariant classical Einstein-canonical-scalar action, with two additional scalar probes minimally coupled to g_i=g+D_i*dphi*dphi. The exactly zero probe background is an Einstein solution for every D_i; nonzero probe propagation is tested to linear order, with backreaction quadratic in probe amplitude. D_A=0, D_B=0.6, t in [1,4], positive probe lapses. This is a counterexample to inferring universal matter geometry from Einstein background matching or covariance alone, not a counterexample to full GR with universal matter coupling already required.',
            'common_background_different_probe_cones':rows,
            'arrival_at_comoving_distance_half':{'A':arrival_a,'B':arrival_b,
                'difference_in_A_comoving_proper_time':arrival_b-arrival_a},
            'coordinate_and_common_field_redefinitions':'The generalized eigenvalues of the two covariant metrics are invariant under a common coordinate transformation and common positive conformal scaling; unequal cones remain unequal.',
            'next_interface':'Make universal operational geometry an explicit selection target. Search for independent microscopic constraints or bounds on species-dependent derivative couplings instead of inferring their absence from the Einstein equation.'}


class Checks(unittest.TestCase):
    def test_01_background_is_exact_einstein_scalar_solution(self):
        for t in np.linspace(1.,4.,31):
            a,h,v=background(t)
            self.assertLess(abs(3*h*h-v*v/2),1e-16)
            self.assertLess(abs(-1/(3*t*t)+v*v/2),1e-16)
            self.assertLess(abs(-v/t+3*h*v),5e-16)

    def test_02_positive_probe_metric_and_hamiltonian_coefficients(self):
        for t in np.linspace(1.,4.,31):
            for d in (0.,.6):
                g=metric(t,d); density=kinetic_density(t,d)
                self.assertEqual(np.count_nonzero(np.linalg.eigvalsh(g)<0),1)
                self.assertGreater(-density[0,0],0.)
                self.assertGreater(np.min(np.diag(density)[1:]),0.)

    def test_03_disformal_formula_matches_direct_metric(self):
        for t in (1.,2.,4.):
            _,_,v=background(t); gradient=np.array([v,0.,0.,0.])
            np.testing.assert_allclose(metric(t,.6),metric(t)+.6*np.outer(gradient,gradient),atol=2e-16)

    def test_04_characteristic_speed_from_independent_principal_matrix(self):
        for t in (1.,2.,4.):
            a,_,_=background(t)
            density=kinetic_density(t,.6)
            omega=np.sqrt(density[1,1]/(-density[0,0]))
            self.assertAlmostEqual(a*omega,cone_speed(t,.6),places=14)
            covector=np.array([-omega,1.,0.,0.])
            self.assertLess(abs(covector@density@covector),1e-15)

    def test_05_shared_coordinate_change_preserves_cone_mismatch(self):
        jac=np.array([[1.2,.2,0.,0.],[.1,.9,0.,0.],[0.,0.,1.1,.1],[0.,0.,0.,.8]])
        first,second=metric(1.),metric(1.,.6)
        before=np.sort(np.linalg.eigvals(np.linalg.solve(first,second)))
        after=np.sort(np.linalg.eigvals(np.linalg.solve(jac.T@first@jac,jac.T@second@jac)))
        np.testing.assert_allclose(after,before,atol=2e-15)
        self.assertGreater(np.ptp(before),.3)

    def test_06_common_conformal_redefinition_preserves_relative_geometry(self):
        first,second=metric(1.),metric(1.,.6)
        np.testing.assert_allclose(np.linalg.solve(7*first,7*second),np.linalg.solve(first,second),atol=2e-16)

    def test_07_probe_source_vanishes_only_at_zero_or_linear_order(self):
        a,_,v=background(1.); d=.6; step=1e-5
        for amplitude in (0.,1e-3,2e-3):
            w=amplitude; gradient_squared=.3*amplitude**2
            derivative=(probe_lagrangian(1+step,a,v,w,gradient_squared,d)-probe_lagrangian(1-step,a,v,w,gradient_squared,d))/(2*step)
            c=1-d*v*v
            expected=w*w/(2*c**1.5)+gradient_squared/(2*a*a*np.sqrt(c))
            self.assertLess(abs(-derivative/a**3-expected),1e-14)
        self.assertEqual(probe_lagrangian(1.,a,v,0.,0.,d),0.)

    def test_08_arrival_difference_is_measured_by_one_clock(self):
        a=arrival(.5,0.); b=arrival(.5,.6)
        self.assertAlmostEqual(a,(1+2*.5/3)**1.5,places=12)
        self.assertGreater(b-a,.05)
        self.assertAlmostEqual(distance(b,.6),.5,places=13)

    def test_09_ray_integration_error_is_below_arrival_effect(self):
        end=arrival(.5,.6)
        self.assertLess(abs(distance(end,.6,2048)-distance(end,.6,4096)),1e-11)

    def test_10_invalid_signature_is_rejected(self):
        with self.assertRaises(ValueError):
            metric(1.,2.)


if __name__=='__main__':
    main(__name__,'disformal_probe_geometry_audit',report)
