"""Round 333: wavelength-dependent response of a stipulated massive scalar envelope."""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main
from conformal_probe_motion_audit import integrate

RHO0=.03
H0=np.sqrt(RHO0/3)
END=20.
DELTA0=1e-4


def background(t):
    z=1+1.5*H0*t
    return z**(2/3),H0/z,RHO0/z**2


def jeans(mass,t=0.):
    a,_,rho=background(t)
    return (2*rho*mass**2*a**4)**.25


def coefficients(t,k,mass):
    a,h,rho=background(t)
    q=k*k/(2*mass*a*a)
    r=mass*a*a*rho/(k*k)
    return h,q,r


def density_rhs(y,k,mass):
    t,delta,velocity=y
    h,q,r=coefficients(t,k,mass)
    return np.array([1.,velocity,-2*h*velocity-(q*q-q*r)*delta])


def envelope_rhs(y,k,mass):
    t,u,v=y
    _,q,r=coefficients(t,k,mass)
    return np.array([1.,q*v,-(q-r)*u])


@lru_cache(maxsize=24)
def solve(k,mass=64.,kind='density',steps=1024):
    if kind=='density':
        initial=np.array([0.,DELTA0,H0*DELTA0])
        return integrate(lambda y:density_rhs(y,k,mass),initial,END,steps)
    _,q,_=coefficients(0.,k,mass)
    initial=np.array([0.,DELTA0/2,H0*DELTA0/(2*q)])
    return integrate(lambda y:envelope_rhs(y,k,mass),initial,END,steps)


def exact_free_frequency(k,mass):
    # Rationalized expression avoids cancellation of two large rest energies.
    return k*k/(np.sqrt(mass*mass+k*k)+mass)


def summarize(k,mass=64.):
    path=solve(k,mass)
    times=path[:,0]
    a,_,_=background(times)
    dust=DELTA0*a
    h,q,r=coefficients(0.,k,mass)
    exact=exact_free_frequency(k,mass)
    return {'mass':mass,'comoving_k':float(k),'initial_k_over_kJ':float(k/jeans(mass)),
            'initial_k_over_M':float(k/mass),'initial_H_over_M':float(h/mass),
            'initial_H_over_physical_k':float(h/k),
            'initial_density_frequency_squared':float(q*q-q*r),
            'final_delta_over_initial':float(path[-1,1]/DELTA0),
            'dust_final_delta_over_initial':float(a[-1]),
            'maximum_relative_to_initial_dust_difference':float(np.max(abs(path[:,1]-dust))/DELTA0),
            'free_envelope_frequency_relative_error':float(abs(k*k/(2*mass)-exact)/exact),
            'linear_envelope_density_equation_agreement':float(np.max(abs(2*solve(k,mass,'envelope')[:,1]-path[:,1])))}


def report():
    return {'round':333,
            'scope':'Classical nonrelativistic scalar-envelope plus Newtonian Poisson equations, linearized around the stipulated pressureless FLRW background of round 332; kappa=1, c=1. The expanding perturbation equation and its complex-envelope form are compared, with an independent exact free Klein-Gordon dispersion check. This is not a full inhomogeneous Einstein-Klein-Gordon matching, a derivation of quantum mechanics, or an identification of dark matter.',
            'mass_64_wavelength_sequence':[summarize(ratio*jeans(64.)) for ratio in (.25,.5,1.,2.)],
            'same_background_different_mass':[summarize(3.,mass) for mass in (16.,64.)],
            'next_interface':'Background agreement does not guarantee agreement of matter response; separately audit common operational geometry and state the assumptions needed for full equivalence.'}


class Checks(unittest.TestCase):
    def test_01_background_satisfies_friedmann_and_dust(self):
        times=np.linspace(0.,END,100); a,h,rho=background(times)
        np.testing.assert_allclose(3*h*h,rho,atol=2e-17)
        np.testing.assert_allclose(a**3*rho,RHO0,atol=3e-17)

    def test_02_poisson_feedback_coefficient(self):
        for t in (0.,3.,20.):
            for k in (1.,3.,7.):
                _,q,r=coefficients(t,k,64.)
                self.assertAlmostEqual(q*r,background(t)[2]/2,places=15)

    def test_03_jeans_zero_and_redshift(self):
        for t in (0.,3.,20.):
            k=jeans(64.,t); _,q,r=coefficients(t,k,64.)
            self.assertLess(abs(q*q-q*r),1e-16)
            self.assertAlmostEqual(k/jeans(64.),background(t)[0]**.25,places=14)

    def test_04_free_kg_dispersion_has_bounded_envelope_error(self):
        for x in (.01,.05,.125,.2):
            mass=64.; k=mass*x; exact=exact_free_frequency(k,mass)
            relative=(k*k/(2*mass)-exact)/exact
            self.assertGreater(relative,0.)
            self.assertLessEqual(relative,x*x/4+1e-15)

    def test_05_declared_scale_separation_window(self):
        for mass,k in [(64.,r*jeans(64.)) for r in (.25,.5,1.,2.)]+[(16.,3.),(64.,3.)]:
            self.assertLess(k/mass,.2)
            self.assertLess(H0/mass,.01)
            self.assertLess(H0/k,.11)

    def test_06_density_and_complex_envelope_forms_agree(self):
        for ratio in (.25,.5,1.,2.):
            r=summarize(ratio*jeans(64.))
            self.assertLess(r['linear_envelope_density_equation_agreement'],2e-11)

    def test_07_dust_growing_mode_is_a_solution_when_gradient_is_removed(self):
        for t in (0.,3.,20.):
            a,h,rho=background(t)
            delta=DELTA0*a; velocity=h*delta
            acceleration=-.5*h*h*delta
            self.assertLess(abs(acceleration+2*h*velocity-rho*delta/2),1e-18)

    def test_08_gradient_response_is_not_pressureless_at_all_k(self):
        low=summarize(.25*jeans(64.)); high=summarize(2*jeans(64.))
        self.assertLess(low['maximum_relative_to_initial_dust_difference'],.02)
        self.assertGreater(high['maximum_relative_to_initial_dust_difference'],1.)
        self.assertLess(low['initial_density_frequency_squared'],0.)
        self.assertGreater(high['initial_density_frequency_squared'],0.)

    def test_09_same_background_does_not_determine_mass_or_response(self):
        light=summarize(3.,16.); heavy=summarize(3.,64.)
        self.assertEqual(light['dust_final_delta_over_initial'],heavy['dust_final_delta_over_initial'])
        self.assertGreater(light['initial_density_frequency_squared'],0.)
        self.assertLess(heavy['initial_density_frequency_squared'],0.)
        self.assertGreater(abs(light['final_delta_over_initial']-heavy['final_delta_over_initial']),.5)

    def test_10_jeans_scale_keeps_an_independent_mass_parameter(self):
        self.assertAlmostEqual(jeans(64.)/jeans(16.),2.,places=14)

    def test_11_linearization_is_consistent_for_all_tested_modes(self):
        for k in [r*jeans(64.) for r in (.25,.5,1.,2.)]:
            mode=solve(k,64.,'envelope')
            self.assertLess(np.max(np.linalg.norm(mode[:,1:],axis=1)),.01)

    def test_12_time_step_error_is_below_wavelength_effect(self):
        k=2*jeans(64.)
        fine=solve(k); coarse=solve(k,steps=512)
        self.assertLess(np.max(abs(fine[::2,1]-coarse[:,1]))/DELTA0,1e-6)


if __name__=='__main__':
    main(__name__,'scalar_spatial_response_audit',report)
