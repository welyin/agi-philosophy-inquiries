"""Round 329: exact convex heavy-field minimization versus a local quartic truncation."""
import unittest
import numpy as np
from growing_stream_audit import main
from closed_species_exchange_audit import ALPHA, MASSES


def moments(chi):
    weight=MASSES**2*np.asarray(chi)**2
    return .5*np.sum(weight),float(ALPHA@weight),float((ALPHA**2)@weight)


def full_potential(phi,chi,mass=4.):
    return .5*mass**2*phi**2+.5*np.sum(MASSES**2*np.asarray(chi)**2*np.exp(2*ALPHA*phi))


def stationary(chi,mass=4.):
    if mass<=0:
        raise ValueError('Strict convexity bound uses mass>0.')
    _,j,_=moments(chi)
    lower,upper=-j/mass**2,0.
    for _ in range(80):
        mid=(lower+upper)/2
        slope=mass**2*mid+np.sum(ALPHA*MASSES**2*np.asarray(chi)**2*np.exp(2*ALPHA*mid))
        if slope>0:
            upper=mid
        else:
            lower=mid
    return (lower+upper)/2


def approximations(chi,mass=4.):
    v0,j,b=moments(chi)
    return {'quartic':v0-j*j/(2*mass**2),
            'sextic':v0-j*j/(2*mass**2)+b*j*j/mass**4,
            'upper_bound':v0-j*j/(2*(mass**2+2*b)),
            'quartic_error_bound':b*j*j/(mass**2*(mass**2+2*b)),
            'amplitude_parameter':2*b/mass**2}


def example(amplitude,mass=4.):
    chi=amplitude*np.array([.8,.6])
    phi=stationary(chi,mass)
    value=full_potential(phi,chi,mass)
    approx=approximations(chi,mass)
    return {'amplitude':amplitude,'phi_minimum':float(phi),'exact_potential':float(value),
            **{k:float(v) for k,v in approx.items()},
            'quartic_error':float(value-approx['quartic']),
            'sextic_error':float(approx['sextic']-value)}


def report():
    return {'round':329,
            'scope':'Pointwise classical potential minimization in the positive-alpha exponential-mass model of round 325, with heavy mass set to 4. Strict convexity and analytic bounds hold in this specified positive-coupling branch. Spatial/time derivatives, gravity backreaction and quantum corrections are not included in this potential-only elimination.',
            'small_amplitude_sequence':[example(a) for a in (.2,.4,.8)],
            'outside_small_amplitude_window':example(8.),
            'next_interface':'Include the induced derivative term as well as the corrected potential; compare reduced light-field dynamics against the full heavy-light model with explicitly prepared heavy initial data.'}


class Checks(unittest.TestCase):
    def test_01_stationary_equation_and_positive_curvature(self):
        for a in (.2,1.,8.):
            chi=a*np.array([.8,.6]); phi=stationary(chi)
            w=MASSES**2*chi**2*np.exp(2*ALPHA*phi)
            self.assertLess(abs(16*phi+ALPHA@w),2e-13)
            self.assertGreaterEqual(16+2*(ALPHA**2)@w,16.)

    def test_02_exact_minimum_within_proved_field_bracket(self):
        for a in (.1,1.,8.):
            chi=a*np.array([.8,.6]); phi=stationary(chi)
            self.assertLessEqual(-moments(chi)[1]/16,phi)
            self.assertLessEqual(phi,0.)

    def test_03_analytic_two_sided_potential_bounds(self):
        rng=np.random.default_rng(329)
        for _ in range(40):
            chi=rng.uniform(-5.,5.,2)
            value=full_potential(stationary(chi),chi)
            r=approximations(chi)
            self.assertGreaterEqual(value,r['quartic']-1e-13)
            self.assertLessEqual(value,r['upper_bound']+1e-13)
            self.assertLessEqual(value-r['quartic'],r['quartic_error_bound']+1e-13)

    def test_04_quartic_error_is_sixth_order_in_amplitude(self):
        errors=[example(a)['quartic_error'] for a in (.2,.4,.8)]
        for ratio in (errors[1]/errors[0],errors[2]/errors[1]):
            self.assertGreater(ratio,55.)
            self.assertLess(ratio,65.)

    def test_05_sextic_improvement_is_eighth_order(self):
        errors=[example(a)['sextic_error'] for a in (.2,.4,.8)]
        for ratio in (errors[1]/errors[0],errors[2]/errors[1]):
            self.assertGreater(ratio,220.)
            self.assertLess(ratio,260.)

    def test_06_negative_truncated_potential_is_outside_its_window(self):
        r=example(8.)
        self.assertLess(r['quartic'],0.)
        self.assertGreater(r['exact_potential'],0.)
        self.assertGreater(r['amplitude_parameter'],1.)

    def test_07_minimized_force_obeys_envelope_derivative(self):
        chi=np.array([.7,-.5]); phi=stationary(chi); h=1e-5
        exact=MASSES**2*np.exp(2*ALPHA*phi)*chi
        numerical=[]
        for d in np.eye(2):
            plus,minus=chi+h*d,chi-h*d
            numerical.append((full_potential(stationary(plus),plus)-full_potential(stationary(minus),minus))/(2*h))
        np.testing.assert_allclose(numerical,exact,atol=2e-11)

    def test_08_retained_sextic_term_improves_inside_window(self):
        for a in (.2,.4,.8):
            r=example(a)
            self.assertLess(r['sextic_error'],r['quartic_error']/10)
            self.assertLess(r['amplitude_parameter'],.04)

    def test_09_zero_matter_has_no_tree_level_potential_shift(self):
        chi=np.zeros(2)
        self.assertEqual(stationary(chi),0.)
        self.assertEqual(full_potential(0.,chi),0.)


if __name__=='__main__':
    main(__name__,'heavy_scalar_potential_audit',report)
