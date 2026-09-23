"""Round 336: finite homogeneous matter occupation changes a critical disformal branch."""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main
from conformal_probe_motion_audit import integrate

BETA=np.sqrt(6.)
D0=.6
Q0=.4


def phase_rhs(y):
    fraction,q,logh=y
    if not (0<=fraction<1 and 0<=q<1):
        raise ValueError('Outside the positive homogeneous phase branch.')
    x=np.sqrt(1-fraction)
    qprime=-q*fraction*(6/(1+x)-3*q/(1-fraction))/(1+fraction*q/(2*(1-fraction)*(1-q)))
    fprime=fraction*(qprime/(2*(1-q))-3*q*fraction)
    return np.array([fprime,qprime,-3+1.5*q*fraction])


@lru_cache(maxsize=16)
def solve_phase(fraction=.1,end=100.,steps=4096):
    return integrate(phase_rhs,np.array([fraction,Q0,0.]),end,steps)


def initial_charge(fraction):
    velocity=np.sqrt(2/3)
    rho_phi=velocity**2/2
    rho_b=rho_phi*fraction/(1-fraction)
    charge_squared=2*rho_b*np.sqrt(1-Q0)
    momentum=velocity*(1+rho_b*D0)
    hubble=np.sqrt((rho_phi+rho_b)/3)
    # State: log(a), phi, scalar momentum, H, proper time, accumulated pressure work.
    return np.array([0.,0.,momentum,hubble,1.,0.]),charge_squared


def velocity_from_momentum(a,phi,momentum,charge_squared):
    coupling=D0*np.exp(BETA*phi)
    if charge_squared==0:
        return momentum/a**3
    low=0.; high=min(momentum/a**3,1/np.sqrt(coupling))
    for _ in range(60):
        v=(low+high)/2
        c=1-coupling*v*v
        recovered=a**3*v+charge_squared*coupling*v/(2*a**3*np.sqrt(c))
        if recovered<momentum:
            low=v
        else:
            high=v
    return (low+high)/2


def state_quantities(y,charge_squared):
    n,phi,momentum,hubble,time,work=y
    a=np.exp(n); coupling=D0*np.exp(BETA*phi)
    v=velocity_from_momentum(a,phi,momentum,charge_squared)
    q=coupling*v*v
    if q>=1:
        raise ValueError('Degenerate matter metric.')
    c=1-q
    rho_b=charge_squared/(2*a**6*np.sqrt(c))
    rho_phi=.5*v*v
    return v,q,rho_phi+rho_b,rho_phi+c*rho_b,rho_b,coupling


def charge_rhs(y,charge_squared,reciprocal_source=True):
    n,phi,momentum,hubble,time,work=y
    v,q,rho,pressure,rho_b,_=state_quantities(y,charge_squared)
    pdot=np.exp(3*n)*rho_b*BETA*q/2 if reciprocal_source else 0.
    return np.array([1.,v/hubble,pdot/hubble,-(rho+pressure)/(2*hubble),
                     1/hubble,3*np.exp(3*n)*pressure])


@lru_cache(maxsize=8)
def solve_charge(fraction=.1,end=1.,steps=2048,reciprocal_source=True):
    initial,charge_squared=initial_charge(fraction)
    path=integrate(lambda y:charge_rhs(y,charge_squared,reciprocal_source),initial,end,steps)
    return path,charge_squared


def routhian(a,phi,v,charge_squared):
    c=1-D0*np.exp(BETA*phi)*v*v
    return a**3*v*v/2-charge_squared*np.sqrt(c)/(2*a**3)


def closure_audit(reciprocal_source=True):
    path,charge=solve_charge(reciprocal_source=reciprocal_source)
    quantities=np.array([state_quantities(y,charge) for y in path])
    rho=quantities[:,2]
    constraint=3*path[:,3]**2-rho
    work=np.exp(3*path[:,0])*rho+path[:,5]-rho[0]
    return {'constraint_max':float(np.max(abs(constraint))),
            'pressure_work_residual_max':float(np.max(abs(work))),
            'final_fraction':float(quantities[-1,4]/rho[-1]),
            'final_q':float(quantities[-1,1]),
            'final_proper_time':float(path[-1,4])}


def report():
    rows=[]
    for fraction in (0.,.01,.1,.5):
        path=solve_phase(fraction)
        rows.append({'initial_matter_fraction':fraction,'final_log_scale_factor':100.,
                     'final_matter_fraction':float(path[-1,0]),'final_q':float(path[-1,1]),
                     'final_matter_metric_cone_ratio':float(np.sqrt(1-path[-1,1])),
                     'minimum_matter_fraction':float(np.min(path[:,0])),
                     'minimum_probe_lapse_squared':float(np.min(1-path[:,1])),
                     'final_log_H_over_initial_H':float(path[-1,2])})
    return {'round':336,
            'scope':'A new stipulated field-dependent coupling D(phi)=0.6*exp(sqrt(6)*phi) in the generally covariant classical Einstein+canonical phi+disformal massless chi action. The zero-chi branch retains q=D*phidot^2=0.4 forever. Nonzero homogeneous chi shift charge is included together with phi exchange and Einstein backreaction. Analytic attraction is restricted to expanding, positive-phidot data 0<f<=1/2, 0<q<=0.4; this is not a proof of full spatial/quantum stability or universal coupling for all matter.',
            'critical_vacuum_branch':{'beta':float(BETA),'q':Q0,'relative_speed':float(np.sqrt(1-Q0))},
            'cone_ratio_caveat':'sqrt(1-q) compares the two prescribed metrics. With nonzero chi background the coupled phi-chi perturbation eigenmodes must be analyzed separately.',
            'finite_occupation_phase_evolution':rows,
            'full_canonical_charge_and_gravity_check':closure_audit(),
            'missing_reciprocal_phi_source_counterexample':closure_audit(False),
            'analytic_bounds':{'uniform_q_upper':.4,'uniform_fraction_upper':.5,
                               'fraction_lower_factor':float(np.exp(-40*.4/3)),
                               'q_log_derivative_upper_coefficient':-.45},
            'next_interface':'Test whether the finite-occupation attraction survives nonzero spatial gradients and multiple species, and establish an independently motivated condition on coupling functions. The zero-occupation and late-time limits do not commute in this restricted branch.'}


class Checks(unittest.TestCase):
    def test_01_critical_coupling_defeats_vacuum_background_dilution(self):
        for t in (1.,10.,1e4):
            phi=np.sqrt(2/3)*np.log(t); v=np.sqrt(2/3)/t
            q=D0*np.exp(BETA*phi)*v*v
            self.assertAlmostEqual(q,Q0,places=14)

    def test_02_routhian_momentum_and_scalar_source_are_independently_varied(self):
        y,charge=initial_charge(.1); a=1.; phi=0.; v=np.sqrt(2/3); h=1e-5
        momentum=(routhian(a,phi,v+h,charge)-routhian(a,phi,v-h,charge))/(2*h)
        source=(routhian(a,phi+h,v,charge)-routhian(a,phi-h,v,charge))/(2*h)
        self.assertLess(abs(momentum-y[2]),2e-10)
        expected=charge*BETA*Q0/(4*np.sqrt(1-Q0))
        self.assertLess(abs(source-expected),2e-10)

    def test_03_momentum_velocity_map_is_regular_and_invertible(self):
        y,charge=initial_charge(.1)
        v=velocity_from_momentum(1.,0.,y[2],charge)
        self.assertAlmostEqual(v,np.sqrt(2/3),places=14)
        derivative=1+charge*D0/(2*(1-Q0)**1.5)
        self.assertGreater(derivative,1.)

    def test_04_initial_gravity_constraint_and_phase_coordinates(self):
        for fraction in (.01,.1,.5):
            y,charge=initial_charge(fraction)
            v,q,rho,p,rho_b,_=state_quantities(y,charge)
            self.assertAlmostEqual(3*y[3]**2,rho,places=14)
            self.assertAlmostEqual(rho_b/rho,fraction,places=14)
            self.assertAlmostEqual(q,Q0,places=14)

    def test_05_full_charge_equations_reproduce_independent_phase_reduction(self):
        path,charge=solve_charge()
        physical=np.array([state_quantities(y,charge) for y in path])
        phase=solve_phase(.1,end=1.,steps=2048)
        np.testing.assert_allclose(physical[:,4]/physical[:,2],phase[:,0],atol=3e-10,rtol=0.)
        np.testing.assert_allclose(physical[:,1],phase[:,1],atol=3e-10,rtol=0.)
        np.testing.assert_allclose(np.log(path[:,3]/path[0,3]),phase[:,2],atol=3e-10,rtol=0.)

    def test_06_full_energy_and_gravity_ledgers_close(self):
        result=closure_audit()
        self.assertLess(result['constraint_max'],2e-10)
        self.assertLess(result['pressure_work_residual_max'],2e-9)

    def test_07_reciprocal_source_cannot_be_omitted(self):
        result=closure_audit(False)
        self.assertGreater(result['constraint_max'],1e-4)
        self.assertGreater(result['pressure_work_residual_max'],1e-3)

    def test_08_invariant_rectangle_and_differential_bounds(self):
        for f in (.001,.01,.1,.5):
            for q in (.001,.1,.4):
                fp,qp,_=phase_rhs([f,q,0.])
                self.assertLess(fp,0.); self.assertLess(qp,0.)
                self.assertLessEqual(qp/q,-.45*f+1e-15)
                self.assertGreaterEqual(fp/f,-6*q*f-1e-15)
                self.assertLessEqual((fp/f)/qp,40/3+1e-13)

    def test_09_nonzero_occupation_lowers_mismatch(self):
        for f in (.01,.1,.5):
            path=solve_phase(f)
            self.assertTrue(np.all(np.diff(path[:,0])<=1e-13))
            self.assertTrue(np.all(np.diff(path[:,1])<=1e-13))
            self.assertGreaterEqual(np.min(path[:,0]),f*np.exp(-40*Q0/3))
            self.assertLess(path[-1,1],Q0)
            self.assertGreaterEqual(np.min(1-path[:,1]),.6-1e-14)

    def test_10_vacuum_limit_stays_nonuniversal(self):
        path=solve_phase(0.)
        np.testing.assert_array_equal(path[:,0],0.)
        np.testing.assert_array_equal(path[:,1],Q0)

    def test_11_step_error_is_below_backreaction_effect(self):
        fine=solve_phase(.1); coarse=solve_phase(.1,steps=2048)
        self.assertLess(np.max(abs(fine[::2,:2]-coarse[:,:2])),1e-8)
        full,_=solve_charge(); coarse_full,_=solve_charge(steps=1024)
        self.assertLess(np.max(abs(full[::2,:3]-coarse_full[:,:3])),1e-8)

    def test_12_small_occupation_does_not_give_a_uniform_convergence_rate(self):
        finite=solve_phase(.1)[-1,1]
        tiny=solve_phase(.0001)[-1,1]
        self.assertLess(finite,.001)
        self.assertGreater(tiny,.35)


if __name__=='__main__':
    main(__name__,'disformal_occupied_attractor_audit',report)
