"""Round 338: an open attraction region in coupling slope and a finite-time failure."""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main
from conformal_probe_motion_audit import integrate
from disformal_occupied_attractor_audit import phase_rhs as critical_rhs, initial_charge
from disformal_characteristics_audit import analytic_squares


def phase_rhs(y, beta):
    fraction,q,logh=y
    if not (0 <= fraction < 1 and 0 < q < 1):
        raise ValueError('Outside the positive homogeneous branch.')
    numerator=beta*np.sqrt(6*(1-fraction))-6+3*q*fraction/(1-fraction)
    denominator=1+q*fraction/(2*(1-fraction)*(1-q))
    qprime=q*numerator/denominator
    return np.array([fraction*(qprime/(2*(1-q))-3*q*fraction),
                     qprime,-3+1.5*q*fraction])


@lru_cache(maxsize=16)
def solve_phase(beta, fraction=.1, q=.4, end=100., steps=4096):
    return integrate(lambda y:phase_rhs(y,beta),np.array([fraction,q,0.]),end,steps)


def canonical_quantities(y, charge, beta):
    n,phi,momentum,hubble,time,work=y
    a=np.exp(n); coupling=.6*np.exp(beta*phi)
    low=0.; high=min(momentum/a**3,1/np.sqrt(coupling))
    for _ in range(60):
        v=(low+high)/2
        c=1-coupling*v*v
        recovered=a**3*v+charge*coupling*v/(2*a**3*np.sqrt(c))
        if recovered < momentum:
            low=v
        else:
            high=v
    v=(low+high)/2; q=coupling*v*v; c=1-q
    matter=charge/(2*a**6*np.sqrt(c))
    rho=.5*v*v+matter
    pressure=.5*v*v+c*matter
    return v,q,matter,rho,pressure


def canonical_rhs(y, charge, beta):
    n,phi,momentum,hubble,time,work=y
    v,q,matter,rho,pressure=canonical_quantities(y,charge,beta)
    return np.array([1.,v/hubble,np.exp(3*n)*matter*beta*q/(2*hubble),
                     -(rho+pressure)/(2*hubble),1/hubble,3*np.exp(3*n)*pressure])


@lru_cache(maxsize=8)
def solve_canonical(beta=10., end=.15, steps=2048):
    initial,charge=initial_charge(.1)
    path=integrate(lambda y:canonical_rhs(y,charge,beta),initial,end,steps)
    return path,charge


def full_audit():
    path,charge=solve_canonical()
    values=np.array([canonical_quantities(y,charge,10.) for y in path])
    phase=solve_phase(10.,end=.15,steps=2048)
    return {'maximum_fraction_difference':float(np.max(abs(values[:,2]/values[:,3]-phase[:,0]))),
            'maximum_q_difference':float(np.max(abs(values[:,1]-phase[:,1]))),
            'maximum_friedmann_residual':float(np.max(abs(3*path[:,3]**2-values[:,3]))),
            'maximum_pressure_work_residual':float(np.max(abs(np.exp(3*path[:,0])*values[:,3]+path[:,5]-values[0,3])))}


def instability_crossing(steps=4096):
    end=.15
    path=solve_phase(10.,end=end,steps=steps)
    r=path[:,0]*path[:,1]/(2*(1-path[:,0]))
    index=int(np.flatnonzero(r>=1)[0])
    weight=(1-r[index-1])/(r[index]-r[index-1])
    point=path[index-1]+weight*(path[index]-path[index-1])
    n=(index-1+weight)*end/steps
    return {'log_scale_factor_at_r_one':float(n),
            'fraction':float(point[0]),'q':float(point[1]),
            'matter_metric_lapse_squared':float(1-point[1]),
            'end_log_scale_factor':end,'end_fraction':float(path[-1,0]),
            'end_q':float(path[-1,1]),
            'end_phi_speed_squared':float(analytic_squares(*path[-1,:2])[0])}


def report():
    rows=[]
    for delta,fraction in ((-.01,.09),(0.,.1),(.01,.09),(.01,.1)):
        beta=(6+delta)/np.sqrt(6)
        path=solve_phase(beta,fraction)
        rows.append({'slope_offset_beta_sqrt6_minus6':delta,'beta':float(beta),
                     'initial_fraction':fraction,'end_log_scale_factor':100.,
                     'end_fraction':float(path[-1,0]),'end_q':float(path[-1,1]),
                     'minimum_fraction':float(np.min(path[:,0]))})
    return {'round':338,
            'scope':'Classical homogeneous zero-potential two-scalar disformal model with D(phi)=0.6 exp(beta phi). A bootstrap proof gives an open interval of slopes and initial matter fractions with stable principal cones approaching the Einstein cone. A beta=10 branch is proved to cross the gradient-stability boundary in finite N while the matter metric remains Lorentzian, with independently verified background evolution. This does not select the coupling from cognition or establish nonlinear, quantum or all-species universality.',
            'analytic_bootstrap':{'initial_fraction_interval':[.09,.1],'initial_q_upper':.4,
                                  'slope_offset_interval':[-.01,.01],
                                  'beta_interval':[(6-.01)/np.sqrt(6),(6+.01)/np.sqrt(6)],
                                  'temporary_fraction_lower':.02,
                                  'proved_fraction_lower':float(.09*np.exp(-1.4)),
                                  'q_log_derivative_bound':'q_prime/q <= -(9/8)*f',
                                  'log_fraction_by_q_upper':3.5,
                                  'uniform_q_decay_rate_in_N':.0225},
            'open_interval_trajectories':rows,
            'supercritical_crossing':instability_crossing(),
            'analytic_instability_deadline_in_N':float(4*np.log(25/3)),
            'crossing_step_refinement_error':abs(instability_crossing()['log_scale_factor_at_r_one']
                                                 -instability_crossing(2048)['log_scale_factor_at_r_one']),
            'independent_canonical_and_einstein_check':full_audit(),
            'next_interface':'Extend the simultaneous attraction and hyperbolicity criterion to finitely many independently occupied matter species; retain the action, initial-data and EFT applicability gaps.'}


class Checks(unittest.TestCase):
    def test_01_general_slope_reduces_to_frozen_critical_equations(self):
        for f,q in ((.01,.4),(.1,.1),(.5,.4)):
            np.testing.assert_allclose(phase_rhs([f,q,0.],np.sqrt(6)),critical_rhs([f,q,0.]),atol=8e-16)

    def test_02_zero_occupation_exact_solution_changes_with_slope(self):
        beta=(6-.1)/np.sqrt(6); end=3.
        path=solve_phase(beta,0.,end=end,steps=512)
        n=np.linspace(0,end,len(path))
        np.testing.assert_allclose(path[:,1],.4*np.exp(-.1*n),atol=2e-14,rtol=2e-14)
        self.assertGreater(phase_rhs([0.,.4,0.],4.)[1],0.)

    def test_03_bootstrap_differential_inequalities(self):
        for delta in (-.01,0.,.01):
            beta=(6+delta)/np.sqrt(6)
            for f in np.linspace(.02,.1,17):
                for q in np.linspace(.001,.4,17):
                    fp,qp,_=phase_rhs([f,q,0.],beta)
                    self.assertLessEqual(qp/q,-9*f/8+1e-14)
                    self.assertLess(fp,0.)
                    self.assertLessEqual((fp/f)/qp,3.5+1e-13)

    def test_04_lower_fraction_bound_closes_bootstrap(self):
        lower=.09*np.exp(-3.5*.4)
        self.assertGreater(lower,.02)
        self.assertLess((6-.01)/np.sqrt(6),np.sqrt(6))
        self.assertGreater((6+.01)/np.sqrt(6),np.sqrt(6))

    def test_05_open_interval_evolutions_obey_bounds_and_characteristics(self):
        for delta,f in ((-.01,.09),(0.,.1),(.01,.09),(.01,.1)):
            path=solve_phase((6+delta)/np.sqrt(6),f)
            n=np.linspace(0,100,len(path))
            self.assertTrue(np.all(path[:,0]>=f*np.exp(-1.4)-1e-12))
            self.assertTrue(np.all(path[:,1]<=.4*np.exp(-.0225*n)+1e-12))
            for fraction,q,_ in path[::64]:
                phi,theta=analytic_squares(fraction,q)
                self.assertGreaterEqual(phi,theta-1e-14)
                self.assertGreaterEqual(theta,.6-1e-14)

    def test_06_routhian_source_for_changed_slope_is_independently_varied(self):
        beta=4.; y,charge=initial_charge(.1); v=np.sqrt(2/3); h=1e-5
        def lagrangian(phi):
            return v*v/2-charge*np.sqrt(1-.6*np.exp(beta*phi)*v*v)/2
        numerical=(lagrangian(h)-lagrangian(-h))/(2*h)
        expected=canonical_rhs(y,charge,beta)[2]*y[3]
        self.assertLess(abs(numerical-expected),2e-10)

    def test_07_full_charge_gravity_equations_match_independent_phase_flow(self):
        result=full_audit()
        self.assertLess(result['maximum_fraction_difference'],2e-8)
        self.assertLess(result['maximum_q_difference'],2e-8)
        self.assertLess(result['maximum_friedmann_residual'],2e-9)
        self.assertLess(result['maximum_pressure_work_residual'],2e-9)

    def test_08_finite_evolution_crosses_gradient_boundary_before_metric_degeneracy(self):
        event=instability_crossing()
        self.assertGreater(event['matter_metric_lapse_squared'],.003)
        self.assertLess(event['log_scale_factor_at_r_one'],.15)
        self.assertLess(event['end_phi_speed_squared'],-.0001)
        self.assertLess(event['end_q'],1.)

    def test_09_crossing_and_trajectories_converge_with_step_refinement(self):
        self.assertLess(abs(instability_crossing()['log_scale_factor_at_r_one']
                            -instability_crossing(2048)['log_scale_factor_at_r_one']),1e-7)
        beta=(6+.01)/np.sqrt(6)
        fine=solve_phase(beta); coarse=solve_phase(beta,steps=2048)
        self.assertLess(np.max(abs(fine[::2,:2]-coarse[:,:2])),1e-8)

    def test_10_canonical_integration_error_is_smaller_than_instability_effect(self):
        fine,_=solve_canonical(); coarse,_=solve_canonical(steps=1024)
        self.assertLess(np.max(abs(fine[::2,:3]-coarse[:,:3])),1e-7)

    def test_11_supercritical_analytic_growth_bound_before_stability_exit(self):
        for q in np.linspace(.4,.999,19):
            upper=2/(q+2)
            for f in np.linspace(.1,upper-1e-7,19):
                fp,qp,_=phase_rhs([f,q,0.],10.)
                self.assertGreater(qp,0.)
                self.assertGreaterEqual(fp/f,.25)
                r=q*f/(2*(1-f))
                self.assertGreater(qp/q+fp/(f*(1-f)),0.)
                self.assertLess(r,1.)


if __name__=='__main__':
    main(__name__,'disformal_slope_robustness_audit',report)
