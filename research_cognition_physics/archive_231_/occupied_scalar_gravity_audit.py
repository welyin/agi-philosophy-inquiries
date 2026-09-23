"""Round 332: retain an occupied heavy mode as an averaged gravitational source."""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main
from conformal_probe_motion_audit import integrate

RHO0=.03
H0=np.sqrt(RHO0/3)
END=10.
PHASES=(0.,np.pi/4,np.pi/2,3*np.pi/4)


def initial(mass,phase=0.):
    amplitude=np.sqrt(2*RHO0)
    return np.array([amplitude*np.cos(phase)/mass,-amplitude*np.sin(phase),0.,H0,0.,0.,0.,0.])


def stress(y,mass):
    phi,v=y[:2]
    return .5*(v*v+mass**2*phi*phi),.5*(v*v-mass**2*phi*phi)


def rhs(y,mass):
    phi,v,loga,hubble=y[:4]
    rho,p=stress(y,mass)
    # Final coordinates integrate pressure work, pressure, density, and H*phi*phidot.
    return np.array([v,-3*hubble*v-mass**2*phi,hubble,-.5*v*v,
                     3*hubble*np.exp(3*loga)*p,p,rho,hubble*phi*v])


@lru_cache(maxsize=24)
def solve(mass=8.,phase=0.,steps=4096):
    return integrate(lambda y:rhs(y,mass),initial(mass,phase),END,steps)


def dust(times):
    z=1+1.5*H0*times
    return z**(2/3),H0/z,RHO0/z**2


def corrected_comoving_energy(y,mass):
    rho,_=stress(y,mass)
    return np.exp(3*y[2])*(rho+1.5*y[3]*y[0]*y[1])


@lru_cache(maxsize=16)
def audit(mass=8.,phase=0.):
    path=solve(mass,phase)
    times=np.linspace(0.,END,len(path))
    a=np.exp(path[:,2]); hubble=path[:,3]
    rho,p=np.array([stress(y,mass) for y in path]).T
    ad,hd,_=dust(times)
    comoving=a**3*rho
    corrected=a**3*(rho+1.5*hubble*path[:,0]*path[:,1])
    virial=2*path[:,5]-(path[:,0]*path[:,1]-path[0,0]*path[0,1])-3*path[:,7]
    # Exact virial identity and |phi*phidot| <= rho/M give this finite-window bound.
    bound=((rho[0]+rho[-1])/(2*mass)+1.5*H0*path[-1,6]/mass)/path[-1,6]
    return {'mass':mass,'phase':float(phase),'initial_H_over_M':float(H0/mass),
            'scale_factor_error':float(np.max(abs(a-ad))),
            'hubble_error':float(np.max(abs(hubble-hd))),
            'final_scale_factor':float(a[-1]),
            'delete_source_scale_factor_error':float(a[-1]-1),
            'comoving_energy_relative_drift':float(np.max(abs(comoving/comoving[0]-1))),
            'corrected_comoving_energy_relative_drift':float(np.max(abs(corrected/corrected[0]-1))),
            'integrated_pressure_over_integrated_density':float(path[-1,5]/path[-1,6]),
            'finite_window_pressure_ratio_bound':float(bound),
            'max_instantaneous_abs_pressure_over_density':float(np.max(abs(p/rho))),
            'constraint_max':float(np.max(abs(3*hubble**2-rho))),
            'pressure_work_residual_max':float(np.max(abs(comoving+path[:,4]-RHO0))),
            'integrated_virial_residual_max':float(np.max(abs(virial)))}


def report():
    rows=[audit(m,phase) for m in (4.,8.,16.) for phase in PHASES]
    return {'round':332,
            'scope':'The exactly invariant chi_1=chi_2=0 sector of the stipulated classical Einstein-scalar theory, with a quadratic occupied heavy field, kappa=1, flat FLRW and Lambda=0. Fixed initial energy 0.03, four phases, t in [0,10]. Phase-averaged dust is a long-time/low-frequency approximation, not an exact pointwise stress replacement or a derivation of cosmological dark matter.',
            'cases':rows,
            'phase_envelopes':[{'mass':m,
                'maximum_scale_factor_error':max(r['scale_factor_error'] for r in rows if r['mass']==m),
                'maximum_absolute_averaged_pressure_ratio':max(abs(r['integrated_pressure_over_integrated_density']) for r in rows if r['mass']==m),
                'maximum_corrected_comoving_drift':max(r['corrected_comoving_energy_relative_drift'] for r in rows if r['mass']==m)}
                for m in (4.,8.,16.)],
            'rejected_stronger_hypothesis':'A phase-independent averaged-pressure ratio below 3 percent already at M=4 is false in this finite window. Use the analytic finite-window bound and keep the phase dependence.',
            'dust_final_scale_factor':float(dust(np.array([END]))[0][0]),
            'next_interface':'Test spatial perturbations and the wavelength window in which the retained heavy source behaves like dust; homogeneous expansion alone does not establish clustering or observed dark matter.'}


class Checks(unittest.TestCase):
    def test_01_all_initial_phases_obey_the_same_constraint(self):
        for m in (4.,8.,16.):
            for phase in PHASES:
                y=initial(m,phase)
                self.assertAlmostEqual(stress(y,m)[0],RHO0,places=15)
                self.assertAlmostEqual(3*y[3]**2,RHO0,places=15)

    def test_02_dust_solution_obeys_constraint_and_comoving_energy(self):
        a,h,rho=dust(np.linspace(0.,END,101))
        np.testing.assert_allclose(3*h*h,rho,atol=2e-17)
        np.testing.assert_allclose(a**3*rho,RHO0,atol=2e-17)

    def test_03_scalar_continuity_identity(self):
        y=initial(8.,.43); derivative=rhs(y,8.); rho,p=stress(y,8.)
        self.assertLess(abs(8.**2*y[0]*derivative[0]+y[1]*derivative[1]+3*y[3]*(rho+p)),1e-15)

    def test_04_integrated_virial_identity(self):
        for m in (4.,8.,16.):
            for phase in PHASES:
                self.assertLess(audit(m,phase)['integrated_virial_residual_max'],2e-8)

    def test_05_corrected_action_variable_has_exact_derivative(self):
        y=initial(8.,.43); step=1e-6
        gradient=np.array([(corrected_comoving_energy(y+step*d,8.)-corrected_comoving_energy(y-step*d,8.))/(2*step) for d in np.eye(8)])
        derivative=rhs(y,8.)
        expected=1.5*np.exp(3*y[2])*derivative[3]*y[0]*y[1]
        self.assertLess(abs(gradient@derivative-expected),2e-11)

    def test_06_full_gravity_and_pressure_work_close(self):
        for m in (4.,8.,16.):
            for phase in PHASES:
                r=audit(m,phase)
                self.assertLess(r['constraint_max'],2e-9)
                self.assertLess(r['pressure_work_residual_max'],2e-8)

    def test_07_retaining_occupied_energy_improves_geometry(self):
        for m in (4.,8.,16.):
            for phase in PHASES:
                r=audit(m,phase)
                self.assertLess(r['scale_factor_error'],r['delete_source_scale_factor_error']/30)

    def test_08_phase_envelope_decreases_in_the_adiabatic_limit(self):
        errors=[max(audit(m,p)['scale_factor_error'] for p in PHASES) for m in (4.,8.,16.)]
        self.assertGreater(errors[0]/errors[1],1.8)
        self.assertGreater(errors[1]/errors[2],1.8)

    def test_09_pressure_converges_after_averaging_not_pointwise(self):
        for m in (4.,8.,16.):
            self.assertAlmostEqual(audit(m,0.)['max_instantaneous_abs_pressure_over_density'],1.,places=14)
            for phase in PHASES:
                r=audit(m,phase)
                self.assertLess(abs(r['integrated_pressure_over_integrated_density']),r['finite_window_pressure_ratio_bound'])
        envelopes=[max(abs(audit(m,p)['integrated_pressure_over_integrated_density']) for p in PHASES) for m in (4.,8.,16.)]
        self.assertGreater(envelopes[0],.03)
        self.assertLess(envelopes[1],envelopes[0])
        self.assertLess(envelopes[2],envelopes[1])

    def test_10_corrected_action_removes_leading_fast_modulation(self):
        for m in (4.,8.,16.):
            for phase in PHASES:
                r=audit(m,phase)
                self.assertLess(r['corrected_comoving_energy_relative_drift'],r['comoving_energy_relative_drift']/5)

    def test_11_time_discretization_is_subdominant(self):
        fine=solve(16.,np.pi/4)
        coarse=solve(16.,np.pi/4,steps=2048)
        error=float(np.max(abs(np.exp(fine[::2,2])-np.exp(coarse[:,2]))))
        self.assertLess(error,audit(16.,np.pi/4)['scale_factor_error']/100)

    def test_12_light_data_alone_do_not_determine_heavy_energy(self):
        # The same chi=chidot=0 admits an empty branch H=0 or this occupied H=0.1 branch.
        for m in (4.,8.,16.):
            occupied=initial(m,0.)
            self.assertAlmostEqual(occupied[3],.1,places=15)
            self.assertGreater(occupied[3],0.)
            self.assertAlmostEqual(stress(occupied,m)[0],RHO0,places=15)


if __name__=='__main__':
    main(__name__,'occupied_scalar_gravity_audit',report)
