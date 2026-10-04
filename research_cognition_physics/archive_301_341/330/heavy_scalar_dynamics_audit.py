"""Round 330: potential and induced kinetic matching for prepared heavy/light dynamics."""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main
from conformal_probe_motion_audit import integrate
from heavy_scalar_potential_audit import ALPHA, MASSES, moments, full_potential

J_HESSIAN=np.diag(2*ALPHA*MASSES**2)


def effective(chi,mass,order=2,kinetic=True):
    v0,j,b=moments(chi)
    gj=2*ALPHA*MASSES**2*chi
    gb=2*ALPHA**2*MASSES**2*chi
    metric=np.eye(2)
    value=v0; force=MASSES**2*chi
    if order>=1:
        value-=j*j/(2*mass**2)
        force=force-j*gj/mass**2
    if order>=2:
        value+=b*j*j/mass**4
        force=force+(gb*j*j+2*b*j*gj)/mass**4
        if kinetic:
            metric=metric+np.outer(gj,gj)/mass**4
    return value,force,metric,gj


def reduced_rhs(y,mass,order=2,kinetic=True):
    chi,v=y[:2],y[2:]
    _,force,metric,gj=effective(chi,mass,order,kinetic)
    correction=(gj*(v@J_HESSIAN@v)/mass**4) if order>=2 and kinetic else np.zeros(2)
    acceleration=np.linalg.solve(metric,-force-correction)
    return np.r_[v,acceleration]


def prepared_heavy_data(chi,v,mass):
    _,j,b=moments(chi)
    jd=2*np.sum(ALPHA*MASSES**2*chi*v)
    jdd=2*np.sum(ALPHA*MASSES**2*(v*v-MASSES**2*chi*chi))
    jddd=-8*np.sum(ALPHA*MASSES**4*chi*v)
    bd=2*np.sum(ALPHA**2*MASSES**2*chi*v)
    return (-j/mass**2+(jdd+2*b*j)/mass**4,
            -jd/mass**2+(jddd+2*bd*j+2*b*jd)/mass**4)


@lru_cache(maxsize=32)
def solve(mass=8.,order=None,steps=2048,kinetic=True):
    chi=np.array([.4,.3]); v=np.array([.1,-.08])
    if order is not None:
        return integrate(lambda y:reduced_rhs(y,mass,order,kinetic),np.r_[chi,v],6.,steps)
    phi,w=prepared_heavy_data(chi,v,mass)
    def rhs(y):
        q,velocity=y[:3],y[3:]
        m2=MASSES**2*np.exp(2*ALPHA*q[0])
        force=np.r_[mass**2*q[0]+np.sum(ALPHA*m2*q[1:]**2),m2*q[1:]]
        return np.r_[velocity,-force]
    return integrate(rhs,np.r_[phi,chi,w,v],6.,steps)


def energy(y,mass,order=None,kinetic=True):
    if order is None:
        return .5*y[3:]@y[3:]+full_potential(y[0],y[1:3],mass)
    value,_,metric,_=effective(y[:2],mass,order,kinetic)
    return .5*y[2:]@metric@y[2:]+value


@lru_cache(maxsize=8)
def audit(mass=8.):
    full=solve(mass)
    target=full[:,[1,2,4,5]]
    errors=[]
    for order in (0,1,2):
        reduced=solve(mass,order)
        errors.append(float(np.max(abs(target-reduced))))
    missing=solve(mass,2,kinetic=False)
    full_energies=np.array([energy(y,mass) for y in full])
    reduced=solve(mass,2)
    reduced_energies=np.array([energy(y,mass,2) for y in reduced])
    return {'mass':mass,'light_phase_space_errors_orders_0_1_2':errors,
            'order_2_missing_kinetic_error':float(np.max(abs(target-missing))),
            'full_energy_drift':float(np.max(abs(full_energies-full_energies[0]))),
            'effective_energy_drift':float(np.max(abs(reduced_energies-reduced_energies[0]))),
            'initial_energy_matching_error':float(abs(full_energies[0]-reduced_energies[0]))}


def occupied_heavy_example(mass,stored_energy=.03):
    amplitude=np.sqrt(2*stored_energy)/mass
    times=np.linspace(0.,2*np.pi/mass,513)
    phi=amplitude*np.cos(mass*times)
    velocity=-mass*amplitude*np.sin(mass*times)
    density=.5*velocity**2+.5*mass**2*phi**2
    return {'mass':mass,'amplitude':float(amplitude),
            'energy_density_min':float(np.min(density)),
            'energy_density_max':float(np.max(density))}


def report():
    return {'round':330,
            'scope':'Tree-level inverse-heavy-mass matching through mass^(-4), including the induced light-field kinetic metric, checked on homogeneous flat-space classical evolution for t in [0,6]. Heavy initial data are prepared through the same order. The initial-state counterexample is a free occupied massive mode with fixed energy. Curved dynamical metric matching and loop corrections remain open.',
            'heavy_mass_sequence':[audit(m) for m in (4.,8.,16.)],
            'occupied_modes_not_in_the_vacuum_effective_branch':[occupied_heavy_example(m) for m in (4.,8.,16.)],
            'next_interface':'Match the effective stress and dynamical geometry together, retaining initial heavy-sector energy when occupied. Separate matter matching from a derivation of the Einstein action.'}


class Checks(unittest.TestCase):
    def test_01_effective_potential_gradient(self):
        chi=np.array([.4,-.3]); h=1e-5
        numerical=[(effective(chi+h*d,8.)[0]-effective(chi-h*d,8.)[0])/(2*h) for d in np.eye(2)]
        np.testing.assert_allclose(numerical,effective(chi,8.)[1],atol=2e-11)

    def test_02_induced_metric_is_positive(self):
        for chi in (np.array([.4,.3]),np.array([2.,-1.])):
            _,_,metric,gj=effective(chi,4.)
            np.testing.assert_allclose(np.linalg.eigvalsh(metric),[1.,1.+gj@gj/4.**4],atol=2e-15)

    def test_03_metric_connection_matches_direct_derivatives(self):
        chi=np.array([.4,-.3]); mass=4.; h=1e-5
        _,_,metric,gj=effective(chi,mass)
        dm=np.array([(effective(chi+h*d,mass)[2]-effective(chi-h*d,mass)[2])/(2*h) for d in np.eye(2)])
        gamma=np.zeros((2,2,2)); inverse=np.linalg.inv(metric)
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    gamma[a,b,c]=sum(inverse[a,d]*(dm[b,d,c]+dm[c,d,b]-dm[d,b,c])/2 for d in range(2))
        analytic=np.einsum('a,bc->abc',inverse@gj,J_HESSIAN)/mass**4
        np.testing.assert_allclose(gamma,analytic,atol=2e-11)

    def test_04_effective_energy_is_conserved_by_its_equations(self):
        y=np.array([.4,.3,.1,-.08]); h=1e-5
        gradient=np.array([(energy(y+h*d,8.,2)-energy(y-h*d,8.,2))/(2*h) for d in np.eye(4)])
        self.assertLess(abs(gradient@reduced_rhs(y,8.)),2e-11)

    def test_05_full_and_effective_integrations_conserve_energy(self):
        r=audit(8.)
        self.assertLess(r['full_energy_drift'],2e-10)
        self.assertLess(r['effective_energy_drift'],2e-10)

    def test_06_matched_terms_improve_light_evolution(self):
        for mass in (4.,8.,16.):
            errors=audit(mass)['light_phase_space_errors_orders_0_1_2']
            self.assertLess(errors[1],errors[0])
            self.assertLess(errors[2],errors[1])

    def test_07_kinetic_term_cannot_be_dropped_at_the_same_order(self):
        r=audit(8.)
        self.assertGreater(r['order_2_missing_kinetic_error'],5*r['light_phase_space_errors_orders_0_1_2'][2])

    def test_08_inverse_mass_error_scaling(self):
        errors=[audit(m)['light_phase_space_errors_orders_0_1_2'] for m in (8.,16.)]
        for n,minimum_ratio in enumerate((3.,12.,40.)):
            self.assertGreater(errors[0][n]/errors[1][n],minimum_ratio)

    def test_09_time_discretization_is_below_matching_error(self):
        fine=solve(8.)[:,[1,2,4,5]]
        coarse=solve(8.,steps=1024)[:,[1,2,4,5]]
        numerical=float(np.max(abs(fine[::2]-coarse)))
        self.assertLess(numerical,audit(8.)['light_phase_space_errors_orders_0_1_2'][2]/100)

    def test_10_mass_gap_does_not_erase_occupied_energy(self):
        for mass in (4.,8.,16.):
            r=occupied_heavy_example(mass)
            self.assertAlmostEqual(r['energy_density_min'],.03,places=14)
            self.assertAlmostEqual(r['energy_density_max'],.03,places=14)
        self.assertEqual(occupied_heavy_example(8.)['amplitude'],2*occupied_heavy_example(16.)['amplitude'])


if __name__=='__main__':
    main(__name__,'heavy_scalar_dynamics_audit',report)
