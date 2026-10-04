"""Round 331: match prepared heavy/light matter and dynamical FLRW geometry."""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main
from conformal_probe_motion_audit import integrate
from heavy_scalar_dynamics_audit import effective, prepared_heavy_data, J_HESSIAN
from heavy_scalar_potential_audit import ALPHA, MASSES, full_potential


def initial(mass, order=None, kinetic=True):
    chi=np.array([.4,.3]); velocity=np.array([.1,-.08])
    if order is not None:
        value,_,metric,_=effective(chi,mass,order,kinetic)
        rho=.5*velocity@metric@velocity+value
        return np.r_[chi,velocity,0.,np.sqrt(rho/3),0.]
    phi,flat_w=prepared_heavy_data(chi,velocity,mass)
    c=12*np.sum(ALPHA*MASSES**2*velocity**2)/mass**4
    rho_base=.5*velocity@velocity+full_potential(phi,chi,mass)+.5*flat_w**2
    coefficient=3-.5*c*c
    if coefficient<=0:
        raise ValueError('Prepared-data expansion is outside the tested mass window.')
    linear=flat_w*c
    hubble=2*rho_base/(linear+np.sqrt(linear**2+4*coefficient*rho_base))
    w=flat_w-c*hubble
    return np.r_[phi,chi,w,velocity,0.,hubble,0.]


def densities(y,mass,order=None,kinetic=True):
    n=3 if order is None else 2
    q,v=y[:n],y[n:2*n]
    if order is None:
        potential=full_potential(q[0],q[1:],mass)
        kinetic_energy=.5*v@v
    else:
        potential,_,metric,_=effective(q,mass,order,kinetic)
        kinetic_energy=.5*v@metric@v
    return kinetic_energy+potential,kinetic_energy-potential


def rhs(y,mass,order=None,kinetic=True):
    n=3 if order is None else 2
    q,v=y[:n],y[n:2*n]
    loga,hubble=y[-3:-1]
    rho,pressure=densities(y,mass,order,kinetic)
    if order is None:
        m2=MASSES**2*np.exp(2*ALPHA*q[0])
        force=np.r_[mass**2*q[0]+np.sum(ALPHA*m2*q[1:]**2),m2*q[1:]]
        acceleration=-3*hubble*v-force
    else:
        _,force,metric,gj=effective(q,mass,order,kinetic)
        connection=(gj*(v@J_HESSIAN@v)/mass**4) if order>=2 and kinetic else np.zeros(2)
        acceleration=-3*hubble*v-np.linalg.solve(metric,force+connection)
    return np.r_[v,acceleration,hubble,-.5*(rho+pressure),3*hubble*np.exp(3*loga)*pressure]


@lru_cache(maxsize=24)
def solve(mass=8.,order=None,steps=2048,kinetic=True):
    return integrate(lambda y:rhs(y,mass,order,kinetic),
                     initial(mass,order,kinetic),6.,steps)


def ledgers(path,mass,order=None,kinetic=True):
    stress=np.array([densities(y,mass,order,kinetic) for y in path])
    constraint=3*path[:,-2]**2-stress[:,0]
    work=np.exp(3*path[:,-3])*stress[:,0]+path[:,-1]-stress[0,0]
    return stress,{'constraint_max':float(np.max(abs(constraint))),
                   'pressure_work_residual_max':float(np.max(abs(work)))}


@lru_cache(maxsize=8)
def audit(mass=8.):
    full=solve(mass)
    stress,full_ledger=ledgers(full,mass)
    rows=[]
    for order in (0,1,2):
        reduced=solve(mass,order)
        reduced_stress,ledger=ledgers(reduced,mass,order)
        rows.append({'order':order,
                     'light_phase_space_error':float(np.max(abs(full[:,[1,2,4,5]]-reduced[:,:4]))),
                     'scale_factor_error':float(np.max(abs(np.exp(full[:,-3])-np.exp(reduced[:,-3])))),
                     'hubble_error':float(np.max(abs(full[:,-2]-reduced[:,-2]))),
                     'density_error':float(np.max(abs(stress[:,0]-reduced_stress[:,0]))),
                     'pressure_error':float(np.max(abs(stress[:,1]-reduced_stress[:,1]))),
                     'initial_hubble_difference':float(abs(full[0,-2]-reduced[0,-2])),
                     **ledger})
    missing=solve(mass,2,kinetic=False)
    return {'mass':mass,'full_final_scale_factor':float(np.exp(full[-1,-3])),
            'full_ledger':full_ledger,'effective_orders':rows,
            'order_2_missing_kinetic_scale_factor_error':float(np.max(abs(np.exp(full[:,-3])-np.exp(missing[:,-3]))))}


def minisuperspace_matter(lapse,scale,chi,coordinate_velocity,mass):
    potential,_,metric,_=effective(chi,mass)
    return scale**3*(coordinate_velocity@metric@coordinate_velocity/(2*lapse)-lapse*potential)


def report():
    return {'round':331,
            'scope':'Classical tree-level effective matter through M^(-4), coupled to the stipulated spatially flat Einstein FLRW branch with kappa=1 and Lambda=0. Each model satisfies its own initial Friedmann constraint for identical initial light data and a=1. Fixed time interval [0,6]; this is a matching audit, not a derivation of the gravitational action.',
            'mass_sequence':[audit(m) for m in (4.,8.,16.)],
            'next_interface':'An occupied heavy mode contributes an independent stress source; test its averaged gravitational effect instead of deleting it.'}


class Checks(unittest.TestCase):
    def test_01_prepared_initial_data_obey_friedmann_constraint(self):
        for mass in (4.,8.,16.):
            for order in (None,0,1,2):
                y=initial(mass,order)
                self.assertLess(abs(3*y[-2]**2-densities(y,mass,order)[0]),1e-15)

    def test_02_common_light_initial_data_are_preserved(self):
        for mass in (4.,8.,16.):
            full=initial(mass)
            for order in (0,1,2):
                reduced=initial(mass,order)
                np.testing.assert_array_equal(full[[1,2,4,5]],reduced[:4])
                self.assertEqual(full[-3],reduced[-3])

    def test_03_lapse_and_scale_variations_recover_stress(self):
        chi=np.array([.4,.3]); velocity=np.array([.1,-.08]); a=1.2; step=1e-5
        y=np.r_[chi,velocity,np.log(a),.2,0.]
        rho,pressure=densities(y,8.,2)
        lapse_derivative=(minisuperspace_matter(1+step,a,chi,velocity,8.)-minisuperspace_matter(1-step,a,chi,velocity,8.))/(2*step)
        scale_derivative=(minisuperspace_matter(1,a+step,chi,velocity,8.)-minisuperspace_matter(1,a-step,chi,velocity,8.))/(2*step)
        self.assertLess(abs(-lapse_derivative/a**3-rho),2e-10)
        self.assertLess(abs(scale_derivative/(3*a*a)-pressure),2e-10)

    def test_04_reduced_continuity_equation_is_local(self):
        y=initial(8.,2); step=1e-5
        gradient=np.array([(densities(y+step*d,8.,2)[0]-densities(y-step*d,8.,2)[0])/(2*step) for d in np.eye(7)])
        rho,p=densities(y,8.,2)
        self.assertLess(abs(gradient@rhs(y,8.,2)+3*y[-2]*(rho+p)),2e-10)

    def test_05_full_and_effective_constraints_and_work_close(self):
        for mass in (4.,8.,16.):
            r=audit(mass)
            for ledger in [r['full_ledger']]+r['effective_orders']:
                self.assertLess(ledger['constraint_max'],2e-10)
                self.assertLess(ledger['pressure_work_residual_max'],2e-9)

    def test_06_matching_improves_geometry_and_matter(self):
        for mass in (4.,8.,16.):
            rows=audit(mass)['effective_orders']
            for key in ('light_phase_space_error','scale_factor_error','hubble_error','density_error'):
                self.assertLess(rows[1][key],rows[0][key])
                self.assertLess(rows[2][key],rows[1][key])

    def test_07_geometry_residual_has_expected_mass_order(self):
        small=audit(8.)['effective_orders']; large=audit(16.)['effective_orders']
        for order,minimum in enumerate((3.,10.,35.)):
            self.assertGreater(small[order]['scale_factor_error']/large[order]['scale_factor_error'],minimum)

    def test_08_integrator_error_is_below_physical_matching_residual(self):
        for mass in (8.,16.):
            for order in (None,2):
                fine=solve(mass,order); coarse=solve(mass,order,steps=1024)
                error=float(np.max(abs(np.exp(fine[::2,-3])-np.exp(coarse[:,-3]))))
                self.assertLess(error,audit(mass)['effective_orders'][2]['scale_factor_error']/100)

    def test_09_flat_preparation_requires_hubble_correction(self):
        y=initial(8.); phi,flat_w=prepared_heavy_data(y[1:3],y[4:6],8.)
        expected=-12*y[-2]*np.sum(ALPHA*MASSES**2*y[4:6]**2)/8.**4
        self.assertAlmostEqual(y[0],phi,places=15)
        self.assertAlmostEqual(y[3]-flat_w,expected,places=15)
        self.assertGreater(abs(expected),1e-7)

    def test_10_dropping_induced_kinetic_term_worsens_geometry(self):
        r=audit(8.)
        self.assertGreater(r['order_2_missing_kinetic_scale_factor_error'],3*r['effective_orders'][2]['scale_factor_error'])


if __name__=='__main__':
    main(__name__,'heavy_scalar_geometry_audit',report)
