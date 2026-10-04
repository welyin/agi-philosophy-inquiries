"""645: one original mass dictionary, thermal reference and geometric sources.

Free continuum fermions in the declared 3+1 branch, not a computation of the
full graph Gauss Gibbs state. Quadrature has an infinite-domain variable map.
"""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_background_contact_matching as vacuum

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_thermal_background_source_results.json'
M = vacuum.matter.original.M
X0 = vacuum.U0.copy()
MASS, DIRAC = vacuum.canonical_mass_data(X0)
DEG = 4 * DIRAC
F0 = M - sum(X0) / 6
Q0 = np.sqrt(6) * np.arctanh(np.sqrt(sum(X0) / (6*M)))


@lru_cache(None)
def quadrature(order):
    v, w = np.polynomial.legendre.leggauss(order)
    u = (v + 1) / 2
    return u / (1-u), w / (2*(1-u)**2)


def thermal(T, masses=MASS, order=192, jets=True):
    y, w = quadrature(order)
    omega = np.sqrt(y[:, None]**2 + (masses[None, :]/T)**2)
    exp = np.exp(-omega)
    n = exp/(1+exp)
    measure = w[:, None]*y[:, None]**2
    factor = DEG/(2*np.pi**2)
    f = -T**4*np.sum(factor*np.sum(measure*np.log1p(exp), axis=0))
    rho = T**4*np.sum(factor*np.sum(measure*omega*n, axis=0))
    p = T**4/3*np.sum(factor*np.sum(measure*y[:, None]**2/omega*n, axis=0))
    theta = T**2*np.sum(factor*masses**2*np.sum(measure*n/omega, axis=0))
    ans = dict(f=float(f), rho=float(rho), p=float(p), theta=float(theta))
    if jets:
        assert np.min(masses)>0
        fz = factor*T*T/2*np.sum(measure*n/omega, axis=0)
        fzz = -factor/4*np.sum(measure*(n/omega**3+n*(1-n)/omega**2), axis=0)
        ans.update(fz=fz, fzz=fzz,
                   F1=float(np.sum(2*masses**2*fz)),
                   F2=float(np.sum(2*masses**2*fz+4*masses**4*fzz)))
    return ans


def scalar_gradient(T):
    t, grad_t, _, _ = vacuum.jordan_eigen_jets(X0)
    dz = grad_t/F0 - t[:, None]*vacuum.FI/F0**2
    return np.sum(thermal(T)['fz'][:, None]*dz, axis=0)


def thermodynamics_check():
    # Original 64-dimensional Nambu spectrum counts 32 positive excitations.
    phi = np.array([0., np.sqrt(X0[0]), 0., 0., np.sqrt(X0[1])])
    h, delta = vacuum.matter.mass_matrices(phi)
    bdg = np.block([[h, delta], [-delta.conj(), -h.T]])
    err_mass = float(np.max(abs(np.linalg.eigvalsh(bdg)[32:]
                               - np.sort(np.repeat(MASS, DEG.astype(int))))))
    assert err_mass < 3e-15 and sum(DEG)==32
    rows=[]
    for T in (.08, .2, .5):
        a = thermal(T); fine = thermal(T, order=256)
        err_quad = max(abs(a[k]-fine[k]) for k in ('f','rho','p','theta','F1','F2'))
        err_identity = max(abs(a['f']+a['p']), abs(a['theta']-a['rho']+3*a['p']),
                           abs(a['theta']-a['F1']))
        dt=T*2e-5
        df=(thermal(T+dt)['f']-thermal(T-dt)['f'])/(2*dt)
        err_rho=abs(a['f']-T*df-a['rho'])
        assert err_quad < 2e-12 and err_identity < 2e-13 and err_rho < 2e-9
        assert a['theta']>0 and a['rho']+a['p']>0
        rows.append(dict(T=T, **{k:a[k] for k in ('f','rho','p','theta','F1','F2')},
                         enthalpy=a['rho']+a['p'], quadrature_error=err_quad,
                         stress_identity_error=err_identity, temperature_derivative_error=err_rho))
    massless=thermal(.2, masses=np.zeros(5), jets=False)
    ideal=-sum(DEG)*7*np.pi**2*.2**4/720
    assert abs(massless['f']-ideal)<2e-15
    return dict(original_squared_background=X0.tolist(), F=float(F0),
                original_masses=MASS.tolist(), excitation_degeneracy=DEG.tolist(),
                original_BdG_mass_error=err_mass, rows=rows,
                massless_finite_temperature_limit_error=float(abs(massless['f']-ideal)),
                no_extra_Nambu_or_spin_multiplicity=True,
                independent_infinite_domain_orders=[192,256])


def hydrostatic_contact_check():
    T=.2; a=thermal(T); f1=a['F1']; f2=a['F2']
    b=1/(np.sqrt(6)*np.tanh(Q0/np.sqrt(6))); c=1/6
    expected=np.array([[b*b*f2+c*f1, b*(f2+f1)],
                       [b*(f2+f1), f2+f1]])
    contacts=f1*np.array([[c,b],[b,1]])
    assert np.max(abs(expected-f2*np.outer([b,1],[b,1])-contacts))<2e-17

    def density(z):
        q,sigma=z
        lam=np.sinh(q/np.sqrt(6))/np.sinh(Q0/np.sqrt(6))
        return np.exp(4*sigma)*thermal(T*np.exp(-sigma), masses=lam*MASS)['f']

    finite_error=0.
    for lam,sigma in ((.91,-.12),(1.04,.09),(1.13,-.07)):
        left=np.exp(4*sigma)*thermal(T*np.exp(-sigma), masses=lam*MASS)['f']
        right=thermal(T, masses=np.exp(sigma)*lam*MASS)['f']
        finite_error=max(finite_error,abs(left-right))
    fd=vacuum.hessian(density,np.array([Q0,0.]),step=1e-4)
    fd_error=float(np.max(abs(fd-expected)))
    step=2e-5
    lapse_fd=((1+step)*thermal(T/(1+step))['f']
              -(1-step)*thermal(T/(1-step))['f'])/(2*step)
    space_fd=(np.exp(3*step)*a['f']-np.exp(-3*step)*a['f'])/(2*step)
    weyl_fd=(density([Q0,step])-density([Q0,-step]))/(2*step)
    assert finite_error<2e-16 and fd_error<3e-9
    assert abs(lapse_fd-a['rho'])<2e-10 and abs(space_fd+3*a['p'])<2e-10
    assert abs(weyl_fd-a['theta'])<2e-10
    # Check one local static Ward direction with gradients of lapse and scalars.
    dlogN=.13; dx=np.array([.03,-.02]); g=scalar_gradient(T)
    def path_p(eps):
        mm,_=vacuum.canonical_mass_data(X0+eps*dx)
        return thermal(T*np.exp(-eps*dlogN), masses=mm)['p']
    ward_fd=(path_p(step)-path_p(-step))/(2*step)
    ward=-(a['rho']+a['p'])*dlogN-np.dot(g,dx)
    assert abs(ward_fd-ward)<1e-10
    return dict(T=T, original_canonical_radius=float(Q0), original_radial_b=float(b),
                shared_q_sigma_hessian=expected.tolist(), contact_matrix=contacts.tolist(),
                finite_Weyl_identity_error=float(finite_error), hessian_difference_error=fd_error,
                lapse_response=float(lapse_fd), volume_response=float(space_fd),
                Weyl_response=float(weyl_fd), wrong_fixed_local_temperature_response=4*a['f'],
                fixed_temperature_Weyl_error=abs(a['theta']-4*a['f']),
                hydrostatic_Ward_error=float(abs(ward_fd-ward)),
                slowly_varying_static_background_zero_derivative_only=True)


def shared_background_check():
    rows=[]
    for T in (.08,.2,.5):
        a=thermal(T); grad=scalar_gradient(T)
        euler=2*np.dot(X0,grad); expected=M/F0*a['theta']
        step=2e-5; fd=[]
        for direction in np.eye(2)*step:
            mp,_=vacuum.canonical_mass_data(X0+direction)
            mn,_=vacuum.canonical_mass_data(X0-direction)
            fd.append((thermal(T,masses=mp)['f']-thermal(T,masses=mn)['f'])/(2*step))
        assert abs(euler-expected)<2e-16
        grad_error=float(np.max(abs(np.array(fd)-grad)))
        assert grad_error<1e-10 and euler>0
        V=-a['rho']
        residual=a['p']-V
        assert abs(residual-(a['rho']+a['p']))<1e-15 and residual>0
        rows.append(dict(T=T, thermal_x_force=grad.tolist(), radial_Euler_force=float(euler),
                         trace_link_error=float(abs(euler-expected)),
                         scalar_gradient_difference_error=grad_error,
                         vacuum_adjustment_cancelling_rho=V,
                         remaining_pressure=residual,
                         original_vacuum_no_longer_stationary_in_this_branch=True))
    return dict(rows=rows, vacuum_matching_from632_held_fixed=True,
                scalar_minimum_shift_formula_not_claimed_new=True,
                obstruction_is_flat_static_constant_field_vacuum_compensation_only=True,
                no_obstruction_to_curved_or_time_dependent_GR=True)


def run():
    deps=('research_note_602.md','research_note_632.md','research_note_635.md',
          'research_note_644.md','joint_background_contact_matching.py',
          'joint_fermion_gauss_completion.py')
    return dict(round=645, tests_run=3, failures=0, errors=0,
                same_spectrum_thermal_stress=thermodynamics_check(),
                shared_thermal_geometry_sources=hydrostatic_contact_check(),
                common_background_obstruction=shared_background_check(),
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope=dict(original_full_generation_free_continuum_fermion_branch=True,
                           specified_3plus1_fixed_background_and_thermal_state_are_inputs=True,
                           same_zero_temperature_renormalization_thermal_difference_UV_finite=True,
                           no_gauge_or_scalar_thermal_loops_or_resummation=True,
                           no_full_graph_Gauss_thermal_or_continuum_mapping_claim=True,
                           no_full_area_law_dynamic_geometry_GR_or_SM_derivation_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
