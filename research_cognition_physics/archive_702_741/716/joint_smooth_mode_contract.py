"""716: original sterile CAR readout, physical volume and geometric profile.

Pointwise full 32-mode coefficients test analytic identities; neither samples
nor spatial quadratures are a simulation of the interacting Gibbs state.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import joint_vertex_shared_evolution as car
import joint_region_energy_gluing as group

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE/'round716_drafts'))
import sterile_mode_entry as entry
TARGET = HERE/'joint_smooth_mode_contract_results.json'
matter = car.matter
geom = matter.original
L, vacuum, _ = geom.lattice.scalar.parameters()
lam = float(np.linalg.det(L)/np.trace(L))
D = float(6*geom.M-vacuum.sum())
A = 32/(lam*D*D)
B0 = 144/(D*D)
C = float(np.sqrt(6*geom.M)*(abs(matter.Y['nu'])+abs(matter.Y['s'])))


def profile(w, spinor):
    z = float(np.sum(w*np.sum(abs(spinor)**2, axis=1)))
    local = np.sqrt(w[:, None])*spinor/np.sqrt(z)
    f = np.zeros((len(w), 32), complex)
    f[:, 30:32] = local
    return f.ravel(), z


def coefficients(points, links, strength):
    n = 32*len(points)
    mass = np.zeros((n, n), complex)
    delta = np.zeros_like(mass)
    hop = np.zeros_like(mass)
    for v, point in enumerate(points):
        s = slice(32*v, 32*v+32)
        mass[s, s], delta[s, s] = matter.mass_matrices(point)
    for i, (g, t) in enumerate(zip(links, strength)):
        j = (i+1) % len(points)
        a, b = slice(32*i, 32*i+32), slice(32*j, 32*j+32)
        block = t*matter.representation(*g)
        hop[a, b] += block
        hop[b, a] += block.conj().T
    return mass, delta, hop


def difference_coeff(h, d, f):
    P = np.outer(f, f.conj())
    U = np.eye(len(f))-2*P
    return (U@h@U-h)/2, (U@d@U.T-d)/2


def exact_norm(dh, dd, span):
    v, s, _ = np.linalg.svd(span, full_matrices=False)
    v = v[:, s > 1e-11]
    hs, ds = v.conj().T@dh@v, v.conj().T@dd@v.conj()
    residual = max(np.linalg.norm(dh-v@hs@v.conj().T),
                   np.linalg.norm(dd-v@ds@v.T))
    n = v.shape[1]
    B = np.zeros((2**n, 2**n), complex)
    for col in range(2**n):
        for row, value in car.quadratic({col: 1.}, hs, ds).items():
            B[row, col] = value
    return float(np.max(abs(np.linalg.eigvalsh(B)))), float(residual), n


def mass_reference_check():
    rng = np.random.default_rng(71621)
    rows = []
    max_error = 0.
    for k in range(15):
        dirs = rng.normal(size=(3, 5))
        dirs /= np.linalg.norm(dirs, axis=1)[:, None]
        fs = np.exp(rng.uniform(np.log(2e-4), np.log(1.9), size=3))
        points = dirs*np.sqrt(6*(geom.M-fs))[:, None]
        w = np.exp(rng.uniform(-4, 1, size=3))
        spinor = rng.normal(size=(3, 2))+1j*rng.normal(size=(3, 2))
        f, z = profile(w, spinor)
        p = np.sum(abs(f.reshape(3, 32))**2, axis=1)
        links = [group.sample(rng) for _ in range(3)]
        t = .23*np.exp(rng.uniform(-.5, .5, size=3))
        mass, d, hop = coefficients(points, links, t)
        h = mass+hop
        P = np.outer(f, f.conj())
        Q = np.eye(len(f))-P
        hf, bf = Q@h@f, Q@d@f.conj()
        dh, dd = difference_coeff(h, d, f)
        norm, residual, dim = exact_norm(dh, dd, np.column_stack((f, hf, bf)))
        F = geom.F(points)
        U = geom.node_potential(points)
        md2 = abs(matter.Y['nu'])**2*np.sum(p*np.sum(points[:, :4]**2, axis=1)/F)
        mm2 = abs(matter.Y['s'])**2*np.sum(p*points[:, 4]**2/F)
        errors = [abs(np.linalg.norm(mass@f)**2-md2)/(1+md2),
                  abs(np.linalg.norm(d@f.conj())**2-mm2)/(1+mm2),
                  abs(np.vdot(f, d@f.conj())), residual/(1+norm)]
        max_error = max(max_error, *errors)
        kappa = float(np.linalg.norm(Q@hop@f))
        Pf = float(p@U)
        r = float(np.max(p/w))
        W = float(w@U)
        local = kappa+C*(A*Pf+B0)**.25
        volume = kappa+C*(A*r*W+B0)**.25
        coefficient_bound = float(np.linalg.norm(hf)+np.linalg.norm(bf))
        assert norm <= coefficient_bound+1e-10 <= local+1e-9 <= volume+1e-8
        assert dim <= 3 and Pf <= r*W+1e-8
        rows.append(dict(sample=k,full_original_modes=len(f),F_min=float(min(F)),
                         exact_CAR_difference_norm=norm,coefficient_bound=coefficient_bound,
                         profile_potential=Pf,local_potential_bound=local,
                         volume_bound=volume,kappa=kappa,profile_density_max=r))
    assert max_error < 1e-10
    return dict(pointwise_samples=15,original_modes=96,
                maximum_relative_identity_error=max_error,rows=rows,
                full_Gibbs_not_numerically_sampled=True)


def volume_check():
    rows = []
    for dimension in (2, 3):
        for n in (12, 20, 32):
            a = 4/n
            axes = [(np.arange(n)+.5)*a-2 for _ in range(dimension)]
            x = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, dimension)
            radius2 = np.sum(x*x, axis=1)
            bump = np.zeros(len(x))
            inside = radius2 < 1
            bump[inside] = np.exp(-1/(1-radius2[inside]))
            spinor = np.column_stack((bump, .2j*bump))
            density = np.exp(.18*np.cos(x[:, 0])-.12*np.sin(x[:, 1]))
            w = a**dimension*density
            z = float(np.sum(w*np.sum(abs(spinor)**2, axis=1)))
            p = w*np.sum(abs(spinor)**2, axis=1)/z
            r = float(np.max(p/w))
            points = np.column_stack((.3+.06*np.cos(x[:, 0]),
                                       .2*np.sin(x[:, 1]), .12*np.ones(len(x)),
                                       .08*np.cos(x[:, 0]+x[:, 1]),
                                       .35+.05*np.sin(x[:, 0])))
            U = geom.node_potential(points)
            F = geom.F(points)
            Pf = float(p@U)
            d2 = abs(matter.Y['nu'])**2*np.sum(p*np.sum(points[:, :4]**2, axis=1)/F)
            m2 = abs(matter.Y['s'])**2*np.sum(p*points[:, 4]**2/F)
            bound = C*(A*Pf+B0)**.25
            assert abs(p.sum()-1) < 1e-13 and np.sqrt(d2)+np.sqrt(m2) <= bound
            assert abs(r-np.max(np.sum(abs(spinor)**2, axis=1))/z) < 1e-12
            rows.append(dict(dimension=dimension,n=n,nodes=len(x),physical_step=a,
                             normalization=z,profile_density_max=r,
                             profile_potential=Pf,mass_action_sum=float(np.sqrt(d2)+np.sqrt(m2)),
                             mass_energy_bound=bound))
    return dict(rows=rows,quadratures_are_not_a_quantum_reference_state=True,
                same_normalization_works_in_both_dimensions=True,
                no_spatial_dimension_or_continuum_dynamics_claim=True)


def geometry_check():
    rng = np.random.default_rng(71631)
    points = np.array([[.4,-.3,.2,.1,.35],[-.2,.15,.1,-.25,.3],
                       [.22,.34,-.12,.08,-.27]])
    links = [group.sample(rng) for _ in range(3)]
    spinor = rng.normal(size=(3, 2))+1j*rng.normal(size=(3, 2))
    w0 = np.array([.7, 1.1, .5])
    slope = np.array([.8,-.4, .25])
    tau = np.array([-.21,.13,-.31])
    t0 = .23*np.exp(np.array([-.24,.14,.05]))
    mass, d, hop = coefficients(points, links, t0)
    _, _, hopdot = coefficients(points, links, t0*tau)
    h = mass+hop
    f, _ = profile(w0, spinor)
    p = np.sum(abs(f.reshape(3, 32))**2, axis=1)
    centered = slope-p@slope
    fdot = (f.reshape(3, 32)*(.5*centered[:, None])).ravel()
    P = np.outer(f, f.conj())
    Pdot = np.outer(fdot, f.conj())+np.outer(f, fdot.conj())
    U, Udot = np.eye(len(f))-2*P, -2*Pdot
    fixed_h, fixed_d = difference_coeff(hopdot, np.zeros_like(d), f)
    instrument_h = (Udot@h@U+U@h@Udot)/2
    instrument_d = (Udot@d@U.T+U@d@Udot.T)/2
    exact_h = fixed_h+instrument_h
    exact_d = fixed_d+instrument_d
    def at(s):
        fs, _ = profile(w0*np.exp(s*slope), spinor)
        _, _, hops = coefficients(points, links, t0*np.exp(s*tau))
        return fs, difference_coeff(mass+hops, d, fs)
    derivative_rows = []
    for e in (2e-3, 7e-4, 2e-4):
        fp, (hp, dp) = at(e)
        fm, (hm, dm) = at(-e)
        error = max(np.linalg.norm((hp-hm)/(2*e)-exact_h),
                    np.linalg.norm((dp-dm)/(2*e)-exact_d),
                    np.linalg.norm((fp-fm)/(2*e)-fdot))
        omission = max(np.linalg.norm((hp-hm)/(2*e)-fixed_h),
                       np.linalg.norm((dp-dm)/(2*e)-fixed_d))
        derivative_rows.append(dict(step=e,error=float(error),omit_profile_error=float(omission)))
    variance_error = float(abs(np.linalg.norm(fdot)**2-.25*np.sum(p*centered**2)))
    uniform_f, _ = profile(w0*np.exp(.7), spinor)
    Q = np.eye(len(f))-P
    ad = -Pdot@hop@f+Q@hopdot@f+Q@hop@fdot
    a = Q@h@f
    b = d@f.conj()
    aprime = ad+mass@fdot
    bprime = d@fdot.conj()
    derivative_norm, residual, dim = exact_norm(
        exact_h, exact_d, np.column_stack((f,fdot,a,b,aprime,bprime)))
    S = float(np.max(abs(centered)))
    Pf = float(p@geom.node_potential(points))
    kappa = float(np.linalg.norm(Q@hop@f))
    kappa_derivative = float(np.linalg.norm(ad))
    bound = S*kappa+2*kappa_derivative+2*S*C*(A*Pf+B0)**.25
    assert derivative_rows[-1]['error'] < 1e-8
    assert derivative_rows[-1]['omit_profile_error'] > .01
    assert variance_error < 1e-14 and np.linalg.norm(uniform_f-f) < 1e-14
    assert derivative_norm <= bound and residual < 1e-10
    assert abs(np.vdot(f,fdot)) < 1e-14
    return dict(full_original_modes=len(f),finite_difference_rows=derivative_rows,
                metric_variance_identity_error=variance_error,
                profile_derivative_norm=float(np.linalg.norm(fdot)),
                exact_source_difference_derivative_norm=derivative_norm,
                source_derivative_bound=bound,active_derivative_modes=dim,
                source_support_error=residual,
                uniform_volume_rescaling_profile_error=float(np.linalg.norm(uniform_f-f)),
                no_thermal_covariance_or_autonomous_device_computed=True)


def run():
    result = dict(round=716,tests_run=3,failures=0,errors=0,
                  constants=dict(A=A,B0=B0,C=C),
                  mass_reference=mass_reference_check(),physical_volume=volume_check(),
                  geometry_profile=geometry_check())
    names = ('research_note_598.md','research_note_623.md','research_note_633.md',
             'research_note_667.md','research_note_714.md',
             'joint_fermion_gauss_completion.py','joint_vertex_shared_evolution.py',
             'round716_drafts/sterile_mode_entry.py')
    result['dependency_hashes'] = {name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                  for name in names}
    result['scope'] = dict(original_full_finite_graph=True,
                           uniform_estimates_conditional_on_profile_and_reference=True,
                           dimension_not_selected=True,spatial_continuum_not_proved=True,
                           autonomous_causal_instrument_not_proved=True,
                           full_interacting_Gibbs_not_numerically_computed=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
