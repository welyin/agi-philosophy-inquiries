"""731: simultaneous initial-constraint response on the original 572 background.

Declared smooth source tests are not a computed quantum state. The analytic
right inverse and first-order application are stated separately in note731.
No Einstein/SM emergence or finite-strength semiclassical solution is claimed.
"""
import argparse
import hashlib
import json
from functools import lru_cache
from pathlib import Path
import numpy as np
import joint_gauss_einstein_initial_data as geo

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_source_constraint_response_results.json'


def maximum(v):
    return float(np.max(np.abs(v)))


def tangent_matrix(q):
    X = q['f']['X']
    TX = 1j * np.einsum('aij,...j->...ia', geo.old.T, X)
    TX = np.concatenate((TX, (3j * X)[..., None]), axis=-1)
    return np.concatenate((TX.real, TX.imag, np.zeros(X.shape[:-1] + (1, 4))), axis=-2)


def covariant(v, i, q):
    return geo.derivative(v, i) - np.cross(q['f']['a'][..., i, :], v)


def gauss(q, k, p, E, E0):
    out = np.einsum('...ag,...a->...g', k, p)
    out[..., :3] += sum(covariant(E[..., i, :], i, q) for i in range(3))
    out[..., 3] += sum(geo.derivative(E0[..., i], i) for i in range(3))
    return out


def curvature(q):
    a, a0 = q['f']['a'], q['f']['a0']
    fw = np.zeros(a.shape[:-2] + (3, 3, 3))
    f0 = np.zeros(a0.shape[:-1] + (3, 3))
    for i in range(3):
        for j in range(3):
            fw[..., i, j, :] = (geo.derivative(a[..., j, :], i)
                                - geo.derivative(a[..., i, :], j)
                                - np.cross(a[..., i, :], a[..., j, :]))
            f0[..., i, j] = geo.derivative(a0[..., j], i) - geo.derivative(a0[..., i], j)
    return fw, f0


def momentum(q, p, E, E0):
    fw, f0 = curvature(q)
    return (np.einsum('...a,...ia->...i', p, q['Dphi'])
            + np.einsum('...ja,...ija->...i', E, fw)
            + np.einsum('...j,...ij->...i', E0, f0))


def kinverse(q, p):
    phi = q['phi']
    return q['F'][..., None] * (p - phi * np.sum(phi * p, axis=-1)[..., None] / 12)


def inverse_gauss(q, k, sigma):
    gram = np.einsum('...ag,...ah->...gh', k, k)

    def operator(v):
        out = np.einsum('...gh,...h->...g', gram, v)
        for i in range(3):
            out[..., :3] -= covariant(covariant(v[..., :3], i, q), i, q)
            out[..., 3] -= geo.derivative(geo.derivative(v[..., 3], i), i)
        return out

    # Odd N avoids the real-collocation Nyquist derivative kernel.
    assert q['N'] % 2 == 1
    potential = gram.mean(axis=(0, 1, 2))
    for i in range(3):
        a = q['f']['a'][..., i, :]
        potential[:3, :3] += (np.sum(a*a, axis=-1)[..., None, None] * np.eye(3)
                              - a[..., :, None] * a[..., None, :]).mean(axis=(0, 1, 2))
    eigen, frame = np.linalg.eigh(potential)
    assert eigen.min() > 0
    k2 = np.sum(geo.waves(q['N'])**2, axis=-1)

    def precondition(v):
        vh = np.einsum('...g,gh->...h', geo.fft(v), frame)
        return geo.ifft(np.einsum('...h,gh->...g', vh / (k2[..., None] + eigen), frame))

    xi, iterations = geo.cg(operator, sigma, precondition, tol=3e-13)
    dp = -np.einsum('...ag,...g->...a', k, xi)
    dE = np.stack([covariant(xi[..., :3], i, q) for i in range(3)], axis=-2)
    dE0 = np.stack([geo.derivative(xi[..., 3], i) for i in range(3)], axis=-1)
    return dp, dE, dE0, dict(iterations=iterations,
        residual=maximum(operator(xi)-sigma), preconditioner_min_eigenvalue=float(eigen.min()))


@lru_cache(maxsize=1)
def setup():
    q = geo.make_source(15)
    psi, tensor, base_info = geo.solve_hamiltonian(q)
    k = tangent_matrix(q)
    assert maximum(gauss(q, k, q['p'], q['f']['E'], q['f']['E0'])) < 2e-13
    assert maximum(momentum(q, q['p'], q['f']['E'], q['f']['E0'])-q['mom']) < 2e-13
    x, y, z = np.moveaxis(q['grid'], -1, 0)
    sigma = np.stack((.012*np.cos(x+y), .008*np.sin(y+z),
                      .009*np.cos(z)-.005, .006+.004*np.sin(x)), axis=-1)
    J = np.stack((.003+.004*np.sin(x), -.002+.003*np.cos(y),
                  .001+.002*np.sin(z)), axis=-1)
    rho = .01+.002*np.cos(x-y)+.003*np.sin(z)
    dp, dE, dE0, solve_info = inverse_gauss(q, k, sigma)
    particular_mom = momentum(q, dp, dE, dE0)
    dh = np.stack([geo.derivative(q['f']['h'], i) for i in range(3)], axis=-1)
    ds = np.stack([geo.derivative(q['f']['s'], i) for i in range(3)], axis=-1)
    menus = []
    for i in (0, 1):
        pm = np.zeros_like(dp); pm[..., 1] = dh[..., i]; pm[..., 4] = ds[..., i]
        menus.append((pm, np.zeros_like(dE), np.zeros_like(dE0)))
    fw, _ = curvature(q)
    weighted = np.sin(2*z)[..., None] * fw[..., 0, 2, :]
    electric = np.zeros_like(dE)
    electric[..., 0, :] = covariant(weighted, 2, q)
    electric[..., 2, :] = -covariant(weighted, 0, q)
    menus.append((np.zeros_like(dp), electric, np.zeros_like(dE0)))
    columns = np.column_stack([q['dx']**3 * momentum(q, *m).sum(axis=(0, 1, 2)) for m in menus])
    before = q['dx']**3 * (particular_mom+J).sum(axis=(0, 1, 2))
    coeff = np.linalg.solve(columns, -before)
    for c, (pm, em, e0m) in zip(coeff, menus):
        dp += c*pm; dE += c*em; dE0 += c*e0m
    dmom = momentum(q, dp, dE, dE0) + J
    # Separate zero-mean color source, at the original zero color configuration.
    sigma_c = np.zeros(q['grid'].shape[:-1] + (8,))
    sigma_c[..., 0] = .004*np.sin(x); sigma_c[..., 1] = .003*np.cos(y)
    k2 = np.sum(geo.waves(q['N'])**2, axis=-1); safe = k2.copy(); safe[0, 0, 0] = 1
    potential = geo.fft(sigma_c) / safe[..., None]; potential[0, 0, 0] = 0
    potential = geo.ifft(potential)
    dEc = np.stack([geo.derivative(potential, i) for i in range(3)], axis=-2)
    color_error = maximum(sum(geo.derivative(dEc[..., i, :], i) for i in range(3))+sigma_c)
    return dict(q=q, psi=psi, tensor=tensor, base_info=base_info, k=k, sigma=sigma,
        J=J, rho=rho, dp=dp, dE=dE, dE0=dE0, dEc=dEc, dmom=dmom,
        columns=columns, coeff=coeff, before=before, solve_info=solve_info,
        color_error=color_error)


def gauss_and_zero_modes_check():
    d = setup(); q = d['q']
    err = maximum(gauss(q, d['k'], d['dp'], d['dE'], d['dE0'])+d['sigma'])
    mean = q['dx']**3*d['dmom'].sum(axis=(0, 1, 2))
    assert err < 3e-12 and maximum(mean) < 3e-12 and d['color_error'] < 1e-13
    assert abs(np.linalg.det(d['columns'])) > 1e-6
    # An explicitly constant nonzero color charge is outside div(E)'s range.
    incompatible_color_integral = (2*np.pi)**3 * .01
    return dict(N=q['N'], prescribed_sources_not_quantum_state=True,
        gauss_solver=d['solve_info'], pointwise_Gauss_error=err,
        zero_mean_color_error=d['color_error'], incompatible_constant_color_integral=incompatible_color_integral,
        momentum_before=d['before'].tolist(), correction_coefficients=d['coeff'].tolist(),
        momentum_after=mean.tolist(), augmented_determinant=float(np.linalg.det(d['columns'])))


def updated_data(eps):
    d = setup(); q = d['q']; z = dict(q)
    p = q['p'] + eps*d['dp']; E = q['f']['E']+eps*d['dE']; E0 = q['f']['E0']+eps*d['dE0']
    bc, bw, b0 = geo.old.PAR['b']
    z['pKp'] = np.sum(p*kinverse(q, p), axis=-1)
    z['Y'] = (q['Y'] + bw*(np.sum(E*E, axis=(-1, -2))-np.sum(q['f']['E']**2, axis=(-1, -2)))
              + b0*(np.sum(E0*E0, axis=-1)-np.sum(q['f']['E0']**2, axis=-1))
              + bc*eps**2*np.sum(d['dEc']**2, axis=(-1, -2)))
    z['C'] = q['C'] - 2*eps*d['rho']
    z['mom'] = momentum(q, p, E, E0)+eps*d['J']
    assert maximum(z['mom']-(q['mom']+eps*d['dmom'])) < 1e-13
    assert np.min(z['C']) > 0 and np.min(z['Y']) > q['Y_lower']
    return z, p, E, E0


def independent_constraint(z, psi, tensor, eps):
    rho = (.5*psi**-12*z['pKp']+.5*psi**-4*z['B']+z['U']
           + psi**-8*z['Y']+eps*setup()['rho'])
    residual = (-8*psi**-5*geo.laplace(psi)-psi**-12*np.sum(tensor*tensor, axis=(-1,-2))
                + 2*z['tau2']/3-2*rho)
    return maximum(residual)


def nonlinear_joint_constraint_check():
    d = setup(); rows = []
    for eps in (.1, .05, -.05):
        z, p, E, E0 = updated_data(eps)
        psi, tensor, info = geo.solve_hamiltonian(z, initial=d['psi'])
        h_error = independent_constraint(z, psi, tensor, eps)
        g_error = maximum(gauss(z, d['k'], p, E, E0)+eps*d['sigma'])
        vector_error = maximum(sum(geo.derivative(tensor[..., :, i], i) for i in range(3))+z['mom'])
        assert h_error < 3e-8 and g_error < 3e-12 and vector_error < 3e-12
        assert abs(info['original_Hamiltonian_residual']-2*abs(eps)*float(d['rho'].max())) < 3e-8
        rows.append(dict(epsilon=eps, Gauss_error=g_error, momentum_error=vector_error,
            Hamiltonian_with_full_source_error=h_error, psi_min=float(psi.min()), psi_max=float(psi.max()),
            psi_change=maximum(psi-d['psi']), C_min=float(z['C'].min()),
            scalar_kinetic_change=maximum(z['pKp']-d['q']['pKp']),
            electric_change=maximum(z['Y']-d['q']['Y']),
            error_if_extra_energy_omitted=info['original_Hamiltonian_residual']))
    return dict(rows=rows, original_magnetic_scalar_configuration_and_tau_retained=True,
                full_nonlin_prescribed_source_only=True)


def linearized_response_check():
    d = setup(); q, psi, tensor = d['q'], d['psi'], d['tensor']
    delta_tensor, _ = geo.solve_momentum(dict(q, mom=d['dmom']))
    A = np.sum(tensor*tensor, axis=(-1,-2))+q['pKp']
    dA = 2*np.sum(tensor*delta_tensor, axis=(-1,-2))+2*np.sum(kinverse(q,q['p'])*d['dp'],axis=-1)
    _, bw, b0 = geo.old.PAR['b']
    dY = 2*bw*np.sum(q['f']['E']*d['dE'],axis=(-1,-2))+2*b0*np.sum(q['f']['E0']*d['dE0'],axis=-1)
    potential = 5*q['C']*psi**4-q['B']+7*A*psi**-8+6*q['Y']*psi**-4
    operator = lambda v: -8*geo.laplace(v)+potential*v
    positive = 4*q['C']*psi**5+8*A*psi**-7+8*q['Y']*psi**-3
    identity_error = maximum(operator(psi)-positive)
    rhs = dA*psi**-7+2*dY*psi**-3+2*d['rho']*psi**5
    k2 = np.sum(geo.waves(q['N'])**2,axis=-1)
    response, it = geo.cg(operator,rhs,lambda v:geo.ifft(geo.fft(v)/(8*k2+potential.mean())))
    assert identity_error < 3e-9 and positive.min()>0 and maximum(operator(response)-rhs)<1e-10
    rows=[]
    for eps in (.02,.01,.005):
        zp,*_=updated_data(eps); zm,*_=updated_data(-eps)
        pp,_,_=geo.solve_hamiltonian(zp, initial=psi); pm,_,_=geo.solve_hamiltonian(zm, initial=psi)
        error=maximum((pp-pm)/(2*eps)-response)
        remainder=maximum(pp-psi-eps*response)
        rows.append(dict(epsilon=eps, centered_derivative_error=error, one_sided_remainder=remainder,
                         remainder_over_epsilon_squared=remainder/eps**2))
    assert rows[-1]['centered_derivative_error'] < rows[0]['centered_derivative_error']/10
    assert max(r['remainder_over_epsilon_squared'] for r in rows) < 1.2*min(r['remainder_over_epsilon_squared'] for r in rows)
    return dict(positive_supersolution_identity_error=identity_error, Lpsi_positive_min=float(positive.min()),
        response_max=maximum(response), cg_iterations=it, response_equation_error=maximum(operator(response)-rhs),
        rows=rows, full_state_dependent_quantum_feedback_not_proven=True)


def run():
    checks=('gauss_and_zero_modes_check','nonlinear_joint_constraint_check','linearized_response_check')
    evidence={name:globals()[name]() for name in checks}
    deps=('joint_gauss_einstein_initial_data.py','joint_gauss_einstein_initial_data_results.json',
          'joint_dynamic_continuum_reference.py','joint_dynamic_continuum_reference_results.json',
          'round731_drafts/gauge_counterflow_entry.py','round731_drafts/gauge_counterflow_entry_results.json')
    return dict(round=731,tests_run=len(checks),failures=0,errors=0,checks=list(checks),evidence=evidence,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        scope='Original572 fixed configurations and CMC on T3. Explicit joint right inverse for prescribed smooth sources with zero total color, full energy retained, nonlinear scalar constraint and first response. Classical conformal framework is input; diagnostic sources are not a computed quantum state. No full quantum Gauss, self-consistent semiclassical dynamics, new spacetime derivation or continuous numerical certificate.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
