"""623: original confining Hamiltonian graph-norm identities.

Numerics check target geometry, integration by parts, and subordinate terms.
The common-domain and thermal weak second-variation results are analytic;
these factor checks do not compute the full Gauss spectrum or continuum.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
import joint_fermion_gauss_completion as fermion

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_operator_domain_completion_results.json'
L, VAC, _ = original.lattice.scalar.parameters()
M = original.M
P = np.eye(2) - VAC[:, None] * np.ones((1, 2)) / (6*M)
HESS = .5 * P.T @ L @ P
LINEAR = -.5 * P.T @ L @ (VAC/M)


def polynomial(a, b):
    ab = np.stack(np.broadcast_arrays(a, b), axis=-1)
    delta = ab @ P.T - VAC/M
    value = np.einsum('...i,ij,...j->...', delta, L, delta)/4
    grad = ab @ HESS.T + LINEAR
    lap = ((8+2*a)*grad[..., 0] + (2+2*b)*grad[..., 1]
           + (4*a+2*a*a/3)*HESS[0, 0]
           + (4*b+2*b*b/3)*HESS[1, 1] + 4*a*b/3*HESS[0, 1])
    return value, grad, lap


def value_x(x):
    return polynomial(np.sum(x[..., :4]**2, axis=-1), x[..., 4]**2)[0]


def ball(x):
    return x * np.sqrt(M/(1+np.sum(x*x, axis=-1)/6))[..., None]


def constants():
    lam = np.linalg.det(L)/np.trace(L)
    eta = 1-np.sum(VAC)/(6*M)
    b0 = np.sum(VAC)/M
    assert lam > 0 and eta > 0 and b0 > 0
    threshold = 2*b0/eta
    kappa = lam*eta**2/32
    ell = float(np.max(abs(LINEAR)))
    hess = float(np.max(abs(HESS)))
    cpoly = 11*ell + 29*hess/3
    control = max(1+threshold**2, 1., 1/kappa)
    return dict(lambda_lower=float(lam), eta=float(eta), b0=float(b0),
                t_threshold=float(threshold), kappa=float(kappa),
                polynomial_constant=cpoly, coercive_control=float(control),
                laplacian_control=float(cpoly*control))


def target_laplacian_check():
    rng = np.random.default_rng(623)
    points = []
    for radius in (0., .5, 3., 10., 30.):
        directions = rng.normal(size=(5, 5))
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        points.extend(radius*directions)
    points = np.array(points)
    a = np.sum(points[:, :4]**2, axis=1)
    b = points[:, 4]**2
    values, _, exact = polynomial(a, b)
    original_values = original.node_potential(ball(points))
    match = float(np.max(abs(values-original_values)/(1+abs(values))))
    c = constants()
    bound_ratio = float(np.max(abs(exact)/(c['laplacian_control']*(1+values))))
    steps = []
    for step in (.01, .005, .0025):
        approximations = []
        for x in points:
            z = np.sqrt(1+x@x/6)
            frame = np.eye(5)+np.outer(x, x)/(6*(z+1))
            xp = (np.cosh(step/np.sqrt(6))*x[:, None]
                  + np.sqrt(6)*np.sinh(step/np.sqrt(6))*frame)
            xm = (np.cosh(step/np.sqrt(6))*x[:, None]
                  - np.sqrt(6)*np.sinh(step/np.sqrt(6))*frame)
            approximations.append(np.sum(value_x(xp.T)+value_x(xm.T)-2*value_x(x))/step**2)
        error = float(np.max(abs(np.array(approximations)-exact)/(1+abs(exact))))
        steps.append(dict(step=step, relative_error=error))
    assert match < 2e-12 and bound_ratio <= 1
    assert steps[-1]['relative_error'] < 2e-5
    assert steps[0]['relative_error'] > 3*steps[1]['relative_error']
    assert steps[1]['relative_error'] > 3*steps[2]['relative_error']
    return dict(points=len(points), original_potential_relative_error=match,
                sampled_bound_ratio=bound_ratio, constants=c, geodesic_steps=steps,
                bound_proof_is_global_polynomial_not_sampling=True)


def gauss_interval(n, low, high):
    x, w = np.polynomial.legendre.leggauss(n)
    return (high+low)/2 + (high-low)*x/2, (high-low)*w/2


def integral_identity(n):
    radius = 3.7
    rho, wr = gauss_interval(n, 0., radius)
    y, wy = gauss_interval(16, -1., 1.)
    wy = wy*.75*(1-y*y)  # normalized S^4 angular density of n_5
    r = np.sqrt(6)*np.sinh(rho/np.sqrt(6))
    aa = r[:, None]**2*(1-y*y)
    bb = r[:, None]**2*y*y
    potential, _, lap = polynomial(aa, bb)
    w = 1.3
    alpha = 1/(2*w)
    exponent = -.45+.27j
    psi = np.exp(exponent*rho*rho)
    dpsi = 2*exponent*rho*psi
    ddpsi = (2*exponent+4*exponent**2*rho*rho)*psi
    kinetic = -alpha*(ddpsi + 4/(np.sqrt(6)*np.tanh(rho/np.sqrt(6)))*dpsi)
    W = w*potential
    measure = (wr*r**4)[:, None]*wy[None, :]
    norm = float(np.sum(measure*abs(psi[:, None])**2))
    kin_norm = float(np.sum(measure*abs(kinetic[:, None])**2))
    pot_norm = float(np.sum(measure*W*W*abs(psi[:, None])**2))
    total_norm = float(np.sum(measure*abs(kinetic[:, None]+W*psi[:, None])**2))
    positive_cross = float(2*alpha*np.sum(measure*W*abs(dpsi[:, None])**2))
    lap_term = float(-alpha*w*np.sum(measure*lap*abs(psi[:, None])**2))
    rR = np.sqrt(6)*np.sinh(radius/np.sqrt(6))
    rpR = np.cosh(radius/np.sqrt(6))
    uR, gradR, _ = polynomial(rR*rR*(1-y*y), rR*rR*y*y)
    Wprime = w*2*rR*rpR*(gradR[:, 0]*(1-y*y)+gradR[:, 1]*y*y)
    psiR = np.exp(exponent*radius**2)
    derivativeR = 2*exponent*radius*psiR
    boundary = float(alpha*rR**4*np.sum(wy*(Wprime*abs(psiR)**2
                         - 2*w*uR*np.real(psiR.conjugate()*derivativeR))))
    rhs = kin_norm+pot_norm+positive_cross+lap_term+boundary
    return dict(radial_points=n, radius=radius, norm=norm,
                h_norm_squared=total_norm, t_norm_squared=kin_norm,
                w_norm_squared=pot_norm, positive_cross=positive_cross,
                laplacian_term=lap_term, finite_radius_boundary=boundary,
                relative_residual=float(abs(total_norm-rhs)/(1+abs(total_norm))),
                residual_if_boundary_omitted=float(abs(total_norm-(rhs-boundary))))


def graph_identity_check():
    rows = [integral_identity(n) for n in (48, 96)]
    assert max(r['relative_residual'] for r in rows) < 2e-12
    assert min(r['residual_if_boundary_omitted'] for r in rows) > .001
    convergence = abs(rows[0]['h_norm_squared']-rows[1]['h_norm_squared'])/(1+rows[1]['h_norm_squared'])
    assert convergence < 2e-12
    return dict(rows=rows, quadrature_relative_difference=float(convergence),
                full_original_angular_potential_retained=True,
                target_factor_identity_not_full_graph_spectral_test=True)


def subordinate_terms_check():
    direction = np.array([1., .5, -.2, .4, .7])
    direction /= np.linalg.norm(direction)
    rows = []
    for rho in (2., 4., 6., 8., 10.):
        x = np.sqrt(6)*np.sinh(rho/np.sqrt(6))*direction
        phis = (ball(x), ball(-x))
        W = float(sum(original.node_potential(p) for p in phis))
        edge = float(original.distance_squared(*phis))
        matrix_bound = 0.
        for phi in phis:
            h, d = fermion.mass_matrices(phi)
            matrix_bound += float(np.sum(abs(h))+np.sum(abs(d)))
        rows.append(dict(rho=rho, onsite=W, edge_distance_squared=edge,
                         geodesic_error=float(abs(edge-4*rho*rho)/(1+4*rho*rho)),
                         full_CAR_mass_upper_bound=matrix_bound,
                         edge_to_onsite=edge/W, mass_to_onsite=matrix_bound/W))
    assert max(r['geodesic_error'] for r in rows) < 5e-12
    for field in ('edge_to_onsite', 'mass_to_onsite'):
        assert all(rows[i+1][field] < rows[i][field] for i in range(len(rows)-1))
        assert rows[-1][field] < rows[0][field]*.002
    return dict(rows=rows, original_32_mode_complex_Dirac_Majorana_retained=True,
                numeric_branch='opposite two-node rays, identity gauge link',
                full_graph_and_all_links_bound_proved_analytically=True)


def run():
    return dict(round=623, tests_run=3, failures=0, errors=0,
                target_laplacian=target_laplacian_check(),
                graph_norm_identity=graph_identity_check(),
                subordinate_terms=subordinate_terms_check(),
                dependency_hashes={n: hashlib.sha256((HERE/n).read_bytes()).hexdigest()
                    for n in ('joint_curved_quantum_source.py','joint_fermion_gauss_completion.py',
                              'research_note_589.md','research_note_598.md','research_note_603.md',
                              'research_note_622.md')},
                scope='Original full fixed finite graph at positive external geometry; fixed target and onsite parameters. Common operator domain and thermal weak second variations proved in note. No continuum, dynamic gravity, source essential self-adjointness, or strong second derivative claimed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf8', newline='\n') as f:
            f.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
