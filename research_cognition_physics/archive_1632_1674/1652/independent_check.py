"""Cross-check identifiability using eigenvectors and a coordinate-fixed square inverse."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(relative):
    return json.loads((HERE / relative).read_bytes())


def evaluate():
    inp = read('inputs.json')
    old = read(inp['texture'])['reconstruction']
    edges = old['normalized_real_edges_A_B_C_E_F'] + list(old['two_real_coefficients'])
    edges[-1] *= np.exp(1j*np.deg2rad(old['eta_example_deg']))
    cells = [(0, 1), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2), (1, 2)]
    basis = np.zeros((7, 3, 3), complex)
    for n, (i, j) in enumerate(cells):
        basis[n, i, j] = edges[n]
    mass = basis.sum(axis=0)
    square_mass, left = np.linalg.eigh(mass @ mass.conj().T)
    right = mass.conj().T @ left / np.sqrt(square_mass)
    np.testing.assert_allclose(right.conj().T @ right, np.eye(3), atol=2e-11)
    matrices = left.conj().T @ basis @ right

    def coordinates(m):
        a, b = m[0, 2] - m[2, 0].conjugate(), m[1, 2] - m[2, 1].conjugate()
        c = m[1, 2] + m[2, 1].conjugate()
        return [a.real, a.imag, b.real, b.imag, c.real, c.imag]

    T = np.column_stack([coordinates(g) for g in matrices])
    np.testing.assert_allclose(T.sum(axis=1), 0, atol=1e-11)
    contrasts = np.array([[0, 0, 0, 0, 0, -1, 1], [0, 0, -1, 1, -1, 0, 1]])
    free = T[:, 1:]  # Gauge A's slope to zero.
    # All seven oriented column minors have the same magnitude because the
    # kernel is the common gain. Cauchy–Binet fixes their normalization.
    singular_values = np.linalg.svd(T, compute_uv=False)
    assert singular_values[-1] > 1.6e-4
    np.testing.assert_allclose(abs(np.linalg.det(free)), np.prod(singular_values)/np.sqrt(7),
                               rtol=1e-8, atol=0)
    inverse = np.linalg.solve(free.T, contrasts[:, 1:].T).T
    np.testing.assert_allclose(inverse @ T, contrasts, atol=1e-10)
    sensitivity = np.linalg.norm(inverse.reshape(2, 3, 2), axis=2)
    assert min(sensitivity[:, 2]) > 0
    for i in range(2):
        unit_phases = inverse[i].reshape(3, 2) / sensitivity[i, :, None]
        for radii in ([1, 0, 0], [.2, .6, .1], [0, 0, 1]):
            test = (np.array(radii)[:, None]*unit_phases).ravel()
            realization = np.r_[0, np.linalg.solve(free, test)]
            np.testing.assert_allclose(contrasts[i] @ realization, sensitivity[i] @ radii,
                                       atol=5e-10)
    rowcol = np.zeros((7, 6))
    for n, (i, j) in enumerate(cells):
        rowcol[n, [i, j+3]] = 1
    assert np.linalg.matrix_rank(T[:4] @ rowcol, tol=1e-8) == 4

    # Independent lower-envelope quadrature at coefficient corners. The
    # interval proof in the main report bounds the whole bin below this value.
    kin = read(inp['semileptonic'])
    params, masses = kin['HPQCD_table_V'], kin['kinematics_GeV']
    q0, q1 = kin['q2_bin_GeV2']
    nodes, weights = np.polynomial.legendre.leggauss(32)
    q = (q0+q1)/2 + (q1-q0)*nodes/2
    mB, mK, mu = (masses[key] for key in ['Bplus','Kplus','mu'])
    threshold = (params['B_average_GeV']+params['K_average_GeV'])**2
    z = (np.sqrt(threshold-q)-math.sqrt(threshold))/(np.sqrt(threshold-q)+math.sqrt(threshold))
    def band(key):
        val, err = params[key]
        return [val-2*err, val+2*err]

    # Rate for the dimensionless parity readout before the independent
    # kappa_mu, 1-rs, m_mu/(2v^2) and lifetime factors.
    box_rates = []
    for a, b, d, L, pole in itertools.product(*(band(key) for key in ['a0','a1','a2','L','pole_mass_GeV'])):
        form = L*(a + b*z + d*z*z)/(1 - q/pole**2)
        lam = ((mB+mK)**2-q)*((mB-mK)**2-q)
        kernel = np.sqrt(lam)*q*(1-4*mu**2/q)**1.5*form**2/(inp['mh_GeV']**2-q)**2
        box_rates.append((q1-q0)/2 * (weights@kernel)*(mB**2-mK**2)**2/(128*math.pi**3*mB**3))
    main = read('results.json')
    np.testing.assert_allclose(main['dual_group_norms'], sensitivity, rtol=1e-8)
    mass_ratio = np.sqrt(square_mass[1] / square_mass[2])
    corrected_min = min(box_rates) / (1 - mass_ratio)**2
    interval_floor = main['form_factor_integral_with_propagator_and_mass_ratio_GeV3']
    assert .995*corrected_min < interval_floor <= corrected_min
    radii = np.sqrt(np.array(main['experimental_rate_upper']) /
                    np.array(main['rate_floors_per_unit_kmu_times_projection_squared']))
    np.testing.assert_allclose(sensitivity @ radii,
                               main['limits']['rate_only']['support_abs_kmu_contrasts_nominal'],
                               rtol=1e-8)
    return {'absolute_det_with_A_fixed': abs(np.linalg.det(free)),
            'sixth_singular_value': np.linalg.svd(T, compute_uv=False)[-1],
            'pair_modulus_sensitivities_in_order_Bd_Bs_K': sensitivity.tolist(),
            'no_plus_readout_leaves_free_dimension_at_fixed_A': 2,
            'K_bin_corner_kernel_minmax_no_quark_mass_or_vev_factors': [min(box_rates), max(box_rates)]}


def tidy(d):
    if isinstance(d, dict):
        return {k: tidy(v) for k, v in d.items()}
    if isinstance(d, list):
        return [tidy(v) for v in d]
    if isinstance(d, (float, np.floating)):
        return float(f'{d:.9g}')
    return d


if __name__ == '__main__':
    cli = argparse.ArgumentParser()
    cli.add_argument('--check', action='store_true')
    args = cli.parse_args()
    result = tidy(evaluate())
    target = HERE / 'independent_results.json'
    if args.check:
        assert json.loads(target.read_bytes()) == result
        print('PASS: alternate mass basis, square inverse, extremal phases, hadronic kernel')
    else:
        with target.open('x', encoding='utf-8') as stream:
            stream.write(json.dumps(result, indent=2) + '\n')
