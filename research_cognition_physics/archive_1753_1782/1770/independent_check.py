"""Test two physical SU(5) currents on identical, threshold-matched pole spectra."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
LOW = runpy.run_path(str(HERE.parent / '1769/check.py'))
LIMITS = json.loads((HERE.parent / '1768/results.json').read_bytes())
TWO_PI = 2 * math.pi
X_INDEX = np.array([5., 3., 2.])
UNBROKEN_INDEX = np.array([0., 2., 3.])
TRIPLET = np.array([.4, 0., 1.])
REFERENCE = 1.e16


def single_group(low_doublet):
    weak = LOW['floor_spectrum'](low_doublet, LOW['NEW_WEAK'])
    ew = LOW['gauge'](weak)
    low_alpha = LOW['CFG']['strong_lattice_5f']['ALPHA_2026'][0]
    a3 = 1 / low_alpha + math.log(LOW['CFG']['top_MS_matching_GeV'][0] / LOW['MZ']) / (3 * math.pi)
    b = LOW['SUSY'] + 4
    a = ew['inverse_alpha_U'] - b * math.log(REFERENCE / ew['crossing_GeV']) / TWO_PI
    a[2] += a3 - ew['a3_SM']
    # Set the comparison single X mass equal to the matching point. Solve the
    # three linear equations for alpha5 inverse and the two chiral mass logs.
    mat = np.column_stack((np.ones(3), TRIPLET / TWO_PI, UNBROKEN_INDEX / TWO_PI))
    a5, log_h, log_sigma = np.linalg.solve(mat, a)
    np.testing.assert_allclose(mat @ [a5, log_h, log_sigma], a, atol=3e-14, rtol=0)
    return a, a5, REFERENCE * math.exp(-log_h), REFERENCE * math.exp(-log_sigma)


def two_groups(low_doublet):
    a, a5, h, sigma = single_group(low_doublet)
    g = math.sqrt(4 * math.pi / a5)
    couplings = np.array([math.sqrt(3) * g, math.sqrt(1.5) * g])
    diag = np.diag(couplings)
    mass0 = REFERENCE
    # Lengths in units of mass0 for numerical conditioning. f sets the heavy
    # unbroken-adjoint mass; v includes the broken generator normalization.
    f, v = 1 / np.linalg.norm(couplings), 1 / g
    link = f**2 * np.array([[1., -1.], [-1., 1.]])
    broken = link + np.diag([v*v, 0.])
    eigen, vec = np.linalg.eigh(diag @ broken @ diag)
    poles = np.sqrt(eigen)
    unbroken, eig_c = np.linalg.eigh(diag @ link @ diag)
    mc = math.sqrt(unbroken[-1])
    np.testing.assert_allclose(unbroken, [0., 1.], atol=1e-15)
    np.testing.assert_allclose(abs(diag @ eig_c[:, 0]), [g, g], atol=1e-15)
    np.testing.assert_allclose(poles, np.sqrt([2-math.sqrt(3), 2+math.sqrt(3)]), atol=1e-15)
    np.testing.assert_allclose(poles.prod() / mc, g*v, atol=1e-15)

    # New vector thresholds minus the reference vector. Full link radial
    # chiral multiplets are fixed at the matching scale in the example.
    vector_shift = -2 * (X_INDEX * (-np.log(poles).sum()) - UNBROKEN_INDEX * math.log(mc)) / TWO_PI
    np.testing.assert_allclose(vector_shift, 0., atol=1e-14)
    np.testing.assert_allclose(X_INDEX + UNBROKEN_INDEX, 5 * np.ones(3))
    tree = np.linalg.inv(broken)
    current_couplings = diag @ vec
    pole_parts = current_couplings**2 / eigen
    np.testing.assert_allclose(pole_parts.sum(axis=1), tree.diagonal(), rtol=2e-15)
    np.testing.assert_allclose(tree / g**2, [[1., 1.], [1., 5.5]], rtol=2e-15)
    kinetic = np.diag(1 / couplings**2)
    for q2 in (0., .01, .6, 5.):
        pure_exchange = (current_couplings / (eigen + q2)) @ current_couplings.T
        np.testing.assert_allclose(pure_exchange, np.linalg.inv(broken + q2*kinetic), rtol=2e-14)
    heavy_sm = (diag @ eig_c[:, -1])**2 / unbroken[-1]
    np.testing.assert_allclose(heavy_sm[0] / heavy_sm[1], 4., atol=3e-15)
    bare = mass0 / np.sqrt(tree.diagonal())
    lowest = min(row['lower_upper_bracket_on_required_scale_GeV'][0] for row in LIMITS.values())
    highest = max(row['lower_upper_bracket_on_required_scale_GeV'][1] for row in LIMITS.values())
    min_w, max_w = (bare / highest)**2, (bare / lowest)**2
    return dict(weak_mass_GeV=low_doublet, common_inverse_alpha_at_1e16=a.tolist(),
                alpha_single_inverse=a5, single_X_equivalent_GeV=mass0,
                shared_triplet_GeV=h, shared_adjoint_chiral_GeV=sigma,
                gA_gB_g=[*couplings, g], f_v_GeV=[f*mass0, v*mass0],
                X_masses_GeV=(poles*mass0).tolist(), unbroken_heavy_mass_GeV=mc*mass0,
                extra_inverse_gauge_thresholds=[float(f'{x:.10f}') for x in vector_shift],
                X_coupling_squared_by_site_and_pole=(current_couplings**2).tolist(),
                bare_current_scales_GeV=bare.tolist(), B_to_A_tree_BNV_amplitude=tree[1,1]/tree[0,0],
                B_to_A_width_at_common_running=(tree[1,1]/tree[0,0])**2,
                max_independent_relative_Wilson_variation_for_any_B_above_A=4.5/6.5,
                heavy_SM_current_squared_by_site=heavy_sm.tolist(),
                A_to_B_heavy_SM_contact=heavy_sm[0]/heavy_sm[1],
                common_Wilson_transfer_all_six_pass_below_by_site=min_w.tolist(),
                common_Wilson_transfer_all_six_fail_above_by_site=max_w.tolist())


def main():
    arg = argparse.ArgumentParser()
    arg.add_argument('--check', action='store_true')
    checking = arg.parse_args().check
    weak = float(LOW['OLD']['weak_fermions_GeV'][0])
    out = {'archived_flavour_spectrum': two_groups(weak),
           'direct_search_spectrum_sensitivity': two_groups(224.)}
    # Finite presentation precision keeps machine-dependent diagonalisation
    # noise below the saved digits while assertions use unrounded values.
    out = json.loads(json.dumps(out), parse_float=lambda x: float(f'{float(x):.11g}'))
    dest = HERE / 'independent_results.json'
    if checking:
        assert out == json.loads(dest.read_bytes())
        print('1770 independent OK: three couplings; two poles; six flavour assignments; localised currents')
    else:
        with dest.open('x', encoding='utf-8') as handle:
            json.dump(out, handle, indent=2, allow_nan=False)
            handle.write('\n')


if __name__ == '__main__':
    main()
