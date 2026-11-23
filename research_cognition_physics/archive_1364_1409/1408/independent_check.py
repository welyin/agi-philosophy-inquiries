"""Verify the two-heavy determinant identity and a flavour-coherence holdout."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def check_complex_five_state_mass_matrices():
    rng = np.random.default_rng(1408)
    h = np.array([2., -1., 1.]) / np.sqrt(6)
    plane = np.eye(3) - np.outer(h, h)
    passed = 0
    for _ in range(64):
        d = plane @ (rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2)))
        raw = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        heavy = np.diag([7., 11.]) + (raw + raw.T) / 2
        mass = np.block([[np.zeros((3, 3)), d], [d.T, heavy]])
        left, singular, _ = np.linalg.svd(mass)
        order = np.argsort(singular)
        masses = singular[order]
        left = left[:, order]
        diagonal = np.diag(left.conj().T @ mass @ left.conj())
        left = left * np.exp(0.5j * np.angle(diagonal))
        np.testing.assert_allclose((left * masses) @ left.T, mass, atol=2e-12)
        s = left[:3, 3:]
        deficit = s @ s.conj().T
        eigenvalues = np.linalg.eigvalsh(deficit)[1:]
        np.testing.assert_allclose(deficit @ h, 0, atol=1e-13)
        np.testing.assert_allclose(np.prod(eigenvalues / (1-eigenvalues)),
                                   np.prod(masses[1:3])/np.prod(masses[3:]), rtol=2e-9)
        for i, j, k in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
            minor = deficit[i, i].real * deficit[j, j].real - abs(deficit[i, j])**2
            np.testing.assert_allclose(minor, np.prod(eigenvalues) * abs(h[k])**2,
                                       atol=2e-15)
        passed += 1
    return passed


def certificate():
    old = json.loads((HERE.parent / '1370' / 'inputs.json').read_text(encoding='utf-8-sig'))
    m2, m3 = [math.sqrt(old['NO'][field][0]) for field in ('solar', 'atm_abs')]
    threshold_eV = 1.e11
    ratio = m2 * m3 / threshold_eV**2
    p = 0.0006
    # The matrices p*P+i*c*J have exactly the same trace and row deficits.
    # Use the determinant directly: subtracting c**2 from p**2 loses precision.
    product = ratio * (1 - 2*p) / (1 - ratio)
    lambda_small = product / (p + math.sqrt(p*p-product))
    lambda_large = 2*p - lambda_small
    de, dmu = p/3, 5*p/6
    appearance_denominator = (1-de)*(1-dmu)
    minor = product/6
    relative_loss = minor/(de*dmu)
    no_imaginary_part = p*p/9/appearance_denominator
    return {
        'random_five_state_takagi_checks': check_complex_five_state_mass_matrices(),
        'central_mass_eV': [m2, m3], 'common_heavy_mass_design_eV': threshold_eV,
        'odds_product': ratio, 'all_real_plane_row_deficits': [de, dmu, dmu],
        'heavy_example_eigenvalues': [lambda_large, lambda_small],
        'heavy_example_emu_minor': minor,
        'fractional_shortfall_from_saturated_abs_Xemu_squared': relative_loss,
        'coherent_zero_length_probability_conditional': (de*dmu-minor)/appearance_denominator,
        'same_diagonals_and_null_real_plane_competitor_probability': no_imaginary_part,
        'real_plane_competitor_inferred_equal_heavy_mass_eV': math.sqrt(m2*m3)*(1-p)/p,
    }


def verify(actual, saved):
    if actual.keys() != saved.keys():
        raise ValueError('Keys differ')
    for name, value in actual.items():
        np.testing.assert_allclose(value, saved[name], rtol=2e-11, atol=0)


if __name__ == '__main__':
    args = argparse.ArgumentParser()
    args.add_argument('--write', action='store_true')
    parsed = args.parse_args()
    target = HERE / 'independent_results.json'
    result = certificate()
    if parsed.write:
        with target.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write('\n')
    else:
        verify(result, json.loads(target.read_text(encoding='utf-8')))
    print(json.dumps(result, allow_nan=False))
