"""Check both flavour moments using a discrete single-spurion construction."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EDGES = [(0, 1), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2), (1, 2)]


def assemble(values):
    mass = np.zeros((3, 3), dtype=complex)
    for cell, value in zip(EDGES, values, strict=True):
        mass[cell] += value
    return mass


def measure(mass, vertex, edges):
    left, singular, vh = np.linalg.svd(mass)
    order = np.argsort(singular)
    left, right, singular = left[:, order], vh.conj().T[:, order], singular[order]
    physical = left.conj().T @ vertex @ right
    p = (physical.diagonal() / singular).imag
    squares = singular**2
    H = mass @ mass.conj().T
    ring = (H[0, 1] * H[1, 2] * H[2, 0]).imag
    trace_part = H[1, 1].real + abs(edges[4])**2
    contrasts = [squares @ p / (edges[5].real * edges[6].imag),
                 H[0, 0].real * ((squares * (squares - trace_part)) @ p) / ring]
    offdiag = physical - np.diag(physical.diagonal())
    assert abs(p.sum()) < 1e-10
    return p, contrasts, float(np.linalg.norm(offdiag)**2 / squares[-1])


def enumeration(N, eta):
    histogram = {-1: 0, 0: 0, 1: 0}
    for l2, l3, r2, r3 in itertools.product(range(N), repeat=4):
        c, e, f, r = [(l + j) % N for l, j in [(l2, r2), (l3, r2), (l3, r3), (l2, r3)]]
        m, remainder = divmod(r + e - c - f, N)
        assert remainder == 0 and m in histogram
        histogram[m] += 1
        assert (abs(np.sin(m * eta)) < 1e-10) == (m == 0)
    return list(histogram.values())


def evaluate():
    cfg = json.loads((HERE / 'inputs.json').read_bytes())
    old = json.loads((HERE / cfg['texture']).read_bytes())['reconstruction']
    bound = json.loads((HERE / cfg['inclusive_bound']).read_bytes())['quadratic_form']
    N, sine = cfg['N'], cfg['sin_theta']
    cosine = np.sqrt(1 - sine**2)
    xi = cfg['v_GeV'] * sine / cfg['f_GeV']
    eta = np.deg2rad(old['eta_example_deg'])
    phi = eta / N
    reals = np.array(old['normalized_real_edges_A_B_C_E_F'] + old['two_real_coefficients'])
    edges = reals.astype(complex)
    edges[-1] *= np.exp(1j * eta)
    mass = assemble(edges)
    powers = np.array([(cfg['left'][i] + cfg['right'][j]) % N for i, j in EDGES])
    powers[-1] += N
    left_phases = np.diag(np.exp(-1j * phi * np.array(cfg['left'])))
    right_phases = np.diag(np.exp(-1j * phi * np.array(cfg['right'])))
    # Real Wilson coefficients in units of the reference bottom Yukawa.
    coeff = reals / cfg['epsilon']**powers

    def perturbed(displacement, extra):
        # Dimensionless displacement is h/v. A neutral insertion on C is
        # normalized to one in the vacuum and has two additional radial powers.
        scale = 1 + xi * displacement
        paths = coeff * (cfg['epsilon'] * np.exp(1j * phi) * scale)**powers
        paths[2] *= scale**extra
        return (1 + cosine * displacement) * left_phases @ assemble(paths) @ right_phases

    table = {}
    for label, extra in [('holomorphic', 0), ('neutral_on_C', 2)]:
        np.testing.assert_allclose(perturbed(0, extra), mass, atol=1e-14)
        step = 1e-5
        numeric = (perturbed(step, extra) - perturbed(-step, extra)) / (2 * step)
        radial = powers.copy()
        radial[2] += extra
        vertex = cosine * mass + xi * assemble(edges * radial)
        np.testing.assert_allclose(numeric, vertex, atol=2e-10)
        p, z, fc = measure(mass, vertex, edges)
        expected = xi * np.array([N, N - extra])
        np.testing.assert_allclose(z, expected, atol=1e-10)
        fc_lower = float(expected @ bound @ expected)
        assert fc >= fc_lower - 1e-10
        table[label] = {'moments': z, 'diagonal_p': p.tolist(),
                        'inclusive_FC_over_SMbb': fc, 'inclusive_lower': fc_lower}
    return {'powers_A_B_C_E_F_R_D': powers.tolist(), 'xi': xi, 'N_xi': N * xi,
            'FC_lower_coefficient_on_equal_contrasts': float(np.sum(bound)),
            'absolute_Wilson_D_over_R': abs(coeff[-1] / coeff[-2]),
            'phase_deg': float(np.rad2deg(phi)), 'samples': table,
            'plaquette_m_minus_zero_plus': {str(j): enumeration(j, eta) for j in cfg['test_N']}}


def tidy(item):
    if isinstance(item, dict):
        return {k: tidy(v) for k, v in item.items()}
    if isinstance(item, list):
        return [tidy(v) for v in item]
    if isinstance(item, (float, np.floating)):
        return float(f'{item:.10g}')
    return item


if __name__ == '__main__':
    cli = argparse.ArgumentParser()
    cli.add_argument('--check', action='store_true')
    args = cli.parse_args()
    result = tidy(evaluate())
    outfile = HERE / 'results.json'
    if args.check:
        assert json.loads(outfile.read_bytes()) == result
        print('PASS: charge windings, physical moments, neutral insertion, inclusive bound')
    else:
        with outfile.open('x', encoding='utf-8') as stream:
            stream.write(json.dumps(result, indent=2) + '\n')
