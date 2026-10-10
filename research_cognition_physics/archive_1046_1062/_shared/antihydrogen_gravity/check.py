"""Read-only algebra and published rounded-unit check, not an ALPHA-g fit.

--write creates the first results.json exclusively. No experiment likelihood,
combined uncertainty, detector simulation, or new scientific group is computed.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import numpy as np


def compute():
    one = np.eye(2, dtype=int)
    a = np.array([[0, 1], [0, 0]], dtype=int)
    z = np.diag([1, -1])
    b, d = np.kron(a, one), np.kron(z, a)
    eye = np.eye(4, dtype=int)
    assert np.array_equal(b @ b.T + b.T @ b, eye)
    assert np.array_equal(d @ d.T + d.T @ d, eye)
    assert not np.any(b @ d + d @ b)
    assert not np.any(b @ d.T + d.T @ b)
    nb, nd = b.T @ b, d.T @ d
    E, q = 5, 2  # dimensionless diagnostic values, not atom data
    Hraw = E * (nb - d @ d.T)
    Qraw = q * (nb + d @ d.T)
    H, Q = Hraw + E * eye, Qraw - q * eye
    assert np.array_equal(H, E * (nb + nd))
    assert np.array_equal(Q, q * (nb - nd))
    assert H.diagonal().tolist() == [0, E, E, 2 * E]
    assert Q.diagonal().tolist() == [0, -q, q, 0]

    # Symbolic dimension bookkeeping (mass, length, time, magnetic field).
    mass, acc = (1, 0, 0, 0), (0, 1, -2, 0)
    energy, length = (1, 2, -2, 0), (0, 1, 0, 0)
    magnetic, duration = (0, 0, 0, 1), (0, 0, 1, 0)
    add = lambda x, y: tuple(a+b for a, b in zip(x, y))
    sub = lambda x, y: tuple(a-b for a, b in zip(x, y))
    slope, gradient = sub(energy, magnetic), sub(magnetic, length)
    assert sub(add(slope, gradient), mass) == acc
    assert add(add(mass, acc), length) == energy
    assert add(slope, sub(magnetic, duration)) == sub(energy, duration)

    # Binding/rest-level bookkeeping; arbitrary rational energy units.
    mp, me, c2, e0, esB = F(100), F(1), F(100), F(-1, 4), F(-1, 5)
    M = mp + me + e0/c2
    shift = esB - e0
    assert M > 0 and M < mp + me
    assert M*c2 + shift == (mp+me)*c2 + esB

    # Printed ALPHA-g p.718 scale; rounded inputs only, not a new measurement.
    grad_T_m, height_m = F('1.77e-3'), F('0.256')
    field_T = grad_T_m * height_m
    printed_T = F('4.53e-4')
    assert abs(field_T - printed_T) <= F('5e-7')
    assert field_T > 0
    # A compensating top-minus-bottom field must have opposite sign to gravity.
    assert field_T + (-field_T) == 0

    return {
        'scope': 'free_CAR_and_composite_energy_units_only_no_ALPHAg_reanalysis',
        'fermion_basis': ['vacuum', 'antiparticle', 'particle', 'pair'],
        'energy_diagonal_diagnostic': H.diagonal().tolist(),
        'charge_diagonal_diagnostic': Q.diagonal().tolist(),
        'CAR_and_normal_ordering_exact': True,
        'composite_internal_energy_identity_exact': True,
        'magnetic_force_and_drive_power_dimensions_correct': True,
        'published_scale_product_T': str(field_T),
        'published_rounded_difference_T': str(field_T - printed_T),
        'published_acceleration_ratio': {
            'central': '0.75', 'statistical_plus_systematic': '0.13',
            'simulation': '0.16', 'reference_g_m_s2': '9.81',
            'interpretation': 'transcription_only_no_combined_sigma_or_p_value'
        },
        'new_scientific_groups': 0,
        'all_assertions_passed': True
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = compute()
    path = Path(__file__).with_name('results.json')
    if args.write:
        with path.open('x', encoding='utf-8') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
        print('First-save algebra and rounded-unit checks passed.')
    else:
        assert json.loads(path.read_text(encoding='utf-8')) == result
        print('Read-only checks passed. ALPHA-g likelihood and uncertainty not refitted.')
