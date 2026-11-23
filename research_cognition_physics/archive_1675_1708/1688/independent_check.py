"""Check discrete anomaly compensation, low-order seesaw zeros, and real UV interpolation."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ARCH = HERE.parent
LEPTONS = (2, 0, 1)
FIELD_CHARGES = (-1, 4, -2, -3)


def monomials(total_charge, budget):
    for powers in itertools.product(range(budget + 1), repeat=4):
        if sum(powers) <= budget and (total_charge + sum(
                q*n for q, n in zip(FIELD_CHARGES, powers))) % 12 == 0:
            yield powers


def allowed(charge, budget):
    return any(True for _ in monomials(charge, budget))


def anomaly_free(extra):
    # Left-handed chiral fermions, all multiplicities of the stated SU(5) reps.
    charges = [1]*10 + [2]*5 + [1]*5 + [3]*10 + [-1, 4, -2, -3, 2, 3]
    charges += list(extra)
    return sum(charges) % 6 == 0 and sum(q**3 for q in charges) % 36 == 0


def texture(pair):
    first = [j for j, q in enumerate(pair) if allowed(LEPTONS[0]+q, 1)]
    bare = [[allowed(l+q, 0) for q in pair] for l in LEPTONS]
    bare_rank = int(np.linalg.matrix_rank(np.array(bare, dtype=float)))
    assert min(sum(any(row[j] for row in bare) for j in range(2)),
               sum(any(row) for row in bare)) <= 1
    assert bare_rank <= 1
    assert len(first) <= 1
    if first:
        other = 1-first[0]
        assert not allowed(2*pair[other], 1)
    diag = [allowed(2*q, 1) for q in pair]
    off = allowed(sum(pair), 1)
    return dict(pair=list(pair), bare_D_rank=bare_rank, first_row_columns=first,
                both_heavy_masses_possible=bool(all(diag) or off))


def interpolation_basis(charge, vevs):
    for budget in range(5):
        options = list(monomials(charge, budget))
        # Best conditioned phase pair, favor the larger amplitudes.
        best = None
        for p, q in itertools.combinations(options, 2):
            u, v = (np.prod(vevs**powers) for powers in (p, q))
            determinant = abs((u.conjugate()*v).imag)
            if determinant > 1e-12 and (best is None or determinant > best[0]):
                best = (determinant, p, q, u, v)
        if best is not None:
            _, p, q, u, v = best
            return budget, (p, q), np.array([[u.real, v.real], [u.imag, v.imag]])
    raise AssertionError(f'No two independent phases for charge {charge}')


def check():
    single = [q for q in range(12) if anomaly_free([q])]
    pairs = [p for p in itertools.combinations_with_replacement(range(12), 2) if anomaly_free(p)]
    assert single == [2, 8] and pairs == [(0, 2), (0, 8), (2, 6), (3, 5), (6, 8), (9, 11)]
    records = [texture(pair) for pair in pairs]

    setup = json.loads((ARCH/'1685/inputs.json').read_bytes())
    phi = setup['scalar_angle_over_pi']*math.pi
    s = setup['scalar_magnitude_over_cutoff']*np.exp(1j*phi)
    x = setup['real_X_over_cutoff']
    vevs = np.array([s, x, s*s, s**3])
    # Pair (0,2): M00 bare, M22 via X^2, M02 has zero coefficient.
    heavy = np.diag([1., x*x])
    simple = np.array([[s*s, 0], [0, s*s], [0, s**3]])
    coarse_mass = -simple @ np.linalg.solve(heavy, simple.T)
    assert np.linalg.matrix_rank(coarse_mass) == 2
    np.testing.assert_allclose(coarse_mass[0, 0], -s**4, atol=2e-15)

    witness = json.loads((ARCH/'1682/independent_results.json').read_bytes())['witnesses'][0]
    phase = np.exp(4j*phi)
    legs = np.array([1, np.exp(2j*phi), witness['neutral_modulus'] *
                     np.exp(3j*phi)*(3+2*phase)/abs(3+2*phase)])
    a = legs * [0, 1, 1]
    b = legs * [1, 3, 1]
    r = witness['source_r']
    eta = 2*math.pi/3
    target = np.outer(a, a) + r*np.exp(1j*eta)*np.outer(b, b)
    columns = np.column_stack([1j*math.sqrt(r)*np.exp(.5j*eta)*b, 1j*x*a])
    actual = np.zeros_like(columns)
    bases = {k: interpolation_basis(k, vevs) for k in range(5)}
    max_coefficient = 0.
    for i, lepton in enumerate(LEPTONS):
        for j, n in enumerate((0, 2)):
            _, powers, real_map = bases[lepton + n]
            desired = columns[i, j]
            coefs = np.linalg.solve(real_map, [desired.real, desired.imag])
            actual[i, j] = sum(c*np.prod(vevs**p) for c, p in zip(coefs, powers))
            max_coefficient = max(max_coefficient, float(max(abs(coefs))))
    mass = -actual @ np.linalg.solve(heavy, actual.T)
    error = float(abs(mass-target).max())
    assert error < 1e-12 and np.linalg.matrix_rank(mass) == 2
    return dict(single_compensators=single, pairs_and_texture=records,
                two_insertion_rank_two_weak_11_modulus=round(float(abs(coarse_mass[0, 0])), 12),
                minimum_insertions_for_arbitrary_complex_entry_by_charge={
                    str(k): v[0] for k, v in bases.items()},
                real_coefficient_fit_max_abs=round(max_coefficient, 9),
                recovered_matrix_max_abs_error_below_1e_12=True,
                frozen_1682_normalized_weak_11_modulus=round(float(abs(target[0, 0])), 12))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    result = check()
    path = HERE/'independent_results.json'
    if parser.parse_args().check:
        assert result == json.loads(path.read_bytes())
        print('PASS: all anomaly pairs, weak-frame zero at one insertion, real two- and four-insertion lifts')
    else:
        with path.open('x', encoding='utf8') as f:
            f.write(json.dumps(result, indent=2)+'\n')
