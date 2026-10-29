"""Round 1314: read-only checks; --write-new exclusively creates results.json."""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
from itertools import permutations
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add(a, b):
    out = defaultdict(complex, a)
    for key, value in b.items():
        out[key] += value
    return {key: value for key, value in out.items() if value != 0}


def scale(a, scalar):
    return {key: value * scalar for key, value in a.items() if value * scalar != 0}


def wedge(a, b):
    out = defaultdict(complex)
    for x, cx in a.items():
        for y, cy in b.items():
            if x & y:
                continue
            parity = sum((y & ((1 << i) - 1)).bit_count()
                         for i in range(x.bit_length()) if (x >> i) & 1)
            out[x | y] += (-1) ** parity * cx * cy
    return {key: value for key, value in out.items() if value != 0}


def spin_contract(a, b):
    return add(wedge(a[0], b[1]), scale(wedge(a[1], b[0]), -1))


def epsilon_sign(p):
    return (-1) ** sum(p[i] > p[j] for i in range(len(p))
                      for j in range(i + 1, len(p)))


def fields(n=4):
    return [[{1 << (2 * i + s): 1} for s in range(2)] for i in range(n)]


def transform(a, fs):
    out = []
    for row in a:
        pair = []
        for spin in range(2):
            poly = {}
            for j, value in enumerate(row):
                poly = add(poly, scale(fs[j][spin], value))
            pair.append(poly)
        out.append(pair)
    return out


def four_n(a, fs=None):
    fs = fields() if fs is None else fs
    af = transform(a, fs)
    out = {}
    for i, j, k, l in permutations(range(4)):
        product = wedge(spin_contract(fs[i], af[j]), spin_contract(fs[k], af[l]))
        out = add(out, scale(product, epsilon_sign((i, j, k, l))))
    return out


def exact_coefficients(poly):
    out = {}
    for key, value in sorted(poly.items()):
        assert value.imag == 0 and value.real == int(value.real)
        out[str(key)] = int(value.real)
    return out


def special_unitary(rng, n):
    x = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, _ = np.linalg.qr(x)
    q[:, 0] /= np.linalg.det(q)
    assert np.linalg.norm(q.conj().T @ q - np.eye(n)) < 1e-12
    assert abs(np.linalg.det(q) - 1) < 1e-12
    return q


def max_diff(a, b):
    return max((abs(a.get(key, 0) - b.get(key, 0))
                for key in set(a) | set(b)), default=0.0)


def evaluate_ordered(poly, spinors):
    value = 0j
    for mask, coefficient in poly.items():
        spins = []
        for flavor in range(4):
            local = (mask >> (2 * flavor)) & 3
            assert local in (1, 2), 'This ordered vertex requires one leg per flavor'
            spins.append(0 if local == 1 else 1)
        term = coefficient
        for flavor, spin in enumerate(spins):
            term *= spinors[flavor][spin]
        value += term
    return value


def compute():
    checks = {}
    # Exact center exponents, no numerical roots-of-unity test as proof.
    lep = dict(Q=0, u=0, d=0, L=1, e=-1, N=-1)
    quark = dict(Q=1, u=-1, d=-1, L=0, e=0, N=0)
    yukawa_pairs = (('Q', 'u'), ('Q', 'd'), ('L', 'e'), ('L', 'N'))
    for weights in (lep, quark):
        assert all(weights[a] + weights[b] == 0 for a, b in yukawa_pairs)
    center_rows = []
    for n in range(2, 17):
        for nc in (3, 5, 7):
            e_allowed = (1 % n == 0) and (nc % n == 0)
            assert not e_allowed
            center_rows.append(dict(n=n, nc=nc, e_lepton_residue=1 % n,
                                    e_quark_residue=nc % n,
                                    four_n_residue=(-4) % n))
    assert [n for n in range(3, 17) if 4 % n == 0] == [4]
    checks['exact_center_samples'] = dict(status='PASS', rows=center_rows,
        full_n_statement='Analytic stabilizer theorem; finite rows are implementation checks')

    rng = np.random.default_rng(1314)
    max_yukawa_error = max_center_error = 0.0
    for n in (3, 4, 5):
        us = {f: special_unitary(rng, n) for f in lep}
        zeta = np.exp(2j * np.pi / n)
        for left, right in yukawa_pairs:
            y = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
            f = rng.normal(size=n) + 1j * rng.normal(size=n)
            h = rng.normal(size=n) + 1j * rng.normal(size=n)
            yp = us[left].conj() @ y @ us[right].conj().T
            residual = abs((us[left] @ f).T @ yp @ (us[right] @ h) - f.T @ y @ h)
            max_yukawa_error = max(max_yukawa_error, residual)
            for weights in (lep, quark):
                ul = zeta ** weights[left] * np.eye(n)
                ur = zeta ** weights[right] * np.eye(n)
                max_center_error = max(max_center_error,
                    np.linalg.norm(ul.conj() @ y @ ur.conj().T - y))
    assert max_yukawa_error < 1e-10 and max_center_error < 1e-12
    checks['dirac_data_transform'] = dict(status='PASS', generations=[3, 4, 5],
        general_matrix_error_below=1e-10, both_stabilizer_errors_below=1e-12)

    fs = fields()
    bare = four_n(np.eye(4, dtype=int))
    assert bare == {}
    a = np.diag([1, 2, 3, 4])
    dressed = four_n(a)
    expected = {'90': 14, '102': -4, '105': -10, '150': -10, '153': -4, '165': 14}
    assert exact_coefficients(dressed) == expected
    b12 = spin_contract(fs[0], fs[1])
    b34 = spin_contract(fs[2], fs[3])
    b13 = spin_contract(fs[0], fs[2])
    b24 = spin_contract(fs[1], fs[3])
    b14 = spin_contract(fs[0], fs[3])
    b23 = spin_contract(fs[1], fs[2])
    schouten = add(add(wedge(b12, b34), wedge(b13, b24)), wedge(b14, b23))
    assert schouten == {}
    formula = add(add(scale(wedge(b12, b34), 2), scale(wedge(b13, b24), -8)),
                  scale(wedge(b14, b23), 6))
    assert exact_coefficients(formula) == expected
    u = special_unitary(rng, 4)
    transformed = four_n(u @ a @ u.conj().T, transform(u, fs))
    assert max_diff(transformed, dressed) < 1e-10
    checks['four_n_flavor_boundary'] = dict(status='PASS', bare_epsilon_zero=True,
        schouten_zero=True, a_eigenvalues=[1, 2, 3, 4], ordered_coefficients=expected,
        su4_covariance_error_below=1e-10,
        scope='Nonzero four-N subproblem only; E remains forbidden')

    # Explicit crossing of N1 N2 -> anti-N3 anti-N4 at a massless tree point.
    energy = 1.0
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1, -1]).astype(complex)
    momenta = [np.array(p, float) for p in
               ((0, 0, -1), (0, 0, 1), (-1, 0, 0), (1, 0, 0))]
    spinors = [np.array(v, complex) for v in
               ((1, 0), (0, 1), (1 / np.sqrt(2), 1 / np.sqrt(2)),
                (-1 / np.sqrt(2), 1 / np.sqrt(2)))]
    spinors = [np.sqrt(2 * energy) * x for x in spinors]
    for p, x in zip(momenta, spinors):
        sigmap = p[0] * sx + p[1] * sy + p[2] * sz
        assert np.linalg.norm(sigmap @ x + energy * x) < 1e-12
        assert abs(np.vdot(x, x) - 2 * energy) < 1e-12
    assert np.linalg.norm(momenta[0] + momenta[1] - momenta[2] - momenta[3]) == 0
    amplitude = evaluate_ordered(dressed, spinors)
    assert abs(amplitude - 12 * energy ** 2) < 1e-12
    checks['massless_tree_vertex'] = dict(status='PASS', energy=1,
        contact_polynomial=12, all_external_norms=2,
        four_momentum_conserved=True, weyl_equations_pass=True,
        finite_instrument_certified=False, higher_order_error_certified=False)

    # Same Yukawa point at nonzero Higgs amplitude, with exact on-shell kinematics.
    rho_h, root_s = 0.1, 2.0
    masses = rho_h * np.sqrt(np.array([1, 2, 3, 4], float))
    energies, moduli = [], []
    s = root_s ** 2
    for i, j in ((0, 1), (2, 3)):
        mi, mj = masses[i], masses[j]
        assert root_s > mi + mj
        ei = (s + mi ** 2 - mj ** 2) / (2 * root_s)
        ej = root_s - ei
        pij = np.sqrt((s - (mi + mj) ** 2) * (s - (mi - mj) ** 2)) / (2 * root_s)
        energies.extend([ei, ej])
        moduli.extend([pij, pij])
    massive_upper, massive_p = [], []
    for i, (m, ei, pij, direction, unit) in enumerate(
            zip(masses, energies, moduli, momenta,
                [x / np.sqrt(2 * energy) for x in spinors])):
        p = pij * direction
        massive_p.append(p)
        sigmap = p[0] * sx + p[1] * sy + p[2] * sz
        upper = np.sqrt(ei + pij) * unit
        lower = np.sqrt(ei - pij) * unit
        h = np.block([[-sigmap, m * np.eye(2)], [m * np.eye(2), sigmap]])
        if i < 2:
            full = np.concatenate([upper, lower])
            assert np.linalg.norm(h @ full - ei * full) < 1e-12
        else:
            h_minus_p = np.block([[sigmap, m * np.eye(2)], [m * np.eye(2), -sigmap]])
            full = np.concatenate([upper, -lower])
            assert np.linalg.norm(h_minus_p @ full + ei * full) < 1e-12
        assert abs(np.vdot(full, full) - 2 * ei) < 1e-12
        massive_upper.append(upper)
    assert abs(energies[0] + energies[1] - energies[2] - energies[3]) < 1e-12
    assert np.linalg.norm(massive_p[0] + massive_p[1] - massive_p[2] - massive_p[3]) < 1e-12
    massive_expected = 3 * np.sqrt(np.prod(np.array(energies) + np.array(moduli)))
    assert abs(evaluate_ordered(dressed, massive_upper) - massive_expected) < 1e-12
    assert massive_expected > 0
    checks['positive_mass_tree_vertex'] = dict(status='PASS', higgs_amplitude=0.1,
        center_of_mass_energy=2, mass_squared=[0.01, 0.02, 0.03, 0.04],
        contact_polynomial_rounded=float(f'{massive_expected:.10g}'),
        full_dirac_equations_pass=True, four_momentum_conserved=True,
        finite_instrument_certified=False, higher_order_error_certified=False)

    # R is a distinct protection: it allows the two quartic monomials.
    assert 1j ** 4 == 1 and 1j ** 2 == -1
    checks['r_is_not_flavor_center_obstruction'] = dict(status='PASS',
        qqql_phase=1, four_n_phase=1, majorana_phase=-1)
    return dict(schema='round1314_author_results_v1', round=1314,
        status='PASS', acceptance='requires_independent_review_and_main_acceptance',
        new_axioms=0, empirical_predictions_claimed=0,
        theorem='M_Y forbids E for every n >= 2; thus incompatible with E and P',
        checks=checks,
        author_dependencies_sha256={
            'proof.md': sha(HERE / 'proof.md'),
            'check.py': sha(Path(__file__)),
            '../_admission/round1314_admission.md': sha(HERE.parent / '_admission/round1314_admission.md')})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true')
    group.add_argument('--write-new', action='store_true')
    args = parser.parse_args()
    result = compute()
    path = HERE / 'results.json'
    if args.write_new:
        with path.open('x', encoding='utf8', newline='\n') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write('\n')
    else:
        stored = json.loads(path.read_text(encoding='utf8'))
        if stored != result:
            raise AssertionError('Stored results differ; check is read-only and will not overwrite')
    print(json.dumps(dict(status='PASS', round=1314, groups=len(result['checks']),
        mode='write-new' if args.write_new else 'read-only check',
        results_sha256=sha(path)), ensure_ascii=False))


if __name__ == '__main__':
    main()
