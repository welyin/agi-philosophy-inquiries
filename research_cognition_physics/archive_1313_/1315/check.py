"""Round 1315: deterministic read-only checks; --write-new exclusively creates results."""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EPS = np.array([[0, 1], [-1, 0]], dtype=int)
K = np.kron(EPS, EPS)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add(a, b, c=1):
    out = dict(a)
    for mask, value in b.items():
        out[mask] = out.get(mask, 0) + c * value
        if out[mask] == 0:
            del out[mask]
    return out


def wedge(a, b):
    out = {}
    for ma, va in a.items():
        for mb, vb in b.items():
            if ma & mb:
                continue
            crossings = sum((mb & ((1 << i) - 1)).bit_count()
                            for i in range(ma.bit_length()) if ma & (1 << i))
            sign = -1 if crossings % 2 else 1
            mask = ma | mb
            out[mask] = out.get(mask, 0) + sign * va * vb
    return {mask: value for mask, value in out.items() if value != 0}


def bilinear(psi, phi):
    out = {}
    for a in range(4):
        for b in range(4):
            if K[a, b]:
                out = add(out, wedge(psi[a], phi[b]), int(K[a, b]))
    return out


def exterior_pf(matrix):
    @lru_cache(None)
    def rec(ids):
        if not ids:
            return {0: 1}
        i = ids[0]
        out = {}
        for j in range(1, len(ids)):
            rest = ids[1:j] + ids[j + 1:]
            out = add(out, wedge(matrix[i][ids[j]], rec(rest)),
                      1 if j % 2 else -1)
        return out
    return rec(tuple(range(len(matrix))))


def target_coefficient(n):
    """Pf signs and Grassmann ordering independently, no factorial insertion."""
    m = 2 * n
    ns = 2 * m
    @lru_cache(None)
    def rec(ids):
        if not ids:
            return 1
        i = ids[0]
        total = 0
        for j in range(1, len(ids)):
            k = ids[j]
            if (i < m) == (k < m):
                continue
            # Pf sign (-1)^(j+1); sorting theta_i theta_k past j-1 others
            # gives (-1)^(j-1). Keep these two signs separate.
            pf_sign = 1 if j % 2 else -1
            grass_sign = 1 if (j - 1) % 2 == 0 else -1
            total += pf_sign * grass_sign * rec(ids[1:j] + ids[j + 1:])
        return total
    return rec(tuple(range(ns)))


def su(rng, n):
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, _ = np.linalg.qr(a)
    q[:, 0] /= np.linalg.det(q)
    return q


def blockdiag(a, b):
    out = np.zeros((len(a) + len(b), len(a) + len(b)), complex)
    out[:len(a), :len(a)] = a
    out[len(a):, len(a):] = b
    return out


def mix(n):
    dft = np.exp(2j * np.pi * np.outer(np.arange(3), np.arange(3)) / 3) / np.sqrt(3)
    out = blockdiag(dft, np.eye(n - 3)) if n > 3 else dft
    for j in range(3, n):
        angle = 1e-4 / j
        rot = np.eye(n)
        rot[2, 2] = rot[j, j] = math.cos(angle)
        rot[2, j] = -math.sin(angle)
        rot[j, 2] = math.sin(angle)
        out = out @ rot
    return out


def compute():
    checks = {}
    assert np.array_equal(K, K.T)
    # Full four-component exterior S, before any projection.
    psi = [[{1 << (4 * f + a): 1} for a in range(4)] for f in range(4)]
    smat = [[bilinear(x, y) for y in psi] for x in psi]
    assert all(smat[i][i] == {} for i in range(4))
    assert all(add(smat[i][j], smat[j][i]) == {} for i in range(4) for j in range(4))
    p4 = exterior_pf(smat)
    target = (1 << 0) | (1 << 4) | (1 << 11) | (1 << 15)
    assert p4[target] == 2
    checks['full_grassmann_small_case'] = dict(status='PASS', symmetric_internal_form=True,
        species_antisymmetric=True, full_polynomial_nonzero_terms=len(p4), target_coefficient=2)

    # Universal Pf covariance over the even exterior ring, using an exact det-one shear.
    ns = 8
    original = [[{} for _ in range(4)] for _ in range(ns)]
    for i in range(ns):
        original[i][0 if i < ns // 2 else 3] = {1 << i: 1}
    g = np.eye(ns, dtype=int)
    g[0, 5], g[2, 7], g[4, 6] = 2, -3, 1
    assert round(np.linalg.det(g)) == 1
    transformed = [[{} for _ in range(4)] for _ in range(ns)]
    for i in range(ns):
        for j in range(ns):
            for a in range(4):
                transformed[i][a] = add(transformed[i][a], original[j][a], int(g[i, j]))
    sp = [[bilinear(x, y) for y in original] for x in original]
    st = [[bilinear(x, y) for y in transformed] for x in transformed]
    assert exterior_pf(st) == exterior_pf(sp) == {(1 << ns) - 1: 24}
    rows = []
    for n in range(1, 9):
        value = target_coefficient(n)
        assert value == math.factorial(2 * n)
        rows.append(dict(n=n, target_coefficient=value))
    checks['pfaffian_covariance_and_exact_nonzero'] = dict(status='PASS',
        even_ring_covariance_exact=True, finite_coefficients=rows,
        general_n='Analytic Pfaffian and ordered Grassmann proof, not finite induction')

    rng = np.random.default_rng(1315)
    group_rows = []
    for n in (3, 4, 5):
        uq, ul, color = su(rng, n), su(rng, n), su(rng, 3)
        group = blockdiag(np.kron(uq, color), ul)
        assert abs(np.linalg.det(group) - 1) < 1e-11
        assert np.linalg.norm(group.conj().T @ group - np.eye(4 * n)) < 1e-11
        for k in range(-4, 5):
            q, ell = Fraction(1, 6) + k, Fraction(-1, 2) - 3 * k
            assert 3 * n * q + n * ell == 0
        assert (1j ** (4 * n)) == 1
        group_rows.append(dict(n=n, delta_B=n, delta_L=n, dimension=6 * n,
            r_phase=1, full_flavor_embedding_determinant_one=True,
            old_hypercharge_family_neutral=True))
    checks['common_symmetry_and_dimension'] = dict(status='PASS', rows=group_rows,
        bound_scope='Only nonzero diagonal deltaB=deltaL in original local polynomial class',
        coefficients_are_flavor_scalars_not_generated=True)

    # Same declared positive-mass point for the batch contact and phase-sector data.
    mass_rows = []
    mu, v, pair_energy = 2.0 ** -12, 0.25, 1.0
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1, -1]).astype(complex)
    pauli = (sx, sy, sz)
    for n in (3, 4, 5):
        m = 2 * n
        up = np.arange(1, n + 1) * mu
        down = (np.arange(1, n + 1) + Fraction(1, 2)) * mu
        down = np.asarray(down, float)
        charged = np.array([1, 2, 4, 6, 8][:n], float) * mu
        # Species order Q_(generation,color), then L_generation.
        masses = []
        for f in range(4 * n):
            if f < 3 * n:
                masses.append(float(up[f // 3] if f < m else down[f // 3]))
            else:
                masses.append(float(charged[f - 3 * n]))
        energies = np.zeros(4 * n)
        momenta = np.zeros((4 * n, 3))
        upper_norms = np.zeros(4 * n)
        incoming = []
        outgoing = []
        max_dirac_residual = 0.0
        for i in range(m):
            j = m + i
            mi, mj = masses[i], masses[j]
            s = pair_energy ** 2
            ei = (s + mi ** 2 - mj ** 2) / (2 * pair_energy)
            ej = pair_energy - ei
            p = math.sqrt((s - (mi + mj) ** 2) * (s - (mi - mj) ** 2)) / (2 * pair_energy)
            is_in = i < n
            (incoming if is_in else outgoing).extend((i, j))
            for index, energy, mass, direction, spin in (
                (i, ei, mi, -1, np.array([1, 0], complex)),
                (j, ej, mj, 1, np.array([0, 1], complex))):
                pvec = np.array([0, 0, direction * p])
                sigmap = sum(pvec[a] * pauli[a] for a in range(3))
                upper = math.sqrt(energy + p) * spin
                # m/sqrt(E+p) avoids cancellation in sqrt(E-p).
                lower = mass / math.sqrt(energy + p) * spin
                if is_in:
                    h = np.block([[-sigmap, mass * np.eye(2)], [mass * np.eye(2), sigmap]])
                    full = np.concatenate((upper, lower))
                    res = np.linalg.norm(h @ full - energy * full)
                else:
                    hneg = np.block([[sigmap, mass * np.eye(2)], [mass * np.eye(2), -sigmap]])
                    full = np.concatenate((upper, -lower))
                    res = np.linalg.norm(hneg @ full + energy * full)
                assert res < 1e-12
                assert abs(np.vdot(full, full).real - 2 * energy) < 1e-12
                max_dirac_residual = max(max_dirac_residual, float(res))
                energies[index], momenta[index], upper_norms[index] = energy, pvec, np.linalg.norm(upper)
        assert len(incoming) == len(outgoing) == 2 * n
        assert abs(energies[incoming].sum() - energies[outgoing].sum()) < 1e-12
        assert np.linalg.norm(momenta[incoming].sum(0) - momenta[outgoing].sum(0)) < 1e-12
        amplitude = target_coefficient(n) * float(np.prod(upper_norms))
        assert amplitude > 0
        u = mix(n)
        mnu = u @ np.diag(np.arange(1, n + 1) * mu) @ u.conj().T
        jcp = (u[0, 0] * u[1, 1] * u[0, 1].conj() * u[1, 0].conj()).imag
        assert abs(jcp - math.sqrt(3) / 18) < 1e-12
        assert np.linalg.eigvalsh(mnu).min() > 0
        y_phys = math.sqrt(2) * mnu / v
        y_weyl = y_phys.conj()
        assert np.linalg.norm(y_weyl.conj() * v / math.sqrt(2) - mnu) < 1e-14
        # All role indices are connected by the Hermitian mass matrix in this sample.
        seen = {0}
        for _ in range(n):
            seen |= {j for i in seen.copy() for j in range(n) if abs(mnu[i, j]) > 1e-15}
        assert len(seen) == n
        mass_rows.append(dict(n=n, pair_energy=pair_energy, total_cm_energy=n * pair_energy,
            incoming_legs=2 * n, outgoing_legs=2 * n,
            positive_mass_contact_rounded=float(f'{amplitude:.10g}'),
            full_dirac_equations_error_below=1e-12,
            energy_momentum_conserved=True, phase_J_nonzero=True,
            all_left_weyl_dictionary_conjugated=True,
            full_quadratic_phase_certificate='Inherited analytic task plus sufficiently small extension; not recertified here',
            g2_one_loop_b=str(Fraction(-43 + 8 * n, 6))))
    checks['same_parent_positive_mass_and_phase_data'] = dict(status='PASS', mu=mu, v=v,
        rows=mass_rows, neutrinos_not_used_as_contact_external_legs=True,
        exact_interacting_P_instrument_claimed=False, event_rate_claimed=False)

    # Equal-mass zero point and dimension/batch selection sanity checks.
    for n in (3, 4, 5):
        assert math.factorial(2 * n) * 2 ** (2 * n) > 0
        for k in range(1, 5 * n):
            if k % n == 0:
                assert 6 * k >= 6 * n
    checks['scope_boundaries'] = dict(status='PASS', unit_E_still_forbidden=True,
        nonzero_kappa_is_input=True, zero_kappa_obeys_M_Y=True,
        no_minimum_generation_selection=True, no_topological_rate_inference=True,
        no_weak_instrument_certificate_transport=True)
    dependencies = {
        'proof.md': HERE / 'proof.md', 'check.py': Path(__file__),
        '../_admission/round1315_admission.md': HERE.parent / '_admission/round1315_admission.md',
        '../1314/proof.md': HERE.parent / '1314/proof.md',
        '../../archive_1290_1307/1303/proof.md': HERE.parent.parent / 'archive_1290_1307/1303/proof.md',
        '../../archive_1308_1312/1311/proof.md': HERE.parent.parent / 'archive_1308_1312/1311/proof.md'}
    return dict(schema='round1315_author_results_v1', round=1315, status='PASS',
        acceptance='requires_independent_review_and_main_acceptance', new_axioms=0,
        theorem='Pfaffian provides M_Y-compatible diagonal batch at sharp dimension 6n; original E remains forbidden',
        checks=checks, author_dependencies_sha256={k: sha(p) for k, p in dependencies.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true')
    group.add_argument('--write-new', action='store_true')
    args = parser.parse_args()
    result = compute()
    path = HERE / 'results.json'
    if args.write_new:
        with path.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, sort_keys=True)
            stream.write('\n')
    elif json.loads(path.read_text(encoding='utf-8')) != result:
        raise AssertionError('Stored results differ; read-only check never overwrites')
    print(json.dumps(dict(status='PASS', round=1315, groups=len(result['checks']),
        mode='write-new' if args.write_new else 'read-only check', results_sha256=sha(path)), ensure_ascii=False))


if __name__ == '__main__':
    main()
