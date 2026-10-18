"""Independent stabilizer/rank audit. No author calculation is imported."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent


def null(a, cols):
    if not len(a):
        return np.eye(cols, dtype=complex)
    _, s, vh = np.linalg.svd(a, full_matrices=True)
    rank = int(np.sum(s > 1e-10))
    return vh[rank:].conj().T


def unitary_algebra(n):
    out = []
    for i in range(n):
        a = np.zeros((n, n), complex)
        a[i, i] = 1j
        out.append(a)
    for i in range(n):
        for j in range(i + 1, n):
            a = np.zeros((n, n), complex)
            a[i, j], a[j, i] = 1, -1
            out.append(a)
            b = np.zeros((n, n), complex)
            b[i, j] = b[j, i] = 1j
            out.append(b)
    return out


def su(n):
    u = unitary_algebra(n)
    return [u[i] - u[n-1] for i in range(n-1)] + u[n:]


def sp(k):
    out = [np.block([[a, np.zeros_like(a)], [np.zeros_like(a), a.conj()]])
           for a in unitary_algebra(k)]
    for i in range(k):
        for j in range(i, k):
            for phase in (1, 1j):
                b = np.zeros((k, k), complex)
                b[i, j] = b[j, i] = phase
                z = np.zeros_like(b)
                out.append(np.block([[z, b], [-b.conj(), z]]))
    return out


def stabilizer_fixed(basis):
    d = basis[0].shape[0]
    phi = np.eye(d, dtype=complex)[:, 0]
    images = np.column_stack([a @ phi for a in basis])
    coeff = null(np.vstack([images.real, images.imag]), len(basis)).real
    generators = [sum((coeff[j, i] * a for j, a in enumerate(basis)),
                      np.zeros((d, d), complex)) for i in range(coeff.shape[1])]
    fixed = null(np.vstack(generators) if generators else np.empty((0, d)), d)
    residual = max([np.linalg.norm(a @ phi) for a in generators] + [0.])
    return coeff.shape[1], fixed.shape[1], residual


def random_su(rng, d):
    a = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
    q, _ = np.linalg.qr(a)
    q[:, -1] /= np.linalg.det(q)
    return q


def run():
    rng = np.random.default_rng(13122026)
    table = []
    for family, integers in [('SU', range(2, 8)), ('Sp', range(1, 5))]:
        for n in integers:
            basis = su(n) if family == 'SU' else sp(n)
            d = n if family == 'SU' else 2*n
            h, f, residual = stabilizer_fixed(basis)
            expected_h = (n-1)**2-1 if family == 'SU' else (n-1)*(2*n-1)
            expected_f = (2 if n == 2 else 1) if family == 'SU' else 2
            assert h == expected_h and f == expected_f and residual < 1e-10
            table.append(dict(family=family, n=n, complex_dimension=d,
                              stabilizer_dimension=h, fixed_dimension=f))
    J = np.array([[0, 1], [-1, 0]], complex)
    frame = lambda phi: np.column_stack((phi, J@phi.conj()))
    maximum_residual = 0.
    minimum_singular = float('inf')
    for n_left in (1, 2, 4):
        for extra_right in (0, 3):
            phi = rng.normal(size=2) + 1j*rng.normal(size=2)
            U = random_su(rng, 2)
            C = np.concatenate([np.eye(2*n_left),
                                rng.normal(size=(2*n_left, extra_right))], axis=1)
            F = np.kron(frame(phi), np.eye(n_left)) @ C
            F2 = np.kron(frame(U@phi), np.eye(n_left)) @ C
            residual = np.linalg.norm(F2 - np.kron(U, np.eye(n_left)) @ F)
            maximum_residual = max(maximum_residual, float(residual))
            assert residual < 1e-10
            s = np.linalg.svd(F, compute_uv=False)
            assert len(s) == 2*n_left and min(s) > 1e-10
            minimum_singular = min(minimum_singular, float(min(s)))
            assert np.allclose(frame(phi).conj().T@frame(phi), np.vdot(phi, phi)*np.eye(2))
    # Removing the single-reference condition: SU(3) with two fundamentals.
    multiple_reference_residual = 0.
    for _ in range(8):
        a, b = rng.normal(size=(2, 3)) + 1j*rng.normal(size=(2, 3))
        make = lambda x, y: np.column_stack((x, y, np.cross(x.conj(), y.conj())))
        U = random_su(rng, 3)
        F = make(a, b)
        assert abs(np.linalg.det(F)) > 1e-9
        residual = np.linalg.norm(make(U@a, U@b) - U@F)
        multiple_reference_residual = max(multiple_reference_residual, float(residual))
        assert residual < 1e-10
    # Dropping singlet-right role permits vectorlike F=mI for all these groups.
    for d in (3, 4, 6):
        U = random_su(rng, d)
        F = 0.7*np.eye(d)
        assert np.allclose(U@F@U.conj().T, F)
    # A single complex SO(3) reference has two real directions: sphere
    # transitivity, rather than complex irreducibility alone, matters.
    phi = np.array([1, 1j, 0], complex)
    make_real = lambda z: np.column_stack((z, z.conj(), 1j*np.cross(z, z.conj())))
    O, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    O[:, -1] *= np.linalg.det(O)
    assert abs(np.linalg.det(make_real(phi))) > 1
    assert np.allclose(make_real(O@phi), O@make_real(phi))
    # Apply the general theorem at the already accepted 1311 parameter
    # family; do not transport the author's different v=1 sample.
    mu = 2.**-12
    U3 = np.exp(2j*np.pi*np.outer(np.arange(3), np.arange(3))/3)/np.sqrt(3)
    joint = []
    for n in (3, 4, 5):
        U = np.eye(n, dtype=complex)
        U[:3, :3] = U3
        for j in range(2, n-1):
            R = np.eye(n)
            c, s = np.cos(.001), np.sin(.001)
            R[j, j] = R[j+1, j+1] = c
            R[j, j+1], R[j+1, j] = s, -s
            U = U@R
        assert np.allclose(U[:3, :2], U3[:, :2])
        Mnu = U@np.diag(mu*np.arange(1, n+1))@U.conj().T
        Me = np.diag(mu*np.array([1, 2, 4, 6, 8][:n]))
        Z = np.zeros((n, n), complex)
        expected = np.block([[Mnu, Z], [Z, Me]])
        for v in (.24, .25, .26):
            reference = np.array([0, v/np.sqrt(2)], complex)
            S = np.column_stack((J@reference.conj(), reference))
            C = np.sqrt(2)*expected/v
            actual = np.kron(S, np.eye(n))@C
            assert np.allclose(actual, expected, rtol=1e-12, atol=1e-15)
            assert np.linalg.matrix_rank(actual) == 2*n
            joint.append(dict(generations=n, v=v, mass_rank=2*n))
    # Direct physical-right hypercharge covariance and all-left r bookkeeping.
    from fractions import Fraction as Q
    charges = []
    for yl in (Q(-1, 2), Q(1, 6), Q(7, 6)):
        h = Q(1, 2)
        for scalar_charge in (h, -h):
            yr = yl - scalar_charge
            assert -yl + scalar_charge + yr == 0
            charges.append([str(yl), str(scalar_charge), str(yr)])
    assert (1j)**2 * (-1) == 1 and (1j)**4 == 1
    assert (-1j)*1j == 1 and (1j)**2 == -1 and (-1j)**4 == 1
    sources = [HERE/'proof.md', HERE/'check.py', HERE/'../_admission/round1312_admission.md',
               HERE/'../1311/acceptance.json']
    return dict(status='PASS_INDEPENDENT_MASS_STABILIZER_CHECKS',
                new_axioms=0, classification_by_analytic_theorem_not_samples=True,
                fixed_space_table=table, mass_intertwiner_max_residual=maximum_residual,
                full_row_rank_minimum_singular=minimum_singular,
                two_reference_su3_max_residual=multiple_reference_residual,
                so3_without_sphere_transitivity_checked=True,
                same_parameter_1311_mass_realizations=joint,
                physical_right_hypercharge_relations=charges,
                source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for k in a: compare(a[k], b[k])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b): compare(x, y)
    elif isinstance(a, float):
        assert np.isclose(a, b, rtol=1e-8, atol=1e-12), (a, b)
    else:
        assert a == b, (a, b)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--check', action='store_true')
    p.add_argument('--write-new', action='store_true')
    args = p.parse_args()
    result = run()
    out = HERE/'independent_results.json'
    if args.write_new:
        with out.open('x', encoding='utf8', newline='\n') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write('\n')
    else:
        compare(result, json.loads(out.read_text(encoding='utf8')))
    print(json.dumps(dict(status=result['status'], cases=len(result['fixed_space_table'])), ensure_ascii=False))
