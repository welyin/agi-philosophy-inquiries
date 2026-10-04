"""Round 407: exact local spectral-rigidity certificates, with explicit scope.

This uses a *given* open-chain two-body class.  No geometry is generated.
The finite-field minors are exact certificates over Gaussian integers; floating
calculations only check derivatives, an isospectral curve, and readout ambiguity.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np

TARGET = Path(__file__).with_name("spectral_tps_audit_results.json")
PRIME = 1009
IMAGINARY_ROOT = 469
PAULI = [np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]], complex), np.diag([1., -1.])]
LETTERS = "IXYZ"


def chain_strings(n):
    strings = [(0,)*n]
    for i in range(n):
        for a in range(1, 4):
            s = [0]*n
            s[i] = a
            strings.append(tuple(s))
    for i in range(n-1):
        for a, b in itertools.product(range(1, 4), repeat=2):
            s = [0]*n
            s[i:i+2] = [a, b]
            strings.append(tuple(s))
    return strings


def tensor_basis(strings, modular=False):
    if modular:
        mats = [np.eye(2, dtype=np.int64), np.array([[0, 1], [1, 0]]),
                np.array([[0, -IMAGINARY_ROOT], [IMAGINARY_ROOT, 0]]),
                np.diag([1, -1])]
    else:
        mats = PAULI
    out = []
    for s in strings:
        m = np.ones((1, 1), dtype=np.int64 if modular else complex)
        for a in s:
            m = np.kron(m, mats[a])
            if modular:
                m %= PRIME
        out.append(m)
    return np.array(out)


def determinant_mod(matrix):
    a = np.array(matrix, dtype=np.int64) % PRIME
    assert a.shape[0] == a.shape[1]
    det = 1
    for j in range(len(a)):
        rows = np.flatnonzero(a[j:, j])
        if not len(rows):
            return 0
        k = j+int(rows[0])
        if k != j:
            a[[j, k]] = a[[k, j]]
            det = -det
        pivot = int(a[j, j])
        det = det*pivot % PRIME
        inverse = pow(pivot, -1, PRIME)
        for k in range(j+1, len(a)):
            factor = int(a[k, j])*inverse % PRIME
            a[k] = (a[k]-factor*a[j]) % PRIME
    return int(det % PRIME)


def rank_certificate(matrix):
    original = np.array(matrix, dtype=np.int64) % PRIME
    a = original.copy()
    row_ids = list(range(a.shape[0]))
    selected_rows, selected_columns = [], []
    r = 0
    for j in range(a.shape[1]):
        rows = np.flatnonzero(a[r:, j])
        if not len(rows):
            continue
        k = r+int(rows[0])
        a[[r, k]] = a[[k, r]]
        row_ids[r], row_ids[k] = row_ids[k], row_ids[r]
        selected_rows.append(row_ids[r])
        selected_columns.append(j)
        a[r] = a[r]*pow(int(a[r, j]), -1, PRIME) % PRIME
        for k in range(r+1, a.shape[0]):
            if a[k, j]:
                a[k] = (a[k]-a[k, j]*a[r]) % PRIME
        r += 1
        if r == a.shape[0]:
            break
    minor = original[np.ix_(selected_rows, selected_columns)]
    value = determinant_mod(minor)
    assert r == 0 or value != 0
    return dict(rank=r, rows=selected_rows, columns=selected_columns,
                minor_determinant_mod_prime=value)


def gauge_coefficients(strings, coefficients):
    """Exact coefficients of i[Q,H] for each one-site Pauli Q."""
    n = len(strings[0])
    index = {s: i for i, s in enumerate(strings)}
    c = np.zeros((len(strings), 3*n), dtype=np.int64)
    positive_cycles = {(1, 2), (2, 3), (3, 1)}
    for site in range(n):
        for a in range(1, 4):
            column = 3*site+a-1
            for s, h in zip(strings, coefficients):
                b = s[site]
                if not h or b in (0, a):
                    continue
                target = list(s)
                target[site] = 6-a-b
                sign = -2 if (a, b) in positive_cycles else 2
                c[index[tuple(target)], column] += sign*int(h)
    return c


def spectral_jacobian(values, vectors, basis):
    # values fixes eigenvalue ordering; simple spectrum is checked separately.
    assert len(values) == vectors.shape[1]
    return np.einsum("ki,akl,li->ia", vectors.conj(), basis, vectors).real


@lru_cache(None)
def model(n):
    strings = chain_strings(n)
    coeff = np.random.default_rng(20260923+n).integers(-5, 6, len(strings))
    coeff[0] = 0
    basis = tensor_basis(strings)
    basis_mod = tensor_basis(strings, True)
    h = np.einsum("a,aij->ij", coeff, basis)
    h_mod = np.einsum("a,aij->ij", coeff, basis_mod) % PRIME
    dim = 2**n
    power = np.eye(dim, dtype=np.int64)
    moments, jac = [], []
    for k in range(2*dim-1):
        moments.append(int(np.trace(power)) % PRIME)
        if k < dim:
            jac.append(np.einsum("ij,aji->a", power, basis_mod) % PRIME)
        power = power@h_mod % PRIME
    jac = np.array(jac)
    hankel = np.array([[moments[i+j] for j in range(dim)] for i in range(dim)])
    gauge = gauge_coefficients(strings, coeff)
    values, vectors = np.linalg.eigh(h)
    numeric_jac = spectral_jacobian(values, vectors, basis)
    singular = np.linalg.svd(numeric_jac, compute_uv=False)
    return dict(n=n, strings=strings, coefficients=coeff, basis=basis, h=h,
                jac=jac, hankel=hankel, gauge=gauge, values=values, vectors=vectors,
                numeric_jac=numeric_jac, singular_values=singular)


def unitary_from_h(h, t):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*t*values))@vectors.conj().T


def isospectral_curve():
    m = model(5)
    h, basis, coeff = m["h"], m["basis"], m["coefficients"]
    q_gauge = np.linalg.qr(m["gauge"].astype(float), mode="reduced")[0]
    stacked = np.vstack([m["numeric_jac"], q_gauge.T])
    _, _, vh = np.linalg.svd(stacked, full_matrices=True)
    free = vh[-5:]
    target = np.array([.02, 0., 0., 0., 0.])
    delta = free.T@target
    iterations = 0
    for iterations in range(1, 13):
        new_h = h+np.einsum("a,aij->ij", delta, basis)
        ev, vec = np.linalg.eigh(new_h)
        residual = np.concatenate([ev-m["values"], q_gauge.T@delta,
                                   free@delta-target])
        if np.max(np.abs(residual)) < 3e-13:
            break
        derivative = np.vstack([spectral_jacobian(ev, vec, basis), q_gauge.T, free])
        delta -= np.linalg.solve(derivative, residual)
    else:
        raise AssertionError("Local isospectral continuation did not converge")
    new_h = h+np.einsum("a,aij->ij", delta, basis)
    ev, vec = np.linalg.eigh(new_h)
    intertwiner = vec@m["vectors"].conj().T
    # Entire operator comparison also controls every unacted-on reference.
    errors = [np.linalg.norm(unitary_from_h(new_h, t)@intertwiner
                            -intertwiner@unitary_from_h(h, t), 2)
              for t in (.1, .7, 2.)]
    return dict(step=float(target[0]), iterations=iterations,
                coefficient_delta=delta.tolist(),
                coefficient_delta_norm=float(np.linalg.norm(delta)),
                local_gauge_slice_residual=float(np.linalg.norm(q_gauge.T@delta)),
                free_coordinate_residual=float(np.linalg.norm(free@delta-target)),
                max_spectrum_error=float(np.max(np.abs(ev-m["values"]))),
                full_operator_intertwining_errors=errors,
                source_coefficients=coeff.tolist())


@lru_cache(None)
def report():
    rows = []
    for n in (5, 6):
        m = model(n)
        jcert = rank_certificate(m["jac"])
        gcert = rank_certificate(m["gauge"])
        hcert = rank_certificate(m["hankel"])
        direct_error = 0.
        for col in range(3*n):
            q = m["basis"][1+col]
            direct = 1j*(q@m["h"]-m["h"]@q)
            via_coeff = np.einsum("a,aij->ij", m["gauge"][:, col], m["basis"])
            direct_error = max(direct_error, float(np.max(np.abs(direct-via_coeff))))
        direction = np.random.default_rng(770+n).normal(size=len(m["strings"]))
        direction /= np.linalg.norm(direction)
        perturbation = np.einsum("a,aij->ij", direction, m["basis"])
        eps = 1e-5
        finite_derivative = (np.linalg.eigvalsh(m["h"]+eps*perturbation)
                             -np.linalg.eigvalsh(m["h"]-eps*perturbation))/(2*eps)
        rank = jcert["rank"]
        rows.append(dict(
            qubits=n, hilbert_dimension=2**n, parameter_dimension=len(m["strings"]),
            pauli_strings=["".join(LETTERS[a] for a in s) for s in m["strings"]],
            coefficients=m["coefficients"].tolist(),
            moment_jacobian_certificate=jcert, gauge_certificate=gcert,
            hankel_certificate=hcert,
            modular_gauge_annihilation=bool(np.all(m["jac"]@m["gauge"] % PRIME == 0)),
            direct_gauge_formula_error=direct_error,
            minimum_eigenvalue_gap=float(np.min(np.diff(m["values"]))),
            numerical_spectral_jacobian_rank=int(np.sum(m["singular_values"] > 1e-9)),
            smallest_retained_singular_value=float(m["singular_values"][rank-1]),
            largest_discarded_singular_value=float(m["singular_values"][rank]) if rank < len(m["singular_values"]) else None,
            derivative_check_error=float(np.max(np.abs(finite_derivative-m["numeric_jac"]@direction))),
            regular_isospectral_level_dimension=len(m["strings"])-rank,
            local_unitary_orbit_dimension=gcert["rank"],
            infinitesimal_quotient_dimension=len(m["strings"])-rank-gcert["rank"]))

    m = model(6)
    t = .13
    w = unitary_from_h(m["h"], t)
    x0 = m["basis"][1]
    z1 = m["basis"][6]
    moved = w@x0@w.conj().T
    singles = m["basis"][:19]
    expansion = np.einsum("aij,ji->a", singles, moved)/64
    residual = moved-np.einsum("a,aij->ij", expansion, singles)
    return dict(
        date="2026-09-23", round=407, scientific_base_through_round=406,
        finite_field_prime=PRIME, image_of_imaginary_unit=IMAGINARY_ROOT,
        gaussian_integer_embedding_valid=IMAGINARY_ROOT**2 % PRIME == PRIME-1,
        locality_class="all on-site and nearest-neighbor Pauli terms on an open chain, plus identity",
        exact_certificates=rows, five_qubit_isospectral_curve=isospectral_curve(),
        commuting_dictionary_example=dict(
            time=t, h_commutator_norm=float(np.linalg.norm(w@m["h"]-m["h"]@w, 2)),
            moved_readout_commutator_with_original_neighbor=float(np.linalg.norm(moved@z1-z1@moved, 2)),
            moved_readout_distance_from_one_body_span=float(np.linalg.norm(residual, "fro")/8)),
        exact_integer_certificates=True, six_qubit_local_rigidity_mod_symmetry=True,
        five_qubit_nontrivial_local_isospectral_family=True,
        generic_claim_only_within_specified_locality_class=True,
        global_unique_tps_proved=False, all_graph_classes_excluded=False,
        actual_access_algebras_identified=False, spatial_dimension_generated=False,
        full_cognition_to_gr_refuted=False, infinite_resource_choice_required=False)


class Audit(unittest.TestCase):
    def test_gaussian_integer_reduction_and_pauli_basis(self):
        self.assertTrue(report()["gaussian_integer_embedding_valid"])
        for n in (5, 6):
            m = model(n)
            gram = np.einsum("aij,bji->ab", m["basis"], m["basis"])/(2**n)
            self.assertLess(np.max(np.abs(gram-np.eye(len(gram)))), 1e-12)
            self.assertEqual(len(gram), 12*n-8)

    def test_exact_simple_spectrum_certificates(self):
        for row in report()["exact_certificates"]:
            self.assertEqual(row["hankel_certificate"]["rank"], row["hilbert_dimension"])
            self.assertNotEqual(row["hankel_certificate"]["minor_determinant_mod_prime"], 0)
            self.assertGreater(row["minimum_eigenvalue_gap"], 1e-5)

    def test_gauge_directions_are_exact_and_independent(self):
        for row in report()["exact_certificates"]:
            self.assertEqual(row["gauge_certificate"]["rank"], 3*row["qubits"])
            self.assertTrue(row["modular_gauge_annihilation"])
            self.assertEqual(row["direct_gauge_formula_error"], 0)

    def test_exact_jacobian_ranks_and_regular_level_dimensions(self):
        rows = report()["exact_certificates"]
        self.assertEqual([r["moment_jacobian_certificate"]["rank"] for r in rows], [32, 46])
        self.assertEqual([r["infinitesimal_quotient_dimension"] for r in rows], [5, 0])
        for row in rows:
            self.assertNotEqual(row["moment_jacobian_certificate"]["minor_determinant_mod_prime"], 0)
            self.assertEqual(row["numerical_spectral_jacobian_rank"],
                             row["moment_jacobian_certificate"]["rank"])

    def test_spectral_derivative_against_independent_finite_difference(self):
        for row in report()["exact_certificates"]:
            self.assertLess(row["derivative_check_error"], 2e-8)

    def test_nonlinear_isospectral_curve_and_full_operator_intertwining(self):
        row = report()["five_qubit_isospectral_curve"]
        self.assertGreaterEqual(row["coefficient_delta_norm"], .0199)
        for key in ("local_gauge_slice_residual", "free_coordinate_residual", "max_spectrum_error"):
            self.assertLess(row[key], 1e-11)
        self.assertLess(max(row["full_operator_intertwining_errors"]), 1e-10)

    def test_spectral_rigidity_does_not_fix_actual_readout_algebras(self):
        row = report()["commuting_dictionary_example"]
        self.assertLess(row["h_commutator_norm"], 1e-11)
        self.assertGreater(row["moved_readout_commutator_with_original_neighbor"], .1)
        self.assertGreater(row["moved_readout_distance_from_one_body_span"], .1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    result = dict(report())
    result["checks"] = dict(run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors))
    result["runtime"] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))

