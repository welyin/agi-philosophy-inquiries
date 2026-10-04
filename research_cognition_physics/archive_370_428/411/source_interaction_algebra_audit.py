"""Round 411: direct interaction algebra from one specified source factor.

Only Python/NumPy; exact finite-field certificates supplement analytic proofs.
No physical position, native environmental factorization, or new axiom is inferred.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / "source_interaction_algebra_audit_results.json"
PRIME = 1009
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
PAULI = (I, X, Y, Z)


def kron(*xs):
    out = np.ones((1, 1), dtype=complex)
    for x in xs:
        out = np.kron(out, x)
    return out


def comm(a, b):
    return a @ b - b @ a


@lru_cache(None)
def labels(n):
    return list(itertools.product(range(4), repeat=n))


@lru_cache(None)
def pauli_rules(n):
    basis = labels(n)
    index = {s: j for j, s in enumerate(basis)}
    out = {}
    for i, a in enumerate(basis):
        for j, b in enumerate(basis):
            phase, word = 0, []
            for u, v in zip(a, b):
                if u == 0 or v == 0:
                    word.append(u or v)
                elif u == v:
                    word.append(0)
                else:
                    word.append(6-u-v)
                    phase += 1 if (u, v) in ((1, 2), (2, 3), (3, 1)) else 3
            out[i, j] = (index[tuple(word)], phase % 4)
    return out


def coeff(n, terms):
    index = {s: j for j, s in enumerate(labels(n))}
    v = [0] * (4**n)
    for word, value in terms:
        v[index[tuple(word)]] += value
    return v


def product(a, b, n, kind):
    """Integer coefficients of [A,B]/(2i) or {A,B}/2."""
    out = [0] * len(a)
    rules = pauli_rules(n)
    for i, u in enumerate(a):
        if not u:
            continue
        for j, v in enumerate(b):
            if not v:
                continue
            k, phase = rules[i, j]
            if kind == "lie" and phase % 2:
                out[k] += u*v*(1 if phase == 1 else -1)
            elif kind == "jordan" and phase % 2 == 0:
                out[k] += u*v*(1 if phase == 0 else -1)
    return out


def add_modular_basis(v, basis):
    w = [int(x) % PRIME for x in v]
    for pivot in sorted(basis):
        factor = w[pivot]
        if factor:
            w = [(x-factor*y) % PRIME for x, y in zip(w, basis[pivot])]
    pivot = next((i for i, x in enumerate(w) if x), None)
    if pivot is None:
        return False
    inv = pow(w[pivot], -1, PRIME)
    basis[pivot] = [(x*inv) % PRIME for x in w]
    return True


def exact_closure(seeds, n, associative=False):
    """Independent integer words; modular nonzero minors certify real rank."""
    rows, tree, basis = [], [], {}
    initial = ([coeff(n, [((0,)*n, 1)])] if associative else []) + seeds
    for s, row in enumerate(initial):
        if add_modular_basis(row, basis):
            rows.append(row)
            tree.append(dict(seed=s))
    cursor = 0
    target = 4**n if associative else 4**n-1
    while cursor < len(rows) and len(rows) < target:
        for s, seed in enumerate(seeds):
            for kind in (("lie", "jordan") if associative else ("lie",)):
                row = product(seed, rows[cursor], n, kind)
                if add_modular_basis(row, basis):
                    rows.append(row)
                    tree.append(dict(parent=cursor, generator=s, product=kind))
        cursor += 1
    # Independently replay the certificate using its integer words.
    replay = []
    for node in tree:
        replay.append(initial[node["seed"]] if "seed" in node else
                      product(seeds[node["generator"]], replay[node["parent"]], n, node["product"]))
    assert replay == rows
    return rows, tree


def rank_mod(rows):
    basis = {}
    for row in rows:
        add_modular_basis(row, basis)
    return len(basis)


def source_seeds(p, q):
    return [coeff(2, [((a, 0), q), ((0, a), p)]) for a in (1, 2, 3)]


def algebra_case(p, q):
    seeds = source_seeds(p, q)
    rows, _ = exact_closure(seeds, 2, True)
    units = [coeff(2, [(word, 1)]) for word in labels(2)]
    equations = []
    for seed in seeds:
        columns = [product(seed, v, 2, "lie") for v in units]
        equations.extend(zip(*columns))
    cr = rank_mod(equations)
    return dict(epsilon_numerator=p, epsilon_denominator=q,
                algebra_dimension_lower_bound_mod_prime=len(rows),
                commutant_dimension_upper_bound_mod_prime=16-cr,
                analytic_algebra_dimension=4 if p == 0 else 10 if p == q else 16,
                analytic_commutant_dimension=4 if p == 0 else 2 if p == q else 1)


def hamiltonian(epsilon):
    return sum(kron(a, a, I)+kron(I, a, a)+epsilon*kron(a, I, a) for a in (X, Y, Z))


def extracted_coefficients(h):
    return [np.trace((kron(a, I, I)@h).reshape(2, 4, 2, 4), axis1=0, axis2=2)/2
            for a in PAULI]


def unitary(h, t):
    w, v = np.linalg.eigh(h)
    return (v*np.exp(-1j*t*w))@v.conj().T


def projector_span(basis):
    columns = np.column_stack([b.ravel() for b in basis])
    q, _ = np.linalg.qr(columns)
    return q @ q.conj().T


def covariance_report():
    cnot = np.zeros((4, 4), complex)
    for b in range(4):
        cnot[b ^ ((b >> 1) & 1), b] = 1
    w = kron(I, cnot)
    h = hamiltonian(0)
    transformed = w@h@w.conj().T
    b = extracted_coefficients(h)
    bp = extracted_coefficients(transformed)
    raw = [kron(a, I) for a in PAULI]
    transported = [cnot@r@cnot.conj().T for r in raw]
    rebuilt = [np.eye(4), *bp[1:]]
    errors = [np.linalg.norm(bp[j]-cnot@b[j]@cnot.conj().T) for j in range(4)]
    # Arbitrary common source-axis rotation does not change the generated span here.
    u = unitary(.31*X+.27*Y-.19*Z, 1)
    hu = kron(u, I, I)@h@kron(u.conj().T, I, I)
    bu = extracted_coefficients(hu)
    return dict(coefficient_transport_error=float(max(errors)),
                generated_algebra_transport_error=float(np.linalg.norm(projector_span(transported)-projector_span(rebuilt))),
                source_axis_independence_error=float(np.linalg.norm(projector_span(raw)-projector_span([np.eye(4), *bu[1:]]))),
                original_and_transported_algebra_dimensions=[4, 4],
                transported_coefficients=["XX", "YX", "ZI"])


def derivative_report():
    h = hamiltonian(.5)+.23*kron(Z, I, I)
    bs = extracted_coefficients(h)
    reconstruction = sum(kron(a, b) for a, b in zip(PAULI, bs))
    source_ops = [kron(a, I, I) for a in (X, Y, Z)]
    from_derivative = 1j/8*sum(comm(a, 1j*comm(h, a)) for a in source_ops)
    nonidentity = h-kron(I, bs[0])
    # Check the double-commutator condition on a genuine silent multiplicity at epsilon=0.
    h0 = hamiltonian(0)
    silent = kron(I, I, X)
    first = max(np.linalg.norm(comm(comm(h0, a), silent)) for a in source_ops)
    later = max(np.linalg.norm(comm(comm(h0, comm(h0, a)), silent)) for a in source_ops)
    return dict(coefficient_reconstruction_error=float(np.linalg.norm(reconstruction-h)),
                derivative_reconstruction_error=float(np.linalg.norm(from_derivative-nonidentity)),
                silent_factor_first_order_double_commutator=float(first),
                silent_factor_second_order_response=float(later))


def controllability_case(p, q):
    terms = []
    for a in (1, 2, 3):
        terms += [((a, a, 0), q), ((0, a, a), q), ((a, 0, a), p)]
    seeds = [coeff(3, terms), coeff(3, [((1, 0, 0), 1)]), coeff(3, [((3, 0, 0), 1)])]
    rows, tree = exact_closure(seeds, 3)
    assert all(row[0] == 0 for row in rows)
    paulis = [kron(*(PAULI[a] for a in word)) for word in labels(3)]
    matrices = [sum(c*b for c, b in zip(seed, paulis)) for seed in seeds]
    direct, direct_error = [], 0.
    for row, node in zip(rows, tree):
        actual = matrices[node["seed"]] if "seed" in node else comm(matrices[node["generator"]], direct[node["parent"]])/(2j)
        reconstructed = sum(c*b for c, b in zip(row, paulis))
        direct_error = max(direct_error, float(np.max(np.abs(actual-reconstructed))))
        direct.append(actual)
    # This is a 63 x 63 minor in the independent traceless Pauli basis.
    matrix = [[x % PRIME for x in row[1:]] for row in rows]
    determinant = 1
    if len(matrix) == 63:
        for j in range(63):
            k = next(k for k in range(j, 63) if matrix[k][j])
            if k != j:
                matrix[k], matrix[j] = matrix[j], matrix[k]
                determinant = -determinant
            pivot = matrix[j][j]
            determinant = determinant*pivot % PRIME
            inv = pow(pivot, -1, PRIME)
            for k in range(j+1, 63):
                factor = matrix[k][j]*inv % PRIME
                matrix[k] = [(a-factor*b) % PRIME for a, b in zip(matrix[k], matrix[j])]
    else:
        determinant = 0
    return dict(epsilon=[p, q], integer_hamiltonian_scale=q,
                lie_dimension_certified=len(rows), prime=PRIME,
                traceless_minor_determinant_mod_prime=determinant,
                independent_matrix_commutator_error=direct_error,
                certificate_tree=tree)


def protocol_kraus(epsilon):
    h = hamiltonian(epsilon)
    times = (.17, .29, .54)
    us = [unitary(h, t) for t in times]
    out = []
    for a in (0, 1):
        # Weak Z measurement, then a recorded, conditional basis choice.
        ka = kron(np.diag(np.sqrt((1+(-1)**a*.6*np.array([1, -1]))/2)), I, I)
        rot = unitary(.37*(X if a == 0 else Y), 1)
        for b in (0, 1):
            lb = kron(rot @ ((I+(-1)**b*X)/2), I, I)
            out.append(us[2]@lb@us[1]@ka@us[0])
    return out, sum(times)


def finite_window_report():
    base, total = protocol_kraus(0)
    rows = []
    for epsilon in (.25, .0625, .00390625, .000244140625):
        other, _ = protocol_kraus(epsilon)
        td = 0.
        for a, b in zip(base, other):
            # Unnormalized pure Choi blocks; exact rank-two trace-norm formula.
            va, vb = a.ravel()/np.sqrt(8), b.ravel()/np.sqrt(8)
            aa, bb = np.vdot(va, va).real, np.vdot(vb, vb).real
            overlap = abs(np.vdot(va, vb))**2
            td += .5*np.sqrt(max(0., (aa+bb)**2-4*overlap))
        normalization = np.linalg.norm(sum(k.conj().T@k for k in other)-np.eye(8))
        rows.append(dict(epsilon=epsilon, total_interaction_time=total,
                         hamiltonian_difference_norm=float(np.linalg.norm(hamiltonian(epsilon)-hamiltonian(0), 2)),
                         recorded_full_reference_choi_trace_distance=float(td),
                         all_protocols_half_diamond_upper_bound=min(1., 3*epsilon*total),
                         complete_instrument_error=float(normalization),
                         exact_direct_algebra_dimension=16))
    return rows


def arbitrary_receiver_report():
    rows = []
    for d in (3, 5):
        diag = np.diag(np.arange(d, dtype=float))
        path = np.eye(d, k=1)+np.eye(d, k=-1)
        projectors = []
        for j in range(d):
            p = np.eye(d)
            for k in range(d):
                if k != j:
                    p = p@(diag-k*np.eye(d))/(j-k)
            projectors.append(p)
        pe = max(np.linalg.norm(p-np.diag(np.eye(d)[j])) for j, p in enumerate(projectors))
        ee = 0.
        for j in range(d-1):
            e = np.zeros((d, d))
            e[j, j+1] = 1
            ee = max(ee, np.linalg.norm(projectors[j]@path@projectors[j+1]-e))
        rows.append(dict(receiver_hilbert_dimension=d, source_hilbert_dimension=2,
                         source_nonzero_traceless_coefficients=2,
                         spectral_projector_error=float(pe), adjacent_matrix_unit_error=float(ee),
                         analytically_generated_algebra_dimension=d*d))
    return rows


@lru_cache(None)
def report():
    return dict(round=411,
                scope="Conditional minimal direct interaction algebra from a specified finite source and Hamiltonian; no spatial neighborhood or dimension derived",
                derivative=derivative_report(), covariance=covariance_report(),
                algebra_cases=[algebra_case(p, q) for p, q in ((0, 1), (1, 2), (1, 1), (1, 16))],
                controllability=[controllability_case(0, 1), controllability_case(1, 2)],
                finite_window=finite_window_report(), arbitrary_receiver=arbitrary_receiver_report(),
                generated_algebra_equals_finite_cost_control_set=False,
                source_only_statistics_identify_global_hamiltonian=False,
                environmental_tensor_sites_given_in_general_theorem=False,
                source_factor_and_generator_still_inputs=True,
                unitary_covariance_counted_as_physical_nonuniqueness=False,
                all_time_silent_region_inferred_from_first_order=False,
                hilbert_dimension_identified_as_spatial_dimension=False,
                spatial_dimension_generated=False, full_cognition_to_gr_refuted=False,
                new_cognitive_axiom_adopted=False)


class Audit(unittest.TestCase):
    def test_first_derivative_and_first_order_scope(self):
        r = report()["derivative"]
        self.assertLess(r["coefficient_reconstruction_error"], 1e-13)
        self.assertLess(r["derivative_reconstruction_error"], 1e-13)
        self.assertEqual(r["silent_factor_first_order_double_commutator"], 0)
        self.assertGreater(r["silent_factor_second_order_response"], 1)

    def test_coordinate_covariance_preserves_canonical_algebra(self):
        r = report()["covariance"]
        for name in ("coefficient_transport_error", "generated_algebra_transport_error", "source_axis_independence_error"):
            self.assertLess(r[name], 1e-13)

    def test_exact_algebra_and_commutant_dimensions(self):
        for r in report()["algebra_cases"]:
            self.assertEqual(r["algebra_dimension_lower_bound_mod_prime"], r["analytic_algebra_dimension"])
            self.assertEqual(r["commutant_dimension_upper_bound_mod_prime"], r["analytic_commutant_dimension"])

    def test_chain_and_perturbed_model_are_fully_controllable(self):
        for r in report()["controllability"]:
            self.assertEqual(r["lie_dimension_certified"], 63)
            self.assertNotEqual(r["traceless_minor_determinant_mod_prime"], 0)
            self.assertEqual(r["independent_matrix_commutator_error"], 0)
            self.assertEqual(len(r["certificate_tree"]), 63)

    def test_small_change_in_complete_recorded_process_despite_algebra_jump(self):
        for r in report()["finite_window"]:
            self.assertAlmostEqual(r["hamiltonian_difference_norm"], 3*r["epsilon"], places=12)
            self.assertLessEqual(r["recorded_full_reference_choi_trace_distance"], r["all_protocols_half_diamond_upper_bound"]+1e-12)
            self.assertLess(r["complete_instrument_error"], 1e-13)

    def test_qubit_source_does_not_bound_receiver_dimension(self):
        for r in report()["arbitrary_receiver"]:
            self.assertLess(r["spectral_projector_error"], 1e-13)
            self.assertLess(r["adjacent_matrix_unit_error"], 1e-13)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    result = dict(report(), checks=dict(run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors)),
                  runtime=dict(python=platform.python_version(), numpy=np.__version__))
    if args.write_results:
        with TARGET.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k != "controllability"}, ensure_ascii=False, indent=2))
