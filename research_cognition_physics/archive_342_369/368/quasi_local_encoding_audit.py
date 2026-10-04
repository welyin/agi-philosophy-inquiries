"""Round 368: static local dictionaries cannot give onsite motion a growing cone."""
import itertools
import unittest
import numpy as np
from growing_stream_audit import main

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
HAD = (X + Z) / np.sqrt(2)
P = (I, X, Y, Z)


def kron(parts):
    out = np.ones((1, 1), dtype=complex)
    for item in parts:
        out = np.kron(out, item)
    return out


def at(n, site, op):
    return kron([op if j == site else I for j in range(n)])


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def comm(a, b):
    return a @ b - b @ a


def bit_table(n):
    return ((np.arange(2 ** n)[:, None] >> np.arange(n - 1, -1, -1)) & 1)


def cluster_dictionary(n):
    bits = bit_table(n)
    parity = np.sum(bits[:, :-1] * bits[:, 1:], axis=1)
    phases = (-1.0) ** parity
    return kron([HAD] * n) * phases


def source_energies(n, frequencies=None):
    return bit_table(n) @ (np.ones(n) if frequencies is None else frequencies)


def source_phases(n, time, frequencies=None):
    frequencies = np.ones(n) if frequencies is None else np.asarray(frequencies)
    local = np.exp(-1j * time * frequencies)
    # Product evaluation preserves the independent-site identity at long times.
    return np.prod(np.where(bit_table(n), local[None, :], 1), axis=1)


def encoded_unitary(w, time, frequencies=None):
    n = int(round(np.log2(len(w))))
    return w.conj().T @ (source_phases(n, time, frequencies)[:, None] * w)


def alpha(w, observable, time, frequencies=None):
    u = encoded_unitary(w, time, frequencies)
    return u.conj().T @ observable @ u


def cluster_terms(n):
    terms = []
    for j in range(n):
        factors = [I] * n
        factors[j] = X
        if j:
            factors[j - 1] = Z
        if j + 1 < n:
            factors[j + 1] = Z
        terms.append(kron(factors))
    return terms


def density(v):
    return np.outer(v, v.conj())


def reduced(rho, dims, keep):
    drop = [j for j in range(len(dims)) if j not in keep]
    order = list(keep) + drop
    da = int(np.prod([dims[j] for j in keep]))
    db = int(np.prod([dims[j] for j in drop]))
    view = rho.reshape(tuple(dims) * 2).transpose(
        order + [j + len(dims) for j in order]).reshape(da, db, da, db)
    return np.einsum("abcb->ac", view)


def distance(a, b):
    return float(np.abs(np.linalg.eigvalsh(a - b)).sum() / 2)


def random_density(d, seed):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    a = a @ a.conj().T
    return a / np.trace(a)


def near_signal():
    w = cluster_dictionary(2)
    receiver = np.array([1, 1j]) / np.sqrt(2)
    u = encoded_unitary(w, np.pi / 2)
    out = [reduced(density(u @ np.kron(np.eye(2)[:, j], receiver)),
                   [2, 2], [1]) for j in (0, 1)]
    return {"receiver_distance": distance(*out),
            "receiver_Z": [float(np.trace(z @ Z).real) for z in out]}


def pair_entanglement():
    u = encoded_unitary(cluster_dictionary(2), np.pi / 2)
    initial = np.eye(4)[:, 0]
    final = u @ initial
    rho = reduced(density(final), [2, 2], [0])
    return float(np.trace(rho @ rho).real)


def random_circuit(n, depth, seed=368):
    rng = np.random.default_rng(seed)
    w = np.eye(2 ** n, dtype=complex)
    for layer in range(depth):
        for site in range(layer % 2, n - 1, 2):
            raw = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
            gate, _ = np.linalg.qr(raw)
            full = kron([np.eye(2 ** site), gate, np.eye(2 ** (n - site - 2))])
            w = full @ w
    return w


def repetition_embedding(cells):
    v = np.zeros((4 ** cells, 2 ** cells), dtype=complex)
    for col, b in enumerate(bit_table(cells)):
        row = 0
        for value in b:
            row = 4 * row + 3 * value
        v[row, col] = 1
    return v


def logical_representative(cells, site, label):
    # pi(A)=tilde(W A W^dagger), with tilde X=XX, Z=ZI, Y=YX.
    local = [np.eye(4, dtype=complex) for _ in range(cells)]
    tx, ty, tz = np.kron(X, X), np.kron(Y, X), np.kron(Z, I)
    if label == 0:
        return kron(local)
    if label == 3:
        local[site] = tx
        return kron(local)
    local[site] = tz if label == 1 else -ty
    for j in (site - 1, site + 1):
        if 0 <= j < cells:
            local[j] = tx
    return kron(local)


def invariant_code(cells):
    w = cluster_dictionary(cells)
    v = repetition_embedding(cells) @ w
    energy = source_energies(2 * cells)
    h = w.conj().T @ (2 * source_energies(cells)[:, None] * w)
    return v, energy, h


def quasi_dictionary(n, theta=.3, q=.5, radius=None):
    z = 1 - 2 * bit_table(n)
    energy = np.zeros(2 ** n)
    for i in range(n):
        for j in range(i + 1, n):
            if radius is None or j - i <= radius:
                energy += theta * q ** (j - i) * z[:, i] * z[:, j]
    return kron([HAD] * n) * np.exp(-1j * energy)


def tail_bound(radius, theta=.3, q=.5):
    return 4 * abs(theta) * q ** (radius + 1) / (1 - q)


def report():
    n, site = 7, 3
    w = cluster_dictionary(n)
    rows = []
    for time in (.2, .9, np.pi / 2, 37.25):
        evolved = alpha(w, at(n, site, X), time)
        rows.append({"time": time,
                     "distance_2_commutator": opnorm(comm(evolved, at(n, 5, X))),
                     "distance_2_analytic": float(2 * abs(np.sin(time))),
                     "distance_3_commutator": opnorm(comm(evolved, at(n, 6, X)))})
    nq = 6
    qw = quasi_dictionary(nq)
    quasi_rows = []
    a, b = at(nq, 0, X), at(nq, 5, X)
    bound = min(2., 2 * (tail_bound(2) + tail_bound(2)))
    for time in (.2, .9, 7.3, 41.):
        quasi_rows.append({"time": time, "distance": 5,
                           "commutator": opnorm(comm(alpha(qw, a, time), b)),
                           "time_independent_upper": bound})
    v, energies, h = invariant_code(3)
    errors = []
    for j in range(3):
        for label in (1, 2, 3):
            pi = logical_representative(3, j, label)
            errors.append(np.linalg.norm(pi @ v - v @ at(3, j, P[label])))
    return {
        "round": 368,
        "scope": {
            "proved": "Static finite-range unitary dictionaries, or invariant isometric codes "
                      "with full local intertwining representatives, cannot turn onsite natural "
                      "motion into propagation beyond a fixed 2R buffer. Quasi-local "
                      "dictionaries give a time-uniform commutator tail.",
            "state_quantifier": "Exact operator identities hold for every state and passive "
                                "reference, including correlated inputs. Numerical states are witnesses.",
            "extra_inputs": ["specified site metric and tensor factors",
                             "onsite natural generator",
                             "static dictionary with declared range or approximation tails",
                             "code invariance and complete local representatives in the isometric case"],
            "not_claimed": ["all cognitive models have these dictionaries",
                            "all globally encoded dynamics have zero influence",
                            "generic interacting many-body localization has a time-independent cone",
                            "no controlled relay protocol can transmit information",
                            "GR or its full cognitive program is disproved"],
        },
        "cluster": {"sites": n, "dictionary_range": 1, "uniform_support_radius": 2,
                    "near_perfect_signal": near_signal(),
                    "two_site_product_input_final_purity": pair_entanglement(),
                    "commutators": rows},
        "quasi_local": {"dictionary": "Hadamard times exp(-i sum theta*q^distance Zi Zj)",
                        "theta": .3, "q": .5, "single_site_tail": "4*theta*q^(r+1)/(1-q)",
                        "commutators": quasi_rows},
        "isometric_code": {"logical_cells": 3, "physical_qubits": 6,
                           "invariance_error": float(np.linalg.norm(energies[:, None] * v - v @ h)),
                           "max_full_local_intertwiner_error": float(max(errors))},
        "continuum_scope": "For fixed positive macroscopic separation, a_N R_N -> 0 gives "
                           "exact zero response eventually, for every time rescaling. "
                           "Quasi-local tails require epsilon_XN(rN)+epsilon_YN(sN)->0, "
                           "including support-size prefactors.",
        "sources": ["https://arxiv.org/pdf/1409.1252",
                    "https://arxiv.org/pdf/quant-ph/0004051",
                    "https://arxiv.org/html/1701.05182"],
    }


class Checks(unittest.TestCase):
    def test_01_cluster_is_a_local_dictionary_of_onsite_energy(self):
        for n in (2, 4, 5):
            w = cluster_dictionary(n)
            np.testing.assert_allclose(w.conj().T @ w, np.eye(2 ** n), atol=2e-15)
            h = w.conj().T @ (source_energies(n)[:, None] * w)
            explicit = sum((np.eye(2 ** n) - k) / 2 for k in cluster_terms(n))
            np.testing.assert_allclose(h, explicit, atol=3e-15)

    def test_02_cluster_terms_commute_and_direct_exponential_agrees(self):
        n = 4
        terms = cluster_terms(n)
        for a, b in itertools.combinations(terms, 2):
            np.testing.assert_allclose(a @ b, b @ a)
        h = sum((np.eye(2 ** n) - k) / 2 for k in terms)
        values, vectors = np.linalg.eigh(h)
        direct = (vectors * np.exp(-.73j * values)) @ vectors.conj().T
        np.testing.assert_allclose(direct, encoded_unitary(cluster_dictionary(n), .73), atol=4e-15)

    def test_03_both_dictionary_directions_have_radius_one(self):
        n = 5
        w = cluster_dictionary(n)
        for j, label in itertools.product(range(n), (X, Y, Z)):
            a = at(n, j, label)
            for mapped in (w @ a @ w.conj().T, w.conj().T @ a @ w):
                for k in range(n):
                    if abs(j - k) > 1:
                        for b in (X, Y, Z):
                            self.assertLess(np.linalg.norm(comm(mapped, at(n, k, b))), 2e-14)

    def test_04_no_operator_tail_beyond_radius_two_at_all_sampled_times(self):
        n, w = 7, cluster_dictionary(7)
        for time in (.1, .9, np.pi / 2, 37.25, 1000.17):
            for label in (X, Y, Z):
                evolved = alpha(w, at(n, 3, label), time)
                for k in (0, 6):
                    for b in (X, Y, Z):
                        self.assertLess(np.linalg.norm(comm(evolved, at(n, k, b))), 8e-14)

    def test_05_radius_two_is_tight_not_zero_dynamics(self):
        n, w = 5, cluster_dictionary(5)
        for time in (.2, .9, np.pi / 2):
            evolved = alpha(w, at(n, 2, X), time)
            self.assertAlmostEqual(opnorm(comm(evolved, at(n, 4, X))),
                                   2 * abs(np.sin(time)), places=13)

    def test_06_neighbor_signal_is_perfect_for_an_explicit_preparation(self):
        row = near_signal()
        self.assertAlmostEqual(row["receiver_distance"], 1, places=13)
        np.testing.assert_allclose(row["receiver_Z"], [-1, 1], atol=3e-15)

    def test_07_cluster_dynamics_really_generates_logical_entanglement(self):
        self.assertAlmostEqual(pair_entanglement(), .5, places=13)

    def test_08_nonclifford_finite_depth_and_incommensurate_onsite_frequencies(self):
        n, depth = 6, 2
        w = random_circuit(n, depth)
        frequencies = np.sqrt(np.arange(1, n + 1))
        for time in (.3, 2.7, 51.1):
            evolved = alpha(w, at(n, 0, X), time, frequencies)
            for b in (X, Y, Z):
                self.assertLess(np.linalg.norm(comm(evolved, at(n, 5, b))), 5e-14)

    def test_09_invariant_code_and_complete_representatives(self):
        cells = 3
        v, energies, h = invariant_code(cells)
        np.testing.assert_allclose(v.conj().T @ v, np.eye(2 ** cells), atol=2e-15)
        np.testing.assert_allclose(energies[:, None] * v, v @ h, atol=3e-15)
        for j, label in itertools.product(range(cells), range(4)):
            pi = logical_representative(cells, j, label)
            np.testing.assert_allclose(pi @ v, v @ at(cells, j, P[label]), atol=2e-15)
            np.testing.assert_allclose(pi.conj().T @ v,
                                       v @ at(cells, j, P[label]).conj().T, atol=2e-15)

    def test_10_code_evolution_and_reference(self):
        v, energies, h = invariant_code(3)
        rng = np.random.default_rng(1368)
        psi = rng.normal(size=(8, 3)) + 1j * rng.normal(size=(8, 3))
        psi /= np.linalg.norm(psi)
        values, vectors = np.linalg.eigh(h)
        for time in (.2, 1.4):
            exact = np.exp(-1j * time * energies)[:, None] * (v @ psi)
            desired = v @ (vectors * np.exp(-1j * time * values)) @ vectors.conj().T @ psi
            np.testing.assert_allclose(exact, desired, atol=4e-15)

    def test_11_invariant_code_far_logical_algebras_do_not_signal(self):
        cells = 4
        v, energies, _ = invariant_code(cells)
        pa = logical_representative(cells, 0, 1)
        pb = logical_representative(cells, 3, 3)
        for time in (.2, 1.1, 9.3):
            phases = source_phases(2 * cells, time)
            evolved = phases.conj()[:, None] * pb * phases[None, :]
            np.testing.assert_allclose(comm(evolved, pa), np.zeros_like(pa), atol=2e-15)
            a = v.conj().T @ pa @ v
            b = v.conj().T @ evolved @ v
            self.assertLess(np.linalg.norm(comm(a, b)), 2e-14)

    def test_12_plain_compression_does_not_preserve_measurement_products(self):
        v = repetition_embedding(1)
        b = np.kron(X, I)
        compressed = v.conj().T @ b @ v
        np.testing.assert_allclose(compressed, np.zeros((2, 2)))
        np.testing.assert_allclose(v.conj().T @ b @ b @ v, I)
        self.assertEqual(opnorm(b @ v - v @ compressed), 1)

    def test_13_quasi_local_dictionary_truncation_bound(self):
        n, w = 6, quasi_dictionary(6)
        for j, radius in itertools.product((0, 2, 5), (0, 1, 2)):
            wr = quasi_dictionary(n, radius=radius)
            a = at(n, j, X)
            exact, approx = w @ a @ w.conj().T, wr @ a @ wr.conj().T
            self.assertLessEqual(opnorm(exact - approx), tail_bound(radius) + 2e-14)
            self.assertAlmostEqual(opnorm(approx), 1, places=13)
            for k in range(n):
                if abs(k - j) > radius:
                    self.assertLess(np.linalg.norm(comm(approx, at(n, k, X))), 2e-14)

    def test_14_quasi_local_commutator_bound_is_uniform_in_time(self):
        n, w = 6, quasi_dictionary(6)
        bound = 2 * (tail_bound(2) + tail_bound(2))
        for time in (.2, .9, 7.3, 41., 1000.17):
            actual = opnorm(comm(alpha(w, at(n, 0, X), time), at(n, 5, X)))
            self.assertLessEqual(actual, bound + 2e-14)

    def test_15_quasi_local_tail_need_not_be_exactly_zero(self):
        n, w = 6, quasi_dictionary(6)
        value = opnorm(comm(alpha(w, at(n, 0, X), .9), at(n, 5, X)))
        self.assertGreater(value, 1e-4)
        self.assertLess(value, .2)

    def test_16_arbitrary_local_cptp_and_correlated_reference_cannot_signal_far(self):
        n, reference = 5, 2
        rho = random_density(2 ** n * reference, 2368)
        eta = .41
        kraus = [np.diag([1, np.sqrt(1 - eta)]), np.array([[0, np.sqrt(eta)], [0, 0]])]
        acted = np.zeros_like(rho)
        for k in kraus:
            full = np.kron(at(n, 0, k), I)
            acted += full @ rho @ full.conj().T
        u = np.kron(encoded_unitary(cluster_dictionary(n), .93), I)
        first = reduced(u @ rho @ u.conj().T, [2] * n + [reference], [4, 5])
        second = reduced(u @ acted @ u.conj().T, [2] * n + [reference], [4, 5])
        self.assertLess(distance(first, second), 2e-14)


if __name__ == "__main__":
    main(__name__, "quasi_local_encoding_audit", report)
