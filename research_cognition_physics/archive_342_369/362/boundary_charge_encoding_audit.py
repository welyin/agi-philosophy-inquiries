"""Round 362: reversible boundary-charge encoding, dressed control, and memory cost.

The endpoint isometry changes the constraint representation. It is NOT a
fixed-Gauss-preserving dynamical admission protocol for arbitrary bare inputs.
"""
import itertools
import unittest
import numpy as np
from growing_stream_audit import main

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
PAULI = (I, X, Y, Z)
PLUS = np.ones(2, dtype=complex) / np.sqrt(2)


def bits(index, n):
    return np.array([(index >> (n - 1 - j)) & 1 for j in range(n)],
                    dtype=np.int64)


def integer(bit_string):
    value = 0
    for bit in bit_string:
        value = 2 * value + int(bit)
    return value


def kron_all(parts):
    result = np.array([1], dtype=complex)
    for part in parts:
        result = np.kron(result, part)
    return result


def density(ket):
    return np.outer(ket, ket.conj())


def random_density(d, seed=362):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    rho = a @ a.conj().T
    return rho / np.trace(rho)


def partial_trace(rho, dims, keep):
    keep = list(keep)
    drop = [i for i in range(len(dims)) if i not in keep]
    order = keep + drop
    k = int(np.prod([dims[i] for i in keep]))
    d = int(np.prod([dims[i] for i in drop]))
    view = rho.reshape(tuple(dims) * 2).transpose(
        order + [i + len(dims) for i in order]).reshape(k, d, k, d)
    return np.einsum("abcb->ac", view)


def distance(a, b):
    return float(np.abs(np.linalg.eigvalsh(a - b)).sum() / 2)


def entropy(rho):
    values = np.linalg.eigvalsh(rho)
    values = values[values > 1e-12]
    return float(-np.sum(values * np.log2(values)))


def pauli(n, assignments):
    return kron_all([assignments.get(j, I) for j in range(n)])


def cnot(n, control, target):
    u = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for j in range(2 ** n):
        out = bits(j, n)
        out[target] ^= out[control]
        u[integer(out), j] = 1
    return u


def incidence(n, edges):
    b = np.zeros((len(edges), n), dtype=np.int64)
    for e, (v, w) in enumerate(edges):
        if v == w or not (0 <= v < n and 0 <= w < n):
            raise ValueError("simple edges with distinct valid endpoints required")
        b[e, v] = b[e, w] = 1
    return b


def row_reduce(b):
    """GF(2) elimination with an invertible row-operation matrix."""
    reduced = b.copy()
    transform = np.eye(len(b), dtype=np.int64)
    rank = 0
    for col in range(b.shape[1]):
        pivots = np.flatnonzero(reduced[rank:, col])
        if not len(pivots):
            continue
        pivot = rank + pivots[0]
        reduced[[rank, pivot]] = reduced[[pivot, rank]]
        transform[[rank, pivot]] = transform[[pivot, rank]]
        for row in range(len(b)):
            if row != rank and reduced[row, col]:
                reduced[row] ^= reduced[rank]
                transform[row] ^= transform[rank]
        rank += 1
        if rank == len(b):
            break
    return reduced, transform, rank


def graph_encoding(n, edges):
    b = incidence(n, edges)
    e = len(edges)
    v = np.zeros((2 ** (n + e), 2 ** n), dtype=complex)
    for a in range(2 ** n):
        bit_string = bits(a, n)
        out = np.concatenate([bit_string, (b @ bit_string) % 2])
        v[integer(out), a] = 1
    gauss = [pauli(n + e, {x: Z, y: Z, n + j: Z})
             for j, (x, y) in enumerate(edges)]
    dressed_x = [pauli(n + e, {j: X, **{
        n + k: X for k, edge in enumerate(edges) if j in edge}})
        for j in range(n)]
    bare_z = [pauli(n + e, {j: Z}) for j in range(n)]
    return {"n": n, "edges": edges, "B": b, "V": v, "G": gauss,
            "X": dressed_x, "Z": bare_z}


def graph_circuit(n, edges):
    u = np.eye(2 ** (n + len(edges)), dtype=complex)
    for e, (v, w) in enumerate(edges):
        u = cnot(n + len(edges), v, n + e) @ u
        u = cnot(n + len(edges), w, n + e) @ u
    return u


def syndrome_channel(rho, b, reference_dimension=1):
    n = b.shape[1]
    labels = [tuple((b @ bits(a, n)) % 2) for a in range(2 ** n)]
    mask = np.array([[x == y for y in labels] for x in labels], dtype=float)
    return rho * np.kron(mask, np.ones((reference_dimension, reference_dimension)))


def linear_permutation(transform):
    n = len(transform)
    u = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for j in range(2 ** n):
        u[integer((transform @ bits(j, n)) % 2), j] = 1
    return u


def matched_memory_encoding(n):
    """Keep the original bare matching constraint, but retain N-1 memory qubits."""
    v = np.zeros((2 ** (2 * n - 1), 2 ** n), dtype=complex)
    for a in range(2 ** n):
        q = bits(a, n)
        out = np.concatenate([np.repeat(q[0], n), q[1:] ^ q[0]])
        v[integer(out), a] = 1
    return v


def bare_matching_projector(n, edges):
    out = np.eye(2 ** n, dtype=complex)
    eye = out.copy()
    for a, b in edges:
        out = out @ (eye + pauli(n, {a: Z, b: Z})) / 2
    return out


def graph_summary(n, edges):
    model = graph_encoding(n, edges)
    b, v = model["B"], model["V"]
    _, _, rank = row_reduce(b)
    psi = v @ kron_all([PLUS] * n)
    charge = partial_trace(density(psi), [2 ** n, 2 ** len(edges)], [1])
    all_syndromes = {tuple((b @ bits(a, n)) % 2) for a in range(2 ** n)}
    return {"vertices": n, "edges": len(edges), "kinematic_qubits": n + len(edges),
            "physical_dimension": v.shape[1], "independent_constraints": len(edges),
            "incidence_rank_GF2": rank, "distinct_relative_parities": len(all_syndromes),
            "charge_entropy_uniform_input_bits": entropy(charge),
            "CNOT_count": 2 * len(edges),
            "dressed_X_support_sizes": [1 + sum(j in edge for edge in edges)
                                       for j in range(n)],
            "bare_matched_dimension": int(round(np.trace(
                bare_matching_projector(n, edges)).real)),
            "minimum_retained_memory_dimension_for_bare_matching": 2 ** (n - 1)}


def report():
    model = graph_encoding(2, [(0, 1)])
    v, g = model["V"], model["G"][0]
    u = graph_circuit(2, [(0, 1)])
    rho = random_density(12)
    vr = np.kron(v, np.eye(3))
    encoded = vr @ rho @ vr.conj().T
    decoded = vr.conj().T @ encoded @ vr
    bare = partial_trace(encoded, [4, 2, 3], [0, 2])
    incoming = kron_all([PLUS, np.array([1, 0])])
    psi = v @ incoming
    local_a = partial_trace(density(psi), [2, 2, 2], [0])
    return {
        "round": 362,
        "scope": {
            "claim": "Exact reversible endpoint encodings with explicit compensating charge "
                     "or retained memory; complete dressed logical algebras and resource bounds.",
            "quantifiers": "Every finite AB (or N-vertex) input state, including every finite "
                           "reference and all instruments represented on the code; graph fixed.",
            "extra_inputs": ["complex qubits", "specified connected finite graph",
                             "specified Z2 matching constraints", "blank charge/memory registers",
                             "access to incident charge qubits", "ideal CNOT and dressed controls"],
            "not_claimed": ["fixed-Gauss-preserving admission from every bare input",
                            "preservation of each bare local marginal or bare local access",
                            "edge count is a universal minimal memory cost",
                            "zero physical time/energy/communication cost",
                            "derivation of gauge laws, space, or Einstein dynamics"],
            "prior_overlap": "Pair parity isometry is round 250 up to Hadamard on its record; "
                             "new audit concerns constrained logical factors and retained resources.",
        },
        "pair": {
            "isometry_error": float(np.linalg.norm(v.conj().T @ v - np.eye(4))),
            "constraint_error": float(np.linalg.norm(g @ v - v)),
            "reference_decoding_trace_distance": distance(decoded, rho),
            "bare_channel_formula_error": float(np.linalg.norm(
                bare - syndrome_channel(rho, model["B"], 3))),
            "plus_zero_bare_A_distance_from_input": distance(local_a, density(PLUS)),
            "plus_zero_bare_X_A_expectation": float(np.vdot(
                psi, pauli(3, {0: X}) @ psi).real),
            "plus_zero_dressed_X_A_expectation": float(np.vdot(
                psi, model["X"][0] @ psi).real),
            "encoder_fixed_G_commutator_operator_norm": float(np.linalg.norm(u @ g - g @ u, 2)),
            "constraint_intertwiner_error": float(np.linalg.norm(
                u @ pauli(3, {2: Z}) @ u.conj().T - g)),
        },
        "graphs": [graph_summary(2, [(0, 1)]),
                   graph_summary(3, [(0, 1), (1, 2)]),
                   graph_summary(3, [(0, 1), (1, 2), (2, 0)]),
                   graph_summary(4, [(0, 1), (1, 2), (2, 3), (3, 0)])],
        "next_interface": "A physical admission process needs a specified incoming charge sector, "
                          "a fixed total constraint, local access and a schedule; endpoint "
                          "isometries alone do not supply those dynamical resources.",
    }


class Checks(unittest.TestCase):
    def test_01_pair_isometry_and_cnot_construction(self):
        v = graph_encoding(2, [(0, 1)])["V"]
        blank = np.kron(np.eye(4), np.array([[1], [0]]))
        np.testing.assert_allclose(graph_circuit(2, [(0, 1)]) @ blank, v)
        np.testing.assert_allclose(v.conj().T @ v, np.eye(4))

    def test_02_full_physical_space_exactly_the_code(self):
        for n, edges in [(2, [(0, 1)]), (3, [(0, 1), (1, 2), (2, 0)])]:
            model = graph_encoding(n, edges)
            v = model["V"]
            p = np.eye(len(v), dtype=complex)
            for g in model["G"]:
                p = p @ (np.eye(len(v)) + g) / 2
            np.testing.assert_allclose(p, v @ v.conj().T)
            self.assertEqual(round(np.trace(p).real), 2 ** n)

    def test_03_arbitrary_mixed_reference_recovery(self):
        v = graph_encoding(2, [(0, 1)])["V"]
        for r in (1, 2, 3):
            rho = random_density(4 * r, 362 + r)
            vr = np.kron(v, np.eye(r))
            encoded = vr @ rho @ vr.conj().T
            np.testing.assert_allclose(vr.conj().T @ encoded @ vr, rho, atol=1e-14)
            np.testing.assert_allclose(partial_trace(encoded, [8, r], [1]),
                                       partial_trace(rho, [4, r], [1]), atol=1e-14)

    def test_04_maximal_reference_entanglement_preserved(self):
        v = graph_encoding(2, [(0, 1)])["V"]
        bell = np.eye(4).reshape(-1) / 2
        out = np.kron(v, np.eye(4)) @ bell
        self.assertAlmostEqual(entropy(partial_trace(density(out), [8, 4], [1])), 2)
        self.assertAlmostEqual(np.linalg.norm(out), 1)

    def test_05_complete_commuting_logical_factors(self):
        model = graph_encoding(2, [(0, 1)])
        v, g = model["V"], model["G"][0]
        local = []
        for j in (0, 1):
            logical = [np.eye(8), model["X"][j],
                       1j * model["X"][j] @ model["Z"][j], model["Z"][j]]
            for op, original in zip(logical, PAULI):
                np.testing.assert_allclose(op @ g, g @ op)
                np.testing.assert_allclose(op @ v, v @ pauli(2, {j: original}))
            local.append(logical)
        compressed = []
        for a, b in itertools.product(local[0], local[1]):
            np.testing.assert_allclose(a @ b, b @ a)
            compressed.append((v.conj().T @ a @ b @ v).ravel())
        self.assertEqual(np.linalg.matrix_rank(np.stack(compressed)), 16)

    def test_06_instrument_branches_and_reference(self):
        model = graph_encoding(2, [(0, 1)])
        v, u, g = model["V"], graph_circuit(2, [(0, 1)]), model["G"][0]
        eta = .37
        kraus = [np.diag([1, np.sqrt(1 - eta)]),
                 np.array([[0, np.sqrt(eta)], [0, 0]])]
        rho = random_density(12, 368)
        vr = np.kron(v, np.eye(3))
        total = np.zeros((8, 8), dtype=complex)
        for k in kraus:
            k_ab = np.kron(k, I)
            dressed = u @ np.kron(k_ab, I) @ u.conj().T
            total += dressed.conj().T @ dressed
            np.testing.assert_allclose(dressed @ g, g @ dressed)
            dr = np.kron(dressed, np.eye(3))
            kr = np.kron(k_ab, np.eye(3))
            np.testing.assert_allclose(dr @ vr @ rho @ vr.conj().T @ dr.conj().T,
                                       vr @ kr @ rho @ kr.conj().T @ vr.conj().T, atol=1e-14)
        np.testing.assert_allclose(total, np.eye(8), atol=1e-14)

    def test_07_bare_partial_trace_is_syndrome_pinching(self):
        for n, edges in [(2, [(0, 1)]), (3, [(0, 1), (1, 2), (2, 0)])]:
            model = graph_encoding(n, edges)
            r = 2
            rho = random_density((2 ** n) * r, 369 + n)
            vr = np.kron(model["V"], np.eye(r))
            bare = partial_trace(vr @ rho @ vr.conj().T,
                                 [2 ** n, 2 ** len(edges), r], [0, 2])
            np.testing.assert_allclose(bare, syndrome_channel(rho, model["B"], r),
                                       atol=1e-14)

    def test_08_bare_ability_lost_dressed_ability_kept(self):
        model = graph_encoding(2, [(0, 1)])
        v = model["V"]
        psi = v @ kron_all([PLUS, np.array([1, 0])])
        a = partial_trace(density(psi), [2, 2, 2], [0])
        np.testing.assert_allclose(a, I / 2, atol=1e-14)
        self.assertAlmostEqual(distance(a, density(PLUS)), .5)
        np.testing.assert_allclose(v.conj().T @ pauli(3, {0: X}) @ v, np.zeros((4, 4)))
        self.assertAlmostEqual(np.vdot(psi, model["X"][0] @ psi).real, 1)

    def test_09_encoder_intertwines_but_does_not_preserve_fixed_constraint(self):
        u = graph_circuit(2, [(0, 1)])
        g = pauli(3, {0: Z, 1: Z, 2: Z})
        np.testing.assert_allclose(u @ pauli(3, {2: Z}) @ u.conj().T, g)
        self.assertAlmostEqual(np.linalg.norm(u @ g - g @ u, 2), 2)
        incoming = np.eye(8)[:, 2]  # |0,1,0>, G=-1
        self.assertAlmostEqual(np.vdot(incoming, g @ incoming).real, -1)
        outgoing = u @ incoming
        self.assertAlmostEqual(np.vdot(outgoing, g @ outgoing).real, 1)

    def test_10_constraint_preserving_unitary_cannot_move_sector(self):
        rng = np.random.default_rng(372)
        raw = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
        h = raw + raw.conj().T
        g = pauli(3, {0: Z, 1: Z, 2: Z})
        h = (h + g @ h @ g) / 2
        vals, vecs = np.linalg.eigh(h)
        w = (vecs * np.exp(-.23j * vals)) @ vecs.conj().T
        psi = w @ np.eye(8)[:, 2]
        self.assertAlmostEqual(np.vdot(psi, (np.eye(8) + g) @ psi).real / 2, 0)
        np.testing.assert_allclose(w @ g, g @ w, atol=1e-14)

    def test_11_graph_dressed_intertwiners_and_shared_support(self):
        n, edges = 3, [(0, 1), (1, 2), (2, 0)]
        model = graph_encoding(n, edges)
        v = model["V"]
        for j in range(n):
            for op, original in [(model["X"][j], X), (model["Z"][j], Z)]:
                np.testing.assert_allclose(op @ v, v @ pauli(n, {j: original}))
                for g in model["G"]:
                    np.testing.assert_allclose(g @ op, op @ g)
        for j in range(n):
            for k in range(j):
                np.testing.assert_allclose(model["X"][j] @ model["X"][k],
                                           model["X"][k] @ model["X"][j])
        self.assertEqual([1 + sum(j in edge for edge in edges) for j in range(n)], [3]*3)

    def test_12_graph_constraint_rank_and_incidence_rank_differ(self):
        for n, edges in [(3, [(0, 1), (1, 2)]),
                         (3, [(0, 1), (1, 2), (2, 0)]),
                         (4, [(0, 1), (1, 2), (2, 3), (3, 0)])]:
            b = incidence(n, edges)
            self.assertEqual(row_reduce(b)[2], n - 1)
            stabilizers = np.concatenate([b, np.eye(len(edges), dtype=np.int64)], axis=1)
            self.assertEqual(row_reduce(stabilizers)[2], len(edges))
            kernel = [a for a in range(2 ** n) if not np.any((b @ bits(a, n)) % 2)]
            self.assertEqual(kernel, [0, 2 ** n - 1])

    def test_13_charge_entropy_has_only_n_minus_one_bits(self):
        for n, edges in [(2, [(0, 1)]), (3, [(0, 1), (1, 2), (2, 0)]),
                         (4, [(0, 1), (1, 2), (2, 3), (3, 0)])]:
            result = graph_summary(n, edges)
            self.assertAlmostEqual(result["charge_entropy_uniform_input_bits"], n - 1)
            self.assertEqual(result["distinct_relative_parities"], 2 ** (n - 1))

    def test_14_linear_charge_compression_is_reversible(self):
        model = graph_encoding(3, [(0, 1), (1, 2), (2, 0)])
        reduced, t, rank = row_reduce(model["B"])
        np.testing.assert_array_equal((t @ model["B"]) % 2, reduced)
        self.assertFalse(np.any(reduced[rank:]))
        uc = linear_permutation(t)
        np.testing.assert_allclose(uc.conj().T @ uc, np.eye(8))
        v_compressed = np.kron(np.eye(8), uc) @ model["V"]
        for col in range(8):
            occupied = np.flatnonzero(np.abs(v_compressed[:, col]) > .5)
            self.assertEqual(len(occupied), 1)
            self.assertEqual(occupied[0] & 1, 0)

    def test_15_old_bare_constraint_has_insufficient_dimension(self):
        for n, edges in [(2, [(0, 1)]), (3, [(0, 1), (1, 2)]),
                         (4, [(0, 1), (1, 2), (2, 3), (3, 0)])]:
            p = bare_matching_projector(n, edges)
            self.assertEqual(round(np.trace(p).real), 2)
            memory = 2 ** (n - 1)
            self.assertEqual(2 * memory, 2 ** n)
            self.assertLess(2 * (memory - 1), 2 ** n)

    def test_16_memory_lower_bound_is_achievable_with_reference(self):
        n = 3
        v = matched_memory_encoding(n)
        np.testing.assert_allclose(v.conj().T @ v, np.eye(2 ** n))
        p = np.kron(bare_matching_projector(n, [(0, 1), (1, 2)]), np.eye(4))
        np.testing.assert_allclose(p @ v, v)
        rho = random_density(16, 378)
        vr = np.kron(v, I)
        np.testing.assert_allclose(vr.conj().T @ (vr @ rho @ vr.conj().T) @ vr, rho)

    def test_17_memory_carries_relocated_original_abilities(self):
        n = 3
        v = matched_memory_encoding(n)
        total = 2 * n - 1
        root_x = pauli(total, {j: X for j in range(total)})
        np.testing.assert_allclose(root_x @ v, v @ pauli(n, {0: X}))
        for j in range(1, n):
            mem = n + j - 1
            np.testing.assert_allclose(pauli(total, {mem: X}) @ v,
                                       v @ pauli(n, {j: X}))
            np.testing.assert_allclose(pauli(total, {0: Z, mem: Z}) @ v,
                                       v @ pauli(n, {j: Z}))

    def test_18_reset_discard_loses_unknown_state_information(self):
        v = matched_memory_encoding(2)
        states = [density(v[:, a]) for a in (0, 1)]
        self.assertAlmostEqual(distance(states[0], states[1]), 1)
        bare = [partial_trace(rho, [4, 2], [0]) for rho in states]
        self.assertAlmostEqual(distance(bare[0], bare[1]), 0)
        original_v = graph_encoding(2, [(0, 1)])["V"]
        rho = density(original_v @ kron_all([PLUS, np.array([1, 0])]))
        charge = partial_trace(rho, [4, 2], [1])
        self.assertAlmostEqual(entropy(charge), 1)


if __name__ == "__main__":
    main(__name__, "boundary_charge_encoding_audit", report)
