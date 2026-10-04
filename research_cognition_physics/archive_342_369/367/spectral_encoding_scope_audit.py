"""Round 367: target-dependent spectral encoding versus physical emergence."""
from itertools import combinations, islice
from math import ceil, comb, log2
import unittest
import numpy as np
from growing_stream_audit import main


def unitary(h, time):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*time*values)) @ vectors.conjugate().T


def compile_spectrum(h, bins):
    if bins < 1 or int(bins) != bins:
        raise ValueError("bins must be a positive integer")
    h = np.asarray(h, complex)
    if h.ndim != 2 or h.shape[0] != h.shape[1] or not np.allclose(h, h.conjugate().T):
        raise ValueError("finite Hermitian target required")
    dimension = len(h)
    values, vectors = np.linalg.eigh(h)
    r = max(1, ceil(log2(dimension)))
    n = bins+2*r
    width = float(values[-1]-values[0])
    delta = width/bins if width > 0 else 1.
    weights = r+np.rint((values-values[0])/delta).astype(int)
    used, indices = {}, []
    for weight in weights:
        weight = int(weight)
        slot = used.get(weight, 0)
        positions = next(islice(combinations(range(n), weight), slot, slot+1))
        indices.append(sum(1 << bit for bit in positions))
        used[weight] = slot+1
    energies = values[0]+delta*(weights-r)
    effective = (vectors*energies) @ vectors.conjugate().T
    return {"D": dimension, "N": n, "r": r, "M": bins, "delta": delta,
            "minimum": float(values[0]), "target_values": values, "vectors": vectors,
            "weights": weights, "indices": indices, "compiled_values": energies,
            "effective": effective, "error": float(np.linalg.norm(effective-h, 2))}


def explicit_isometry(plan):
    if plan["N"] > 12:
        raise ValueError("Only small full-Hilbert-space cross-checks are allocated")
    v = np.zeros((2**plan["N"], plan["D"]), complex)
    v[plan["indices"], :] = plan["vectors"].conjugate().T
    return v


def physical_energies(plan):
    return plan["minimum"]+plan["delta"]*np.array(
        [index.bit_count()-plan["r"] for index in range(2**plan["N"])])


def pure_distance(a, b):
    return float(np.sqrt(max(0., 1-abs(np.vdot(a, b))**2)))


def trace_distance(a, b):
    delta = a-b
    return float(np.sum(np.abs(np.linalg.eigvalsh((delta+delta.conjugate().T)/2)))/2)


def target_matrix():
    rng = np.random.default_rng(367)
    raw = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    q, _ = np.linalg.qr(raw)
    return (q*np.array([-.25, .12, .68, 1.32])) @ q.conjugate().T


def one_qubit_purity(state, qubits, position):
    matrix = np.moveaxis(state.reshape((2,)*qubits), position, 0).reshape(2, -1)
    rho = matrix @ matrix.conjugate().T
    return float(np.trace(rho @ rho).real)


def encoded_entangling_example():
    z = np.diag([1., -1.])
    h = np.kron(z, z)
    plan = compile_spectrum(h, 4)
    v = explicit_isometry(plan)
    initial = np.full(4, .5, complex)
    encoded = v @ initial
    time = np.pi/4
    final = np.exp(-1j*time*physical_energies(plan))*encoded
    decoded = v.conjugate().T @ final
    return {"N": plan["N"], "weights": plan["weights"].tolist(),
            "physical_time": plan["delta"]*time,
            "logical_purity_initial": one_qubit_purity(initial, 2, 0),
            "logical_purity_final": one_qubit_purity(decoded, 2, 0),
            "physical_one_qubit_purity_max_change": max(abs(
                one_qubit_purity(encoded, plan["N"], i)-
                one_qubit_purity(final, plan["N"], i)) for i in range(plan["N"])),
            "target_vector_error": float(np.linalg.norm(decoded-unitary(h, time) @ initial)),
            "minimum_excitation_gap_for_logical_A_flip": int(plan["M"])}


def periodic_channel_witness():
    period = 2*np.pi
    h = np.diag([0., np.sqrt(2)])
    plus = np.ones(2)/np.sqrt(2)
    return {"fixed_source_period": period, "target_gap": float(np.sqrt(2)),
            "target_plus_distance_at_period": pure_distance(plus, unitary(h, period) @ plus),
            "epsilon_0_plus_epsilon_period_lower": float(abs(np.sin(np.pi*np.sqrt(2))))}


def report():
    h, time = target_matrix(), 1.3
    rows = []
    for bins in (8, 16, 32, 64, 128):
        plan = compile_spectrum(h, bins)
        desired, actual = unitary(h, time), unitary(plan["effective"], time)
        reference_distance = pure_distance(desired.reshape(-1)/2, actual.reshape(-1)/2)
        rows.append({"M_bins": bins, "physical_qubits_N": plan["N"],
                     "spectral_spacing_delta": plan["delta"],
                     "excitation_weights": plan["weights"].tolist(),
                     "distinct_codewords": len(set(plan["indices"])),
                     "generator_error": plan["error"],
                     "generator_error_bound": plan["delta"]/2,
                     "reference_witness_distance_at_t": reference_distance,
                     "half_diamond_upper_at_t": min(1., time*plan["delta"]/2),
                     "physical_elapsed_time": time*plan["delta"],
                     "minimum_bucket_capacity": str(comb(plan["N"], plan["r"]))})
    return {"round": 367,
            "scope": "Invariant-code dynamics only, with a target-dependent global encoding and scale-dependent time identification. Not a local low-energy Hamiltonian simulation, thermodynamic equivalence, autonomous controller, or generation of GR.",
            "hypothesis_tested": "Can unrestricted static encodings and freely calibrated clocks make a noninteracting source reproduce every prescribed finite Hermitian dynamics?",
            "compiled_target_dimension": 4, "comparison_time": time,
            "spectral_compilation": rows,
            "encoded_entangling_gate": encoded_entangling_example(),
            "fixed_clock_periodic_obstruction": periodic_channel_witness(),
            "proved": ["Every finite Hermitian target admits an invariant spectral code with reference-uniform finite-time error.",
                       "The code dimension is supported by bulk binomial degeneracy, including rounded eigenvalue collisions.",
                       "For an exact two-level target, nonzero logical transitions across excitation gap M require support at least M.",
                       "Any fixed CPTP encoding and decoding preserves the source's common physical-time period."],
            "resource_account": {"source": "H0_N=sum_i |1><1|_i; no physical intersite interactions",
                                 "qubits": "N=M+2*max(1,ceil(log2 D))",
                                 "source_norm": "N", "logical_time_identification": "tau=(W/M)*t for W>0",
                                 "unpaid_physical_interfaces": ["target eigensystem computation",
                                     "global state encoding and observable decoding",
                                     "selection and protection of the middle-spectrum code",
                                     "clock calibration and composition across targets"]},
            "sources": ["https://arxiv.org/html/1701.05182",
                        "https://link.aps.org/accepted/10.1103/PhysRevLett.114.031104"]}


class Checks(unittest.TestCase):
    def test_01_bulk_binomial_capacity_includes_all_degeneracies(self):
        for dimension in range(1, 17):
            r = max(1, ceil(log2(dimension)))
            for bins in (1, 2, 7):
                n = bins+2*r
                self.assertGreaterEqual(comb(n, r), 2**r)
                for k in range(r, n-r+1):
                    self.assertGreaterEqual(comb(n, k), dimension)

    def test_02_collisions_and_constant_target_keep_distinct_codewords(self):
        for h in (np.diag([0., 0., .01, 1.]), 2*np.eye(5)):
            plan = compile_spectrum(h, 2)
            self.assertEqual(len(set(plan["indices"])), len(h))
            self.assertEqual([x.bit_count() for x in plan["indices"]], plan["weights"].tolist())
            self.assertLessEqual(plan["error"], plan["delta"]/2+1e-14)
        self.assertLess(compile_spectrum(2*np.eye(5), 2)["error"], 1e-14)

    def test_03_full_space_isometry_and_invariant_generator(self):
        h = np.array([[.4, .2j], [-.2j, 1.1]])
        plan = compile_spectrum(h, 4); v = explicit_isometry(plan)
        energies = physical_energies(plan)
        np.testing.assert_allclose(v.conjugate().T @ v, np.eye(2), atol=1e-14)
        np.testing.assert_allclose(energies[:, None]*v, v @ plan["effective"], atol=1e-14)
        np.testing.assert_allclose(v.conjugate().T @ (energies[:, None]*v),
                                   plan["effective"], atol=1e-14)

    def test_04_general_hermitian_spectral_error_bound(self):
        for bins in (3, 8, 32, 128):
            plan = compile_spectrum(target_matrix(), bins)
            self.assertLessEqual(plan["error"], plan["delta"]/2+2e-14)
            self.assertEqual(len(set(plan["indices"])), 4)

    def test_05_full_space_unknown_reference_evolution(self):
        h = np.array([[.4, .2j], [-.2j, 1.1]])
        plan = compile_spectrum(h, 4); v = explicit_isometry(plan)
        rng = np.random.default_rng(1367)
        psi = rng.normal(size=(2, 3))+1j*rng.normal(size=(2, 3))
        psi /= np.linalg.norm(psi)
        for time in (-1.2, .5, 2.):
            actual = np.exp(-1j*time*physical_energies(plan))[:, None]*(v @ psi)
            ideal = v @ (unitary(h, time) @ psi)
            self.assertLessEqual(np.linalg.norm(actual-ideal), abs(time)*plan["delta"]/2+1e-14)
            self.assertLessEqual(pure_distance(actual.ravel(), ideal.ravel()),
                                 abs(time)*plan["delta"]/2+1e-12)

    def test_06_choi_witness_and_uniform_time_bound(self):
        h = target_matrix()
        for bins in (8, 16, 32):
            plan = compile_spectrum(h, bins)
            for time in (.2, .7, 1.3):
                u = unitary(h, time); w = unitary(plan["effective"], time)
                self.assertLessEqual(np.linalg.norm(u-w, 2), abs(time)*plan["delta"]/2+1e-14)
                self.assertLessEqual(pure_distance(u.ravel()/2, w.ravel()/2),
                                     abs(time)*plan["delta"]/2+1e-12)

    def test_07_encoded_entanglement_without_physical_cut_entanglement_growth(self):
        row = encoded_entangling_example()
        self.assertAlmostEqual(row["logical_purity_initial"], 1., places=14)
        self.assertAlmostEqual(row["logical_purity_final"], .5, places=14)
        self.assertLess(row["physical_one_qubit_purity_max_change"], 2e-15)
        self.assertLess(row["target_vector_error"], 2e-15)

    def test_08_logical_flip_needs_extensive_excitation_support(self):
        for bins in (2, 4, 8, 16):
            plan = compile_spectrum(np.diag([-1., 1.]), bins)
            a, b = plan["indices"]
            self.assertEqual(abs(a.bit_count()-b.bit_count()), bins)
            self.assertGreaterEqual((a ^ b).bit_count(), bins)
        # Every one/two-site Pauli flip has zero matrix element for M=4.
        plan = compile_spectrum(np.diag([-1., 1.]), 4)
        a, b = plan["indices"]
        for size in (0, 1, 2, 3):
            for positions in combinations(range(plan["N"]), size):
                flip = sum(1 << p for p in positions)
                self.assertNotEqual(a ^ flip, b)

    def test_09_nonlocal_CPTP_maps_cannot_remove_fixed_period(self):
        rng = np.random.default_rng(2367)
        raw = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
        u, _ = np.linalg.qr(raw); encoder = u[:, :2]
        rho = np.array([[.6, .1j], [-.1j, .4]])
        encoded = encoder @ rho @ encoder.conjugate().T
        phase = np.exp(-2j*np.pi*np.array([0, 1, 1, 2]))
        evolved = phase[:, None]*encoded*phase.conjugate()[None, :]
        # Decode by tracing out the second physical qubit: generally not inverse.
        first = np.trace(encoded.reshape(2, 2, 2, 2), axis1=1, axis2=3)
        second = np.trace(evolved.reshape(2, 2, 2, 2), axis1=1, axis2=3)
        np.testing.assert_allclose(first, second, atol=1e-14)

    def test_10_fixed_clock_target_obstruction_is_finite(self):
        row = periodic_channel_witness()
        self.assertAlmostEqual(row["target_plus_distance_at_period"],
                               row["epsilon_0_plus_epsilon_period_lower"], places=14)
        self.assertGreater(row["target_plus_distance_at_period"], .96)

    def test_11_changing_time_scale_changes_the_logical_period(self):
        for bins in (8, 16, 32):
            plan = compile_spectrum(target_matrix(), bins)
            period = 2*np.pi/plan["delta"]
            relative = plan["compiled_values"]-plan["minimum"]
            np.testing.assert_allclose(np.exp(-1j*period*relative), 1., atol=2e-13)
        self.assertGreater(2*np.pi/compile_spectrum(target_matrix(), 32)["delta"], 100.)

    def test_12_small_full_tensor_matches_bare_product_evolution(self):
        plan = compile_spectrum(np.diag([-.2, .7]), 4)
        n, tau = plan["N"], .43
        local = np.diag([1., np.exp(-1j*tau)])
        full = np.array([[1.+0j]])
        for _ in range(n):
            full = np.kron(full, local)
        diagonal = np.exp(-1j*tau*np.array([k.bit_count() for k in range(2**n)]))
        np.testing.assert_allclose(full, np.diag(diagonal), atol=2e-15)


if __name__ == "__main__":
    main(__name__, "spectral_encoding_scope_audit", report)
