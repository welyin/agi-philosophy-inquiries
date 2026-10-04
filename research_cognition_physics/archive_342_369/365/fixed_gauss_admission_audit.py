"""Round 365: a genuine fixed-Gauss admission/interaction path and its clock.

The graph, matter/link encoding, neutral bridge, couplings and initial clock
are inputs. Autonomous execution reuses the frozen round-346 construction;
it neither selects a natural graph nor closes preparation/readout costs.
"""
from functools import lru_cache
import itertools
import unittest
import numpy as np
from growing_stream_audit import main
from boundary_charge_encoding_audit import (
    I, X, Y, Z, PLUS, bits, integer, kron_all, density, random_density,
    partial_trace, distance, entropy, pauli, incidence)
from autonomous_history_clock_audit import (
    exponential, hop_weights, bare_clock, capped_omega)


ZERO = np.array([1., 0.], complex)
ONE = np.array([0., 1.], complex)
EDGES = ((0, 1), (2, 3), (1, 2))  # logical A, B, bridge C


@lru_cache(maxsize=None)
def gauge_model(vertices=4, edges=EDGES):
    if len(set(tuple(sorted(edge)) for edge in edges)) != len(edges):
        raise ValueError("The two-body obstruction assumes no parallel edges.")
    b = incidence(vertices, edges).T
    nlinks, total = len(edges), vertices + len(edges)
    embedding = np.zeros((2 ** total, 2 ** nlinks), complex)
    for a in range(2 ** nlinks):
        link = bits(a, nlinks)
        embedding[integer(np.r_[(b @ link) % 2, link]), a] = 1.
    gauss = [pauli(total, {v: Z, **{
        vertices + e: Z for e, edge in enumerate(edges) if v in edge}})
        for v in range(vertices)]
    flip = [pauli(total, {v: X, w: X, vertices + e: X})
            for e, (v, w) in enumerate(edges)]
    phase = [pauli(total, {vertices + e: Z}) for e in range(nlinks)]
    projector = np.eye(2 ** total, dtype=complex)
    for g in gauss:
        projector = projector @ (np.eye(2 ** total) + g) / 2
    return {"vertices": vertices, "edges": edges, "incidence": b,
            "V": embedding, "G": gauss, "X": flip, "Z": phase,
            "P": projector, "qubits": total}


def blank_embedding(model=None, bit=0):
    model = gauge_model() if model is None else model
    return model["V"] @ np.kron(np.eye(4), np.array([[1], [0]]) if bit == 0
                                else np.array([[0], [1]]))


@lru_cache(maxsize=None)
def generators(logical=False):
    if logical:
        flip = [pauli(3, {e: X}) for e in range(3)]
        phase = [pauli(3, {e: Z}) for e in range(3)]
        eye = np.eye(8)
    else:
        model = gauge_model()
        flip, phase = model["X"], model["Z"]
        eye = np.eye(len(model["V"]))
    h = (flip[2] + phase[2]) / np.sqrt(2)
    d = (eye - phase[0]) @ (eye - phase[2]) / 4
    q = phase[2] @ phase[1]
    return {"eye": eye, "H": h, "D": d, "Q": q,
            "hamiltonians": (eye - h, d, q)}


def pulse(kind, angle, logical=False, spectral=False):
    gen = generators(logical)
    if spectral:
        return exponential(gen["hamiltonians"][kind], angle)
    eye = gen["eye"]
    if kind == 0:
        return np.exp(-1j * angle) * (
            np.cos(angle) * eye + 1j * np.sin(angle) * gen["H"])
    if kind == 1:
        return eye + np.expm1(-1j * angle) * gen["D"]
    return np.cos(angle) * eye - 1j * np.sin(angle) * gen["Q"]


def schedule(theta):
    return ((0, np.pi / 2), (1, np.pi), (0, np.pi / 2), (2, theta),
            (0, np.pi / 2), (1, np.pi), (0, np.pi / 2))


def gates(theta, logical=False):
    return [pulse(kind, angle, logical) for kind, angle in schedule(theta)]


def cumulative(gate_list):
    values = [np.eye(len(gate_list[0]), dtype=complex)]
    for gate in gate_list:
        values.append(gate @ values[-1])
    return values


def old_gate(theta):
    return np.diag(np.exp(-1j * theta * np.diag(np.kron(Z, Z))))


def stage_metrics(theta=np.pi / 4):
    model = gauge_model()
    initial = blank_embedding() @ kron_all([PLUS, PLUS])
    result = []
    for stage, u in enumerate(cumulative(gates(theta))):
        physical = u @ initial
        logical = model["V"].conj().T @ physical
        link = partial_trace(density(physical), [2] * 7, [6])
        logical_link = partial_trace(density(logical), [4, 2], [1])
        result.append({
            "stage": stage,
            "gauss_leakage_probability": float(np.linalg.norm(
                (np.eye(128) - model["P"]) @ physical) ** 2),
            "bare_bridge_entropy_bits": entropy(link),
            "logical_bridge_entropy_bits": entropy(logical_link),
            "bridge_Z": float(np.vdot(physical, model["Z"][2] @ physical).real),
            "matter_Z": [float(np.vdot(physical,
                pauli(7, {v: Z}) @ physical).real) for v in range(4)]})
    return result


def two_body_allowed(vertices=4, edges=EDGES):
    total = vertices + len(edges)
    b = incidence(vertices, edges).T
    accepted, examined = [], 0
    for weight in range(3):
        for positions in itertools.combinations(range(total), weight):
            for values in itertools.product((1, 2, 3), repeat=weight):
                pattern = np.zeros(total, dtype=int)
                pattern[list(positions)] = values
                flips = np.isin(pattern, (1, 2)).astype(int)
                syndrome = (flips[:vertices] + b @ flips[vertices:]) % 2
                examined += 1
                if not np.any(syndrome):
                    accepted.append(pattern.tolist())
    return examined, accepted


@lru_cache(maxsize=None)
def history(theta=.37, omega=.5, logical=True):
    word = gates(theta, logical)
    size, length = len(word[0]), len(word)
    weights = hop_weights(length, omega)
    h = np.zeros(((length + 1) * size,) * 2, complex)
    path = cumulative(word)
    for j, gate in enumerate(word):
        lower = slice(j * size, (j + 1) * size)
        upper = slice((j + 1) * size, (j + 2) * size)
        h[upper, lower] = weights[j] * gate
        h[lower, upper] = weights[j] * gate.conj().T
    return {"H": h, "weights": weights, "word": word,
            "cumulative": path, "omega": omega, "length": length,
            "data_dimension": size}


def clock_injection(position, data_dimension=8):
    clock = np.eye(8)[:, [position]]
    return np.kron(clock, np.eye(data_dimension))


def history_report(theta=.37, cap=1.):
    omega = capped_omega(7, cap)
    model = history(theta, omega)
    total_time = np.pi / omega
    initial = clock_injection(0)
    target = clock_injection(7) @ model["cumulative"][-1]
    actual = exponential(model["H"], total_time) @ initial
    times = (0., total_time / 2, total_time, 2 * total_time)
    psi0 = initial @ kron_all([PLUS, PLUS, ZERO])
    snapshots = []
    for time in times:
        state = exponential(model["H"], time) @ psi0
        hstate = model["H"] @ state
        snapshots.append({
            "time": float(time),
            "clock_endpoint_probability": float(
                np.linalg.norm(state.reshape(8, 8)[7]) ** 2),
            "total_energy_mean": float(np.vdot(state, hstate).real),
            "total_energy_second_moment": float(np.vdot(hstate, hstate).real)})
    return {"length": 7, "clock_levels": 8, "fixed_word_no_program_register": True,
            "full_clock_data_dimension": 8 * 128,
            "invariant_sector_dimension_used_for_exponential": 8 * 8,
            "max_clock_hop": float(max(model["weights"])),
            "omega": float(omega), "arrival_time": float(total_time),
            "Hamiltonian_operator_norm": float(omega * 7 / 2),
            "initial_energy_variance": float(omega ** 2 * 7 / 4),
            "physical_sector_full_input_operator_transfer_error": float(np.linalg.norm(
                actual - (-1j) ** 7 * target)),
            "return_operator_error": float(np.linalg.norm(
                exponential(model["H"], 2 * total_time) + np.eye(64))),
            "snapshots": snapshots,
            "timing_offset_0_1_endpoint_probability": float(
                np.cos(omega * .1 / 2) ** 14)}


def report():
    model = gauge_model()
    theta = np.pi / 4
    u = cumulative(gates(theta))[-1]
    w = blank_embedding()
    expected = w @ old_gate(theta)
    output = old_gate(theta) @ kron_all([PLUS, PLUS])
    selected = [0., .19, np.pi / 8, np.pi / 4]
    reference = random_density(12, 365)
    actual_map = np.kron(u @ w, np.eye(3))
    target_map = np.kron(expected, np.eye(3))
    reference_error = distance(actual_map @ reference @ actual_map.conj().T,
                               target_map @ reference @ target_map.conj().T)
    examined, allowed = two_body_allowed()
    return {
        "round": 365,
        "scope": {
            "claim": "A prepared neutral bridge admits exact continuous interactions "
                     "between existing unknown logical objects without leaving one "
                     "fixed Z2 Gauss sector; an internal history clock can execute the word.",
            "quantifiers": "Any finite old logical density matrix and finite reference; "
                           "general finite simple input graph for the algebraic encoding, "
                           "with selected old edges adjacent to a supplied bridge.",
            "inputs": ["specified graph, matter/link qubits and Gauss generators",
                       "old inputs already encoded in their physical sectors",
                       "one initially blank neutral bridge link",
                       "three-body gauge-invariant controls and adjacent-link diagonal controls",
                       "pulse angles or a fixed hardwired seven-gate history processor",
                       "prepared eight-level clock and engineered hopping weights",
                       "arrival-time calibration and future readout/storage remain inputs"],
            "not_claimed": ["admission of arbitrary bare nonphysical inputs",
                            "the same constraint representation as round 362",
                            "unchanged old marginals after intended interaction",
                            "free disposal of an occupied bridge",
                            "a two-body spatially local autonomous implementation",
                            "autonomous preparation, permanent output or free processor reuse",
                            "selection of graph, Z2 gauge law, dimension or Einstein dynamics"],
            "prior_relation": "Addresses round-362 fixed-constraint gap; reuses known "
                              "gauge-matter dressing and round-346 clock compilation, "
                              "not a new claim to invent either construction."},
        "model": {"vertices": 4, "old_links": 2, "added_bridge_links": 1,
                  "kinematic_qubits": 7, "kinematic_dimension": 128,
                  "independent_Gauss_constraints": 4, "physical_dimension": 8,
                  "old_input_dimension": 4, "edges_in_A_B_C_order": EDGES},
        "interaction": {
            "theta": float(theta),
            "endpoint_isometry_error": float(np.linalg.norm(u @ w - expected)),
            "arbitrary_mixed_reference_trace_distance": reference_error,
            "old_pair_entanglement_entropy_bits": entropy(
                partial_trace(density(output), [2, 2], [0])),
            "old_A_Y_for_B_Z_plus_minus": [float(np.sin(2 * theta)),
                                           float(-np.sin(2 * theta))],
            "angle_checks": [{"theta": float(t),
                              "old_A_reduced_purity": float((1 + np.cos(2 * t) ** 2) / 2)}
                             for t in selected]},
        "stage_data_plus_plus": stage_metrics(),
        "scheduled_resources": {
            "primitive_pulses": 7,
            "Hamiltonian_norms_at_unit_strength_H_D_phase": [2., 1., 1.],
            "Hadamard_pulse_support": ["matter_1", "matter_2", "bridge_C"],
            "diagonal_pulse_support": [["old_A_link", "bridge_C"],
                                       ["bridge_C", "old_B_link"]],
            "largest_Pauli_term_weight": 3,
            "duration_at_unit_strength": float(4 * np.pi + abs(theta)),
            "integrated_operator_norm_at_unit_strength": float(6 * np.pi + abs(theta)),
            "no_time_or_energy_optimality_claim": True},
        "two_body_exact_constraint_obstruction": {
            "Pauli_patterns_of_weight_at_most_two_examined": examined,
            "commuting_patterns": len(allowed),
            "all_commuting_patterns_are_I_Z": all(
                not any(p in (1, 2) for p in pattern) for pattern in allowed),
            "scope": "Forbids off-diagonal logical flips with exact two-body generators; "
                     "does NOT forbid all interactions or prove the seven-pulse protocol "
                     "optimal. In this minimal graph Z_m1 Z_m2 acts as Z_A Z_B."},
        "autonomous_execution": history_report(),
        "sources": [
            "https://arxiv.org/html/2205.08541",
            "https://arxiv.org/html/1712.07395",
            "https://arxiv.org/html/quant-ph/0309131"],
        "numerical_scope": "Direct finite matrices, independent closed-form/spectral pulse "
                           "exponentials, physical-sector intertwining and a 64-dimensional "
                           "clock-sector exponential. Round-346 exact transfer identity "
                           "supplies all-input/all-reference quantifiers."}


class Checks(unittest.TestCase):
    def test_01_physical_sector_dimension_and_isometry(self):
        for vertices, edges in ((2, ((0, 1),)), (3, ((0, 1), (1, 2))),
                                (3, ((0, 1), (1, 2), (2, 0))), (4, EDGES)):
            model = gauge_model(vertices, edges)
            v = model["V"]
            np.testing.assert_allclose(v.conj().T @ v, np.eye(2 ** len(edges)))
            np.testing.assert_allclose(v @ v.conj().T, model["P"])
            self.assertEqual(round(np.trace(model["P"]).real), 2 ** len(edges))

    def test_02_full_logical_Pauli_ability_and_Gauss(self):
        model = gauge_model()
        for e in range(3):
            for physical, logical in ((model["X"][e], X), (model["Z"][e], Z),
                    (1j * model["X"][e] @ model["Z"][e], Y)):
                np.testing.assert_allclose(physical @ model["V"],
                    model["V"] @ pauli(3, {e: logical}))
                for g in model["G"]:
                    np.testing.assert_allclose(physical @ g, g @ physical)

    def test_03_initial_old_objects_and_neutral_bridge_are_already_valid(self):
        w = blank_embedding()
        for a, b in itertools.product((0, 1), repeat=2):
            independent = np.zeros(128, complex)
            independent[integer((a, a, b, b, a, b, 0))] = 1
            np.testing.assert_array_equal(w[:, 2 * a + b], independent)
        for g in gauge_model()["G"]:
            np.testing.assert_allclose(g @ w, w)

    def test_04_every_Hamiltonian_commutes_with_fixed_constraints(self):
        for h in generators()["hamiltonians"]:
            np.testing.assert_allclose(h, h.conj().T)
            for g in gauge_model()["G"]:
                np.testing.assert_allclose(h @ g, g @ h)

    def test_05_partial_pulse_closed_formula_vs_direct_exponential(self):
        for kind, angle in itertools.product(range(3), (-.37, .23, 1.1)):
            np.testing.assert_allclose(pulse(kind, angle),
                pulse(kind, angle, spectral=True), atol=3e-14)
            u = pulse(kind, angle)
            np.testing.assert_allclose(u.conj().T @ u, np.eye(128), atol=3e-14)

    def test_06_seven_gate_identity_for_whole_data_space(self):
        triple = pauli(7, {4: Z, 5: Z, 6: Z})
        for theta in (-.7, 0., .23, np.pi / 4, np.pi / 2):
            actual = cumulative(gates(theta))[-1]
            target = np.cos(theta) * np.eye(128) - 1j * np.sin(theta) * triple
            np.testing.assert_allclose(actual, target, atol=3e-14)

    def test_07_every_partial_time_preserves_whole_input_sector(self):
        cumulative_u = np.eye(128, dtype=complex)
        model = gauge_model()
        for kind, angle in schedule(.61):
            for fraction in (0., .13, .47, .88, 1.):
                u = pulse(kind, fraction * angle) @ cumulative_u
                np.testing.assert_allclose(model["P"] @ u @ model["V"],
                                           u @ model["V"], atol=3e-14)
            cumulative_u = pulse(kind, angle) @ cumulative_u

    def test_08_mixed_reference_endpoint_and_inverse(self):
        w = blank_embedding()
        u = cumulative(gates(.47))[-1]
        for reference_dim in (1, 2, 3):
            rho = random_density(4 * reference_dim, 36508 + reference_dim)
            actual = np.kron(u @ w, np.eye(reference_dim))
            target = np.kron(w @ old_gate(.47), np.eye(reference_dim))
            out = actual @ rho @ actual.conj().T
            np.testing.assert_allclose(out, target @ rho @ target.conj().T, atol=2e-14)
            decoder = np.kron(w.conj().T @ u.conj().T, np.eye(reference_dim))
            np.testing.assert_allclose(decoder @ out @ decoder.conj().T, rho,
                                       atol=2e-14)
            np.testing.assert_allclose(partial_trace(out, [128, reference_dim], [1]),
                partial_trace(rho, [4, reference_dim], [1]), atol=2e-14)

    def test_09_entanglement_involves_old_objects_not_just_bridge(self):
        w = blank_embedding()
        for theta in (.19, np.pi / 8, np.pi / 4):
            physical = cumulative(gates(theta))[-1] @ w @ kron_all([PLUS, PLUS])
            old = w.conj().T @ physical
            reduced = partial_trace(density(old), [2, 2], [0])
            self.assertAlmostEqual(float(np.trace(reduced @ reduced).real),
                                   (1 + np.cos(2 * theta) ** 2) / 2)
        self.assertAlmostEqual(entropy(reduced), 1.)
        bridge_only = pulse(0, .43) @ w @ kron_all([PLUS, PLUS])
        logical = gauge_model()["V"].conj().T @ bridge_only
        old_only = partial_trace(density(logical), [4, 2], [0])
        np.testing.assert_allclose(old_only, density(kron_all([PLUS, PLUS])), atol=1e-14)

    def test_10_old_information_causes_actual_other_old_response(self):
        model = gauge_model()
        y_a = 1j * model["X"][0] @ model["Z"][0]
        theta = .37
        u = cumulative(gates(theta))[-1]
        for b, sign in ((ZERO, 1), (ONE, -1)):
            out = u @ blank_embedding() @ kron_all([PLUS, b])
            self.assertAlmostEqual(np.vdot(out, y_a @ out).real,
                                   sign * np.sin(2 * theta))

    def test_11_occupied_bridge_must_not_be_discarded_or_reset(self):
        model = gauge_model()
        compute = cumulative(gates(.37))[3]
        state = compute @ blank_embedding() @ kron_all([PLUS, ZERO])
        logical = model["V"].conj().T @ state
        old = partial_trace(density(logical), [4, 2], [0])
        self.assertAlmostEqual(distance(old, density(kron_all([PLUS, ZERO]))), .5)
        bridge = partial_trace(density(state), [64, 2], [1])
        self.assertAlmostEqual(entropy(bridge), 1.)
        kept = partial_trace(density(state), [64, 2], [0])
        illegal_reset = np.kron(kept, density(ZERO))
        self.assertAlmostEqual(np.trace(model["P"] @ illegal_reset).real, .5)
        final = cumulative(gates(.37))[-1] @ blank_embedding()
        np.testing.assert_allclose(final, blank_embedding() @ old_gate(.37), atol=1e-14)

    def test_12_nonblank_logical_bridge_changes_old_operation(self):
        u = cumulative(gates(.37))[-1]
        w1 = blank_embedding(bit=1)
        np.testing.assert_allclose(u @ w1, w1 @ old_gate(-.37), atol=1e-14)
        self.assertGreater(np.linalg.norm(u @ w1 - w1 @ old_gate(.37)), .5)

    def test_13_postinteraction_instrument_and_reference_abilities(self):
        model = gauge_model()
        eta = .31
        k0 = np.diag([1., np.sqrt(1 - eta)])
        k1 = np.array([[0., np.sqrt(eta)], [0., 0.]])
        xa, za = model["X"][0], model["Z"][0]
        ya = 1j * xa @ za
        dressed = [((1 + np.sqrt(1 - eta)) * np.eye(128) +
                    (1 - np.sqrt(1 - eta)) * za) / 2,
                   np.sqrt(eta) * (xa + 1j * ya) / 2]
        w = blank_embedding()
        rho = random_density(8, 36513)
        source = np.kron(w @ old_gate(.37), I)
        for d, k in zip(dressed, (k0, k1)):
            for g in model["G"]:
                np.testing.assert_allclose(d @ g, g @ d)
            actual = np.kron(d, I) @ source
            target = np.kron(w @ np.kron(k, I) @ old_gate(.37), I)
            np.testing.assert_allclose(actual @ rho @ actual.conj().T,
                                       target @ rho @ target.conj().T, atol=2e-14)
        np.testing.assert_allclose(sum(d.conj().T @ d for d in dressed),
                                   np.eye(128), atol=1e-14)

    def test_14_two_body_exact_Gauss_commutant_is_only_diagonal(self):
        for vertices, edges in ((2, ((0, 1),)),
                                (3, ((0, 1), (1, 2), (2, 0))), (4, EDGES)):
            examined, allowed = two_body_allowed(vertices, edges)
            total = vertices + len(edges)
            self.assertEqual(examined, 1 + 3 * total + 9 * total * (total - 1) // 2)
            self.assertEqual(len(allowed), 1 + total + total * (total - 1) // 2)
            self.assertTrue(all(not any(p in (1, 2) for p in row) for row in allowed))
        bare_bridge_flip = pauli(7, {6: X})
        self.assertAlmostEqual(np.linalg.norm(
            (np.eye(128) - gauge_model()["P"]) @ bare_bridge_flip
            @ blank_embedding()[:, 0]), 1.)

    def test_15_pulse_support_and_energy_norm_budget(self):
        gh, gd, gq = generators()["hamiltonians"]
        for h, expected in ((gh, 2), (gd, 1), (gq, 1)):
            self.assertAlmostEqual(max(abs(np.linalg.eigvalsh(h))), expected)
        independent_h = np.eye(128) - (
            pauli(7, {1: X, 2: X, 6: X}) + pauli(7, {6: Z})) / np.sqrt(2)
        independent_d = (np.eye(128) - pauli(7, {4: Z}) -
                         pauli(7, {6: Z}) + pauli(7, {4: Z, 6: Z})) / 4
        np.testing.assert_allclose(gh, independent_h)
        np.testing.assert_allclose(gd, independent_d)
        np.testing.assert_allclose(gq, pauli(7, {5: Z, 6: Z}))

    def test_16_full_history_Gauss_and_exact_invariant_sector(self):
        full, reduced = history(logical=False), history()
        model = gauge_model()
        for g in model["G"]:
            diagonal = np.tile(np.diag(g), 8)
            self.assertLess(np.max(np.abs(full["H"] *
                (diagonal[:, None] - diagonal[None, :]))), 1e-14)
        for physical_gate, logical_gate in zip(full["word"], reduced["word"]):
            np.testing.assert_allclose(physical_gate @ model["V"],
                                       model["V"] @ logical_gate, atol=1e-14)

    def test_17_history_unitary_conjugacy_norm_and_spectrum(self):
        model = history()
        transform = np.zeros((64, 64), complex)
        for j, u in enumerate(model["cumulative"]):
            transform[8*j:8*(j+1), 8*j:8*(j+1)] = u
        expected = transform @ np.kron(bare_clock(model["weights"]), np.eye(8)) \
                   @ transform.conj().T
        np.testing.assert_allclose(model["H"], expected, atol=1e-14)
        spectrum = np.linalg.eigvalsh(model["H"])
        np.testing.assert_allclose(spectrum, np.repeat(.5*np.arange(-3.5, 4, 1), 8),
                                   atol=2e-14)

    def test_18_exact_clock_arrival_preserves_all_inputs_and_reference(self):
        model = history()
        u = exponential(model["H"], 2*np.pi)
        incoming = clock_injection(0)
        target = clock_injection(7) @ model["cumulative"][-1]
        np.testing.assert_allclose(u @ incoming, (-1j)**7 * target, atol=2e-14)
        rho = random_density(16, 36518)
        actual = np.kron(u @ incoming, I)
        ideal = np.kron(target, I)
        np.testing.assert_allclose(actual @ rho @ actual.conj().T,
                                   ideal @ rho @ ideal.conj().T, atol=2e-14)

    def test_19_internal_clock_energy_and_Gauss_at_intermediate_times(self):
        model = history()
        data0 = kron_all([PLUS, PLUS, ZERO])
        initial = clock_injection(0) @ data0
        v = gauge_model()["V"]
        for t in (0., .19, np.pi, 2*np.pi):
            state = exponential(model["H"], t) @ initial
            hstate = model["H"] @ state
            self.assertAlmostEqual(np.vdot(state, hstate).real, 0.)
            self.assertAlmostEqual(np.vdot(hstate, hstate).real, .4375)
            physical_blocks = state.reshape(8, 8) @ v.T
            for g in gauge_model()["G"]:
                np.testing.assert_allclose(physical_blocks @ g.T,
                                           physical_blocks, atol=2e-14)

    def test_20_clock_readout_is_timed_and_returns_instead_of_halting(self):
        model = history()
        initial = clock_injection(0) @ kron_all([PLUS, PLUS, ZERO])
        for t in (np.pi, 2*np.pi + .1):
            state = exponential(model["H"], t) @ initial
            prob = np.linalg.norm(state.reshape(8, 8)[7])**2
            self.assertAlmostEqual(prob, np.sin(.5*t/2)**14)
        np.testing.assert_allclose(exponential(model["H"], 4*np.pi),
                                   -np.eye(64), atol=3e-14)
        self.assertAlmostEqual(capped_omega(7, 1.), .5)
        self.assertAlmostEqual(max(hop_weights(7, .5)), 1.)

    def test_21_two_body_obstruction_does_not_forbid_all_interaction(self):
        model = gauge_model()
        direct = pauli(7, {1: Z, 2: Z})
        np.testing.assert_allclose(direct @ model["V"],
            model["V"] @ pauli(3, {0: Z, 1: Z}))
        for g in model["G"]:
            np.testing.assert_allclose(direct @ g, g @ direct)
        theta = .37
        np.testing.assert_allclose(exponential(direct, theta) @ blank_embedding(),
                                   blank_embedding() @ old_gate(theta), atol=1e-14)


if __name__ == "__main__":
    main(__name__, "fixed_gauss_admission_audit", report)
