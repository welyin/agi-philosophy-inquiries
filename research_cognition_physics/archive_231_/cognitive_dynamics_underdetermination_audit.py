"""Round 342: complete finite quantum contracts do not select natural dynamics.

All finite instruments are primitives of the inherited operational contract.
This is not an autonomous implementation of every control by one Hamiltonian.
"""
import itertools
import unittest
import numpy as np
from growing_stream_audit import main


def pure(vector):
    vector = np.asarray(vector, dtype=complex)
    vector = vector / np.linalg.norm(vector)
    return np.outer(vector, vector.conj())


def random_density(dimension, rng):
    matrix = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(
        size=(dimension, dimension))
    value = matrix @ matrix.conj().T
    return value / np.trace(value)


def spectral(matrix, function):
    values, vectors = np.linalg.eigh(matrix)
    return (vectors * function(values)) @ vectors.conj().T


def evolution(hamiltonian, time):
    return spectral(hamiltonian, lambda values: np.exp(-1j * time * values))


def channel(kraus, rho):
    return sum(operator @ rho @ operator.conj().T for operator in kraus)


def reduced(rho, dimensions, keep):
    remaining = list(range(len(dimensions)))
    tensor = rho.reshape(tuple(dimensions) * 2)
    for label in reversed(range(len(dimensions))):
        if label not in keep:
            position = remaining.index(label)
            tensor = np.trace(tensor, axis1=position,
                              axis2=position + len(remaining))
            remaining.remove(label)
    size = int(np.prod([dimensions[i] for i in remaining]))
    return tensor.reshape(size, size)


def trace_distance(left, right):
    return float(np.sum(np.abs(np.linalg.eigvalsh(left - right))) / 2)


def steering_resource(rho):
    root = spectral(rho, lambda values: np.sqrt(np.maximum(values, 0)))
    return pure(root.reshape(-1))


def steering_povm(rho, decomposition):
    values, vectors = np.linalg.eigh(rho)
    support = values > 1e-11
    inverse = (vectors[:, support] / np.sqrt(values[support])) @ (
        vectors[:, support].conj().T)
    projector = vectors[:, support] @ vectors[:, support].conj().T
    effects = [(inverse @ piece @ inverse).T for piece in decomposition]
    effects[0] = effects[0] + (np.eye(len(rho)) - projector).T
    return effects


def conditional_target(resource, effect):
    d = len(effect)
    # Independent density-tensor contraction; target comes before reference.
    return np.einsum("abcd,db->ac", resource.reshape(d, d, d, d), effect)


def decomposition_case(d, rank, count, rng):
    basis, _ = np.linalg.qr(rng.normal(size=(d, rank)) +
                            1j * rng.normal(size=(d, rank)))
    pieces = []
    for _ in range(count):
        matrix = rng.normal(size=(rank, rank)) + 1j * rng.normal(
            size=(rank, rank))
        pieces.append(basis @ matrix @ matrix.conj().T @ basis.conj().T)
    total = sum(float(np.trace(piece).real) for piece in pieces)
    pieces = [piece / total for piece in pieces]
    return sum(pieces), pieces


def steering_error(rho, pieces):
    resource = steering_resource(rho)
    effects = steering_povm(rho, pieces)
    return max(float(np.linalg.norm(conditional_target(resource, effect) - piece))
               for effect, piece in zip(effects, pieces))


def pure_path_generator(psi, phi):
    psi, phi = psi / np.linalg.norm(psi), phi / np.linalg.norm(phi)
    overlap = np.vdot(psi, phi)
    aligned = phi * np.exp(-1j * np.angle(overlap))
    cosine = min(1., abs(overlap))
    if 1 - cosine < 1e-14:
        return np.zeros((len(psi), len(psi)), dtype=complex)
    normal = (aligned - cosine * psi) / np.sqrt(1 - cosine**2)
    return 1j * np.arccos(cosine) * (
        np.outer(normal, psi.conj()) - np.outer(psi, normal.conj()))


def tomographic_states(d):
    basis = np.eye(d, dtype=complex)
    states = [pure(row) for row in basis]
    for i, j in itertools.combinations(range(d), 2):
        states += [pure(basis[i] + basis[j]),
                   pure(basis[i] + 1j * basis[j])]
    return states


def random_instrument(source, target, outcomes, rng):
    raw = rng.normal(size=(target * outcomes, source)) + 1j * rng.normal(
        size=(target * outcomes, source))
    isometry, _ = np.linalg.qr(raw)
    return [isometry[i * target:(i + 1) * target] for i in range(outcomes)]


def hamiltonian(dimensions, edges=(), coupling=True):
    labels = np.array(list(itertools.product(
        *(range(d) for d in dimensions))), dtype=float)
    energies = np.sum(labels, axis=1)
    if coupling:
        for first, second, weight in edges:
            energies += weight * labels[:, first] * labels[:, second]
    return np.diag(energies)


def pair_witness(coupling, time=np.pi):
    h = hamiltonian((2, 2), ((0, 1, 1.),), coupling)
    unitary = evolution(h, time)
    plus = pure([1, 1])
    outputs = []
    for bit in (0, 1):
        rho = np.kron(pure(np.eye(2)[bit]), plus)
        outputs.append(reduced(channel([unitary], rho), (2, 2), (1,)))
    entangled = channel([unitary], np.kron(plus, plus))
    local = reduced(entangled, (2, 2), (0,))
    return {"receiver_trace_distance": trace_distance(*outputs),
            "receiver_optimal_equal_prior_success":
                (1 + trace_distance(*outputs)) / 2,
            "product_input_final_local_purity":
                float(np.trace(local @ local).real),
            "energy_spectrum": np.diag(h).tolist()}


class Checks(unittest.TestCase):
    def test_01_full_rank_universal_steering(self):
        rng = np.random.default_rng(34201)
        for d in (2, 3, 5):
            for count in (1, 2, 7):
                rho, pieces = decomposition_case(d, d, count, rng)
                effects = steering_povm(rho, pieces)
                self.assertLess(steering_error(rho, pieces), 2e-12)
                self.assertLess(np.linalg.norm(sum(effects) - np.eye(d)), 3e-11)
                self.assertGreater(min(np.linalg.eigvalsh(e)[0] for e in effects),
                                   -2e-12)

    def test_02_rank_deficient_and_zero_branch_steering(self):
        rng = np.random.default_rng(34202)
        for d in (2, 3, 5):
            for rank in range(1, d):
                rho, pieces = decomposition_case(d, rank, 4, rng)
                pieces.append(np.zeros_like(rho))
                effects = steering_povm(rho, pieces)
                self.assertLess(steering_error(rho, pieces), 2e-12)
                self.assertLess(np.linalg.norm(sum(effects) - np.eye(d)), 3e-11)
                self.assertGreater(min(np.linalg.eigvalsh(e)[0] for e in effects),
                                   -2e-12)

    def test_03_same_resource_for_different_ensembles(self):
        rho = np.eye(2) / 2
        computational = [pure([1, 0]) / 2, pure([0, 1]) / 2]
        phase = [pure([1, 1j]) / 2, pure([1, -1j]) / 2]
        resource = steering_resource(rho)
        for pieces in (computational, phase):
            for effect, piece in zip(steering_povm(rho, pieces), pieces):
                self.assertLess(np.linalg.norm(
                    conditional_target(resource, effect) - piece), 1e-13)
        self.assertLess(np.linalg.norm(reduced(resource, (2, 2), (0,)) - rho),
                        1e-13)

    def test_04_connected_pure_transitivity(self):
        rng = np.random.default_rng(34204)
        for d in (2, 3, 6, 9):
            psi = rng.normal(size=d) + 1j * rng.normal(size=d)
            phi = rng.normal(size=d) + 1j * rng.normal(size=d)
            generator = pure_path_generator(psi, phi)
            unitary = evolution(generator, 1.)
            self.assertLess(np.linalg.norm(channel([unitary], pure(psi)) -
                                           pure(phi)), 2e-13)
            self.assertLess(np.linalg.norm(evolution(generator, .3) @
                evolution(generator, .7) - unitary), 2e-13)

    def test_05_local_tomographic_spanning_for_every_sampled_type(self):
        for d in (2, 3, 4):
            states = tomographic_states(d)
            self.assertEqual(len(states), d*d)
            self.assertEqual(np.linalg.matrix_rank(np.array(
                [rho.reshape(-1) for rho in states])), d*d)
        product = [np.kron(a, b).reshape(-1) for a in tomographic_states(2)
                   for b in tomographic_states(3)]
        self.assertEqual(np.linalg.matrix_rank(product), 36)

    def test_06_P_kraus_equivalence_with_entangled_reference(self):
        gamma = .37
        kraus = [np.diag([1, np.sqrt(1-gamma)]),
                 np.array([[0, np.sqrt(gamma)], [0, 0]])]
        mixed = [(kraus[0] + 1j*kraus[1])/np.sqrt(2),
                 (1j*kraus[0] + kraus[1])/np.sqrt(2)]
        bell = pure([1, 0, 0, 1])
        for rho in tomographic_states(2):
            self.assertLess(np.linalg.norm(channel(kraus, rho) -
                                           channel(mixed, rho)), 1e-13)
        self.assertLess(np.linalg.norm(channel(
            [np.kron(k, np.eye(2)) for k in kraus], bell) - channel(
            [np.kron(k, np.eye(2)) for k in mixed], bell)), 1e-13)

    def test_07_full_instrument_composition_and_reference(self):
        rng = np.random.default_rng(34207)
        first = random_instrument(2, 3, 3, rng)
        second = random_instrument(3, 2, 4, rng)
        rho = random_density(4, rng)
        composite = [b @ a for a in first for b in second]
        for operators, d in ((first, 2), (second, 3), (composite, 2)):
            self.assertLess(np.linalg.norm(sum(k.conj().T @ k for k in operators)
                                           - np.eye(d)), 2e-13)
        output = channel([np.kron(k, np.eye(2)) for k in composite], rho)
        self.assertGreater(np.linalg.eigvalsh(output)[0], -1e-13)
        self.assertLess(np.linalg.norm(reduced(output, (2, 2), (1,)) -
            reduced(rho, (2, 2), (1,))), 2e-13)

    def test_08_nontrivial_Time_group_in_both_models(self):
        for coupling in (False, True):
            h = hamiltonian((2, 3), ((0, 1, .7),), coupling)
            self.assertLess(np.linalg.norm(evolution(h, .2) @ evolution(h, .6) -
                                           evolution(h, .8)), 1e-13)
            probe = pure(np.ones(6))
            self.assertGreater(trace_distance(channel([evolution(h, .4)],
                                                       probe), probe), .1)

    def test_09_independent_composition_all_hamiltonian_families(self):
        for coupling in (False, True):
            first = hamiltonian((2, 3), ((0, 1, .7),), coupling)
            second = hamiltonian((2, 2), ((0, 1, .3),), coupling)
            joint = hamiltonian((2, 3, 2, 2),
                                ((0, 1, .7), (2, 3, .3)), coupling)
            self.assertLess(np.linalg.norm(joint -
                np.kron(first, np.eye(4)) - np.kron(np.eye(6), second)), 1e-13)
            self.assertLess(np.linalg.norm(evolution(joint, .4) -
                np.kron(evolution(first, .4), evolution(second, .4))), 3e-13)

    def test_10_extension_preserves_unknown_state_and_reference(self):
        rng = np.random.default_rng(34210)
        rho_ar = random_density(8, rng)
        sigma_b = random_density(3, rng)
        for coupling in (False, True):
            ua = evolution(hamiltonian((2, 2), ((0, 1, .7),), coupling), .31)
            ub = evolution(hamiltonian((3,), (), coupling), .31)
            output = channel([np.kron(np.kron(ua, np.eye(2)), ub)],
                             np.kron(rho_ar, sigma_b))
            expected = channel([np.kron(ua, np.eye(2))], rho_ar)
            self.assertLess(np.linalg.norm(reduced(output, (4, 2, 3), (0, 1)) -
                                           expected), 2e-13)

    def test_11_exact_natural_signaling_difference(self):
        self.assertLess(pair_witness(False)["receiver_trace_distance"], 1e-13)
        self.assertAlmostEqual(pair_witness(True)["receiver_trace_distance"],
                               1., places=13)

    def test_12_entanglement_generation_distinguishes_models(self):
        self.assertAlmostEqual(pair_witness(False)[
            "product_input_final_local_purity"], 1., places=13)
        self.assertAlmostEqual(pair_witness(True)[
            "product_input_final_local_purity"], .5, places=13)

    def test_13_product_time_no_signal_for_correlated_initial_states(self):
        rng = np.random.default_rng(34213)
        rho = random_density(12, rng)
        local = random_instrument(2, 2, 3, rng)
        changed = channel([np.kron(k, np.eye(6)) for k in local], rho)
        for time in (0., .7, 2.1):
            u = evolution(hamiltonian((2, 3, 2), (), False), time)
            outputs = [reduced(channel([u], value), (2, 3, 2), (1, 2))
                       for value in (rho, changed)]
            self.assertLess(trace_distance(*outputs), 3e-13)

    def test_14_block_local_coding_cannot_create_natural_influence(self):
        rng = np.random.default_rng(34214)
        # Two microscopic qubits per block; receiver has arbitrary local
        # encoding/decoding and a possibly initially correlated partner.
        rho = random_density(16, rng)
        ua = evolution(random_density(4, rng), .6)
        local_decoder = random_instrument(4, 3, 3, rng)
        changed = channel([np.kron(ua, np.eye(4))], rho)
        natural = evolution(hamiltonian((2, 2, 2, 2), (), False), .9)
        outputs = [channel(local_decoder,
            reduced(channel([natural], state), (4, 4), (1,)))
                   for state in (rho, changed)]
        self.assertLess(trace_distance(*outputs), 3e-13)

    def test_15_spectra_not_time_unit_or_energy_origin_change(self):
        left = np.linalg.eigvalsh(hamiltonian((2, 2), (), False))
        right = np.linalg.eigvalsh(hamiltonian((2, 2), ((0, 1, 1.),), True))
        # Positive rescaling and scalar shift must match the first two
        # distinct levels, fixing scale=1 and shift=0, but fail at the last.
        scale = (right[1] - right[0]) / (left[1] - left[0])
        shift = right[0] - scale * left[0]
        self.assertEqual(float(np.max(np.abs(right-scale*left-shift))), 1.)


def report():
    rng = np.random.default_rng(34242)
    cases = []
    for d, rank, count in ((2, 2, 4), (3, 2, 5), (5, 1, 7), (5, 5, 8)):
        rho, pieces = decomposition_case(d, rank, count, rng)
        effects = steering_povm(rho, pieces)
        cases.append({"dimension": d, "rank": rank, "outcomes": count,
                      "conditional_state_error": steering_error(rho, pieces),
                      "POVM_completeness_error":
                          float(np.linalg.norm(sum(effects)-np.eye(d))),
                      "resource_purity":
                          float(np.trace(steering_resource(rho) @
                                         steering_resource(rho)).real)})
    return {"round": 342,
            "scope": "All finite complex quantum types and instruments; two "
                     "composition-compatible nontrivial Time expansions. "
                     "Controls are primitives, not autonomously compiled.",
            "contract": ["F", "U_all_states_all_finite_ensembles_same_resource",
                         "C_on_every_composite", "P", "R_seed",
                         "nontrivial_Time", "independent_finite_composition"],
            "steering_checks": cases,
            "pair_time": float(np.pi),
            "product_natural_time": pair_witness(False),
            "interacting_natural_time": pair_witness(True),
            "analytic_claim": "The inherited formal contract does not entail "
                "nonzero natural influence across independent microscopic "
                "factors, even after fixed block-local coarse graining.",
            "not_proved": ["No-go for all SoCA semantics",
                          "No-go for controlled GR simulations",
                          "Autonomous realization of all controls by one "
                          "closed Hamiltonian and internal programs",
                          "No-go for every nonlocal or state-dependent "
                          "emergent encoding",
                          "Impossibility of deriving GR from strengthened "
                          "cognitive principles"]}


if __name__ == "__main__":
    main(__name__, "cognitive_dynamics_underdetermination_audit", report)
