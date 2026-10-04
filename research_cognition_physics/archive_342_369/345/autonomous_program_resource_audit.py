"""Round 345: finite autonomous programming, exact no-go and finite accuracy.

The fixed-H constructions have an assumed initial program and sampling time.
They do not supply an autonomous clock, initial preparation, erasure, or a
physical energy cost from Hilbert-space dimension alone. The phase-grid error
is exact for nearest-single-catalogue selection, not optimal among processors.
"""
import math
import unittest
import numpy as np
from growing_stream_audit import main

I2 = np.eye(2, dtype=complex)
X = np.array([[0., 1.], [1., 0.]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)
HADAMARD = (X+Z)/np.sqrt(2)


def projector(psi):
    return np.outer(psi, np.conjugate(psi))


def evolve(hamiltonian, time=1.):
    values, vectors = np.linalg.eigh(hamiltonian)
    return (vectors*np.exp(-1j*time*values)) @ vectors.conj().T


def block_diagonal(blocks):
    d = blocks[0].shape[0]
    result = np.zeros((len(blocks)*d, len(blocks)*d), dtype=complex)
    for j, block in enumerate(blocks):
        result[j*d:(j+1)*d, j*d:(j+1)*d] = block
    return result


def exact_catalogue():
    """Fixed noncommuting catalogue, with one common positive energy shift."""
    unitaries = [I2, X, HADAMARD, np.diag([1., 1j])]
    logarithms = [np.zeros((2, 2)), np.pi*(I2-X)/2,
                  np.pi*(I2-HADAMARD)/2, np.diag([0., -np.pi/2])]
    shift = np.pi/2
    hamiltonian = block_diagonal([h+shift*I2 for h in logarithms])
    return unitaries, hamiltonian, np.exp(-1j*shift)


def phase_unitary(theta):
    return np.diag(np.exp(-.5j*theta*np.array([1., -1.])))


def phase_catalogue(k):
    angles = 2*np.pi*np.arange(k)/k
    hamiltonian = block_diagonal([np.pi*I2+.5*theta*Z for theta in angles])
    return angles, hamiltonian


def wrap_angle(theta):
    return (theta+np.pi) % (2*np.pi)-np.pi


def nearest_phase(theta, k):
    index = int(np.floor((theta % (2*np.pi))*k/(2*np.pi)+.5)) % k
    selected = 2*np.pi*index/k
    return index, selected, float(wrap_angle(theta-selected))


def phase_distance(theta, phi):
    """Half diamond distance, proved analytically for the qubit phase family."""
    return float(abs(np.sin(wrap_angle(theta-phi)/2)))


def trace_distance(rho, sigma):
    difference = rho-sigma
    difference = (difference+difference.conj().T)/2
    return float(.5*np.sum(abs(np.linalg.eigvalsh(difference))))


def random_density(d, rng, pure=False):
    if pure:
        psi = rng.normal(size=d)+1j*rng.normal(size=d)
        return projector(psi/np.linalg.norm(psi))
    mat = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
    rho = mat @ mat.conj().T
    return rho/np.trace(rho)


def apply_data_gate(unitary, rho_dr, reference_dim=2):
    joint = np.kron(unitary, np.eye(reference_dim))
    return joint @ rho_dr @ joint.conj().T


def program_output(rho, program_dim, remainder_dim):
    array = rho.reshape(program_dim, remainder_dim, program_dim, remainder_dim)
    return np.einsum('aibi->ab', array)


def discard_program(rho, program_dim, remainder_dim):
    array = rho.reshape(program_dim, remainder_dim, program_dim, remainder_dim)
    return np.einsum('aiaj->ij', array)


def inner_product_obstruction(u, v, overlap):
    """Minimum over complex output overlap b of ||a I - b U^dagger V||_F."""
    relative = u.conj().T @ v
    target = overlap*np.eye(u.shape[0])
    b = np.vdot(relative, target)/np.vdot(relative, relative)
    return float(np.linalg.norm(target-b*relative)), complex(b)


def resource_row(k):
    return {'catalogue_size': k, 'program_dimension': k,
            'qubit_embedding_size': int(math.ceil(math.log2(k))),
            'processor_dimension_P_times_D': 2*k,
            'including_idle_qubit_reference_dimension': 4*k,
            'nearest_single_gate_worst_half_diamond': float(np.sin(np.pi/(2*k))),
            'midpoint_adjacent_mixture_half_diamond': float(np.sin(np.pi/(2*k))**2),
            'fixed_H_spectral_lower_bound': 0.,
            'fixed_H_spectral_upper_bound': float(2*np.pi)}


def sufficient_catalogue_size(epsilon):
    if not 0 < epsilon <= 1:
        raise ValueError('epsilon must lie in (0,1].')
    # Guarantee applies only to nearest-single-unitary grid selection.
    return max(1, int(math.ceil(np.pi/(2*np.arcsin(epsilon)))))


def report():
    unitaries, hamiltonian, phase = exact_catalogue()
    w = evolve(hamiltonian)
    rng = np.random.default_rng(345)
    reference_error = 0.
    program_error = 0.
    energy_error = 0.
    for j, u in enumerate(unitaries):
        p = projector(np.eye(4)[j])
        rho = random_density(4, rng)
        initial = np.kron(p, rho)
        wr = np.kron(w, I2)
        out = wr @ initial @ wr.conj().T
        target = np.kron(p, apply_data_gate(u, rho))
        reference_error = max(reference_error, trace_distance(out, target))
        program_error = max(program_error, trace_distance(program_output(out, 4, 4), p))
    initial = random_density(8, rng)
    for time in (0., .23, .71, 1., 2.3):
        wt = evolve(hamiltonian, time)
        out = wt @ initial @ wt.conj().T
        energy_error = max(energy_error, abs(np.trace(hamiltonian @ (out-initial))))
    obstruction, b = inner_product_obstruction(I2, Z, 1/np.sqrt(2))
    eps_rows = []
    for epsilon in (.1, .01, .001):
        k = sufficient_catalogue_size(epsilon)
        eps_rows.append({'requested_half_diamond': epsilon,
                         'sufficient_K_for_this_scheme': k,
                         'achieved_worst_half_diamond': float(np.sin(np.pi/(2*k))),
                         'program_qubit_embedding': int(math.ceil(math.log2(k)))})
    return {'round': 345,
            'scope': 'Specified fixed finite-dimensional deterministic processor at one common time; the exact no-programming theorem is applied conditionally. Finite catalogue and approximate phase-grid constructions keep all registers and the fixed H, but assume initial preparation and time/readout. No cognitive derivation of spacetime, autonomous clock or thermodynamic cost is claimed.',
            'baseline': 'Frozen round 280 and phase-II operation results, with research through 341.',
            'question': 'Does resource closure permit a single fixed finite-dimensional autonomous processor that exactly and deterministically implements every unitary at one fixed time?',
            'answer': 'No under those added quantifiers. Resource closure itself does not impose fixed finite hardware, exactness, deterministic success, or a universal common completion time.',
            'primary_source': 'https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.79.321/fulltext',
            'no_programming': {'program_independent_of_unknown_data': True,
                               'hardware_and_completion_time_fixed': True,
                               'unitaries_compared_modulo_global_phase': True,
                               'orthogonal_programs_needed_for_distinct_exact_unitaries': True,
                               'exact_catalogue_dimension_lower_bound': 'dim(P) >= number of distinct unitary channels',
                               'program_return_required_for_the_no_go': False,
                               'I_vs_Z_attempt_input_overlap': float(1/np.sqrt(2)),
                               'least_possible_inner_product_identity_Frobenius_residual': obstruction,
                               'best_output_overlap_real': b.real,
                               'best_output_overlap_imag': b.imag},
            'autonomous_finite_catalogue': {'gates': ['I', 'X', 'Hadamard', 'phase S'],
                                           'fixed_time': 1., 'program_dimension': 4,
                                           'data_dimension': 2, 'idle_reference_dimension': 2,
                                           'processor_dimension': 8,
                                           'whole_P_D_R_dimension': 16,
                                           'fresh_environment_registers_for_this_unitary_catalogue': 0,
                                           'H_spectrum': np.linalg.eigvalsh(hamiltonian).tolist(),
                                           'fixed_H_evolution_vs_controlled_gate_max_error': float(
                                               np.max(abs(w-phase*block_diagonal(unitaries)))),
                                           'correlated_reference_output_max_trace_error': reference_error,
                                           'returned_program_max_trace_error': program_error,
                                           'closed_total_H_energy_max_drift': float(energy_error)},
            'phase_grid_resources': [resource_row(k) for k in (2, 4, 8, 16, 64)],
            'accuracy_examples': eps_rows,
            'exact_error_convention': 'epsilon=one-half diamond norm; nearest grid angle error <= pi/K, hence epsilon_K=sin(pi/(2K)). The ordinary diamond norm is twice this.',
            'not_an_optimality_claim': 'At a cell midpoint, an equal mixture of its two endpoint channels already has half-diamond error sin(pi/(2K))**2. Approximate quantum program lower bounds are not proved here.',
            'closed_resource_account': ['P, D, any correlated idle R, and the full fixed H are explicit.',
                                        'Basis programs return and may repeat the same catalogue gate at integer sampling times. Classical mixtures retain their program marginal, but can become correlated with data; coherent programs need not return.',
                                        'Preparing/selecting a program, compiling H, specifying time units and realizing a clock/readout remain inputs.',
                                        'No switching apparatus or replacement of H occurs during the modeled evolution.',
                                        'Changing target programs or performing measurement/reset requires its own resources; it is not supplied by this reusable-gate example.',
                                        'Hilbert dimension and formal H conservation do not determine thermodynamic preparation or maintenance costs.'],
            'scope_for_GR': 'The result blocks an over-strong finite-universe control interpretation; it neither derives spacetime/Einstein dynamics nor falsifies the original abstract operation contract.',
            'next_interface': 'Specify whether cognitive closure requires finite growing controllers, bounded errors and explicit clocks, before using implementability as a selector of spacetime or gravity.'}


class Checks(unittest.TestCase):
    def test_01_catalogue_gates_are_unitary_and_distinct_up_to_phase(self):
        gates, _, _ = exact_catalogue()
        for u in gates:
            np.testing.assert_allclose(u.conj().T @ u, I2, atol=1e-14)
        for j, u in enumerate(gates):
            for v in gates[j+1:]:
                self.assertLess(abs(np.trace(u.conj().T @ v)), 2.-1e-10)

    def test_02_one_positive_time_independent_H_realizes_the_whole_catalogue(self):
        gates, h, phase = exact_catalogue()
        np.testing.assert_allclose(h, h.conj().T, atol=1e-14)
        self.assertGreaterEqual(np.min(np.linalg.eigvalsh(h)), -1e-14)
        np.testing.assert_allclose(evolve(h), phase*block_diagonal(gates), atol=2e-14)
        np.testing.assert_allclose(evolve(h).conj().T @ evolve(h), np.eye(8), atol=2e-14)

    def test_03_unknown_data_and_correlated_reference_are_preserved_as_the_target_gate_requires(self):
        rng = np.random.default_rng(3345)
        gates, h, _ = exact_catalogue()
        for rdim in (1, 2, 3):
            w = np.kron(evolve(h), np.eye(rdim))
            for j, u in enumerate(gates):
                p = projector(np.eye(4)[j])
                for pure in (False, True):
                    rho = random_density(2*rdim, rng, pure)
                    initial = np.kron(p, rho)
                    out = w @ initial @ w.conj().T
                    expected = np.kron(p, apply_data_gate(u, rho, rdim))
                    np.testing.assert_allclose(out, expected, atol=2e-14)
                    np.testing.assert_allclose(program_output(out, 4, 2*rdim), p, atol=2e-14)

    def test_04_reusable_basis_program_gives_the_same_gate_powers_without_reset(self):
        gates, h, phase = exact_catalogue()
        for n in (2, 3, 7):
            target = phase**n*block_diagonal([np.linalg.matrix_power(u, n) for u in gates])
            np.testing.assert_allclose(evolve(h, n), target, atol=3e-14)

    def test_05_all_program_projectors_and_total_H_energy_are_conserved(self):
        gates, h, _ = exact_catalogue()
        for j in range(len(gates)):
            p = np.kron(projector(np.eye(4)[j]), I2)
            np.testing.assert_allclose(h @ p-p @ h, 0., atol=1e-14)
        rho = random_density(8, np.random.default_rng(5345))
        for time in np.linspace(0., 3., 17):
            w = evolve(h, time)
            self.assertLess(abs(np.trace(h @ (w @ rho @ w.conj().T-rho))), 1e-13)

    def test_06_nonorthogonal_exact_programs_conflict_with_inner_product_preservation(self):
        residual, best_b = inner_product_obstruction(I2, Z, 1/np.sqrt(2))
        self.assertAlmostEqual(residual, 1., places=14)
        self.assertAlmostEqual(abs(best_b), 0., places=14)
        for u, v in ((I2, X), (I2, phase_unitary(.8)), (X, HADAMARD)):
            self.assertGreater(inner_product_obstruction(u, v, .5)[0], .1)
            self.assertEqual(inner_product_obstruction(u, v, 0.)[0], 0.)
        residual, _ = inner_product_obstruction(I2, np.exp(.7j)*I2, .5)
        self.assertLess(residual, 1e-14)

    def test_07_global_input_and_output_Gram_matrices_match_for_orthogonal_programs(self):
        gates, h, _ = exact_catalogue()
        rng = np.random.default_rng(7345)
        columns = []
        for j in range(4):
            for _ in range(3):
                psi = rng.normal(size=2)+1j*rng.normal(size=2)
                psi /= np.linalg.norm(psi)
                columns.append(np.kron(np.eye(4)[j], psi))
        initial = np.column_stack(columns)
        out = evolve(h) @ initial
        np.testing.assert_allclose(out.conj().T @ out, initial.conj().T @ initial, atol=2e-14)

    def test_08_half_diamond_phase_formula_has_a_saturating_input_and_reference_bound(self):
        plus = projector(np.array([1., 1.])/np.sqrt(2))
        rng = np.random.default_rng(8345)
        for delta in (0., .01, .7, 1.9, np.pi, 2*np.pi-.2):
            u = phase_unitary(delta)
            expected = phase_distance(delta, 0.)
            self.assertAlmostEqual(trace_distance(u @ plus @ u.conj().T, plus), expected, places=14)
            for pure in (False, True):
                for _ in range(8):
                    rho = random_density(4, rng, pure)
                    self.assertLessEqual(trace_distance(apply_data_gate(u, rho), rho), expected+2e-14)

    def test_09_nearest_single_grid_gate_has_the_exact_uniform_worst_error(self):
        for k in (2, 4, 8, 16, 64):
            bound = np.sin(np.pi/(2*k))
            for theta in np.linspace(-2*np.pi, 4*np.pi, 1001):
                _, chosen, delta = nearest_phase(theta, k)
                self.assertLessEqual(abs(delta), np.pi/k+2e-14)
                self.assertLessEqual(phase_distance(theta, chosen), bound+2e-14)
            for j in range(k):
                theta = 2*np.pi*(j+.5)/k
                _, chosen, _ = nearest_phase(theta, k)
                self.assertAlmostEqual(phase_distance(theta, chosen), bound, places=13)

    def test_10_phase_catalogue_size_and_energy_spectrum_account_are_explicit(self):
        for k in (2, 4, 8, 16):
            angles, h = phase_catalogue(k)
            np.testing.assert_allclose(evolve(h), -block_diagonal(
                [phase_unitary(theta) for theta in angles]), atol=2e-14)
            self.assertGreaterEqual(np.min(np.linalg.eigvalsh(h)), 0.)
            self.assertLess(np.max(np.linalg.eigvalsh(h)), 2*np.pi)
            self.assertEqual(h.shape, (2*k, 2*k))
        for epsilon in (.1, .01, .001):
            k = sufficient_catalogue_size(epsilon)
            self.assertLessEqual(np.sin(np.pi/(2*k)), epsilon)
            self.assertGreater(np.sin(np.pi/(2*(k-1))), epsilon)

    def test_11_coherent_program_is_not_a_free_returned_program_for_a_unitary_superposition(self):
        _, h, _ = exact_catalogue()
        p = np.array([1., 1., 0., 0.])/np.sqrt(2)
        psi = np.kron(p, np.array([1., 0.]))
        out = projector(evolve(h) @ psi)
        pout = program_output(out, 4, 2)
        dout = discard_program(out, 4, 2)
        np.testing.assert_allclose(dout, .5*I2, atol=2e-14)
        self.assertGreater(trace_distance(pout, projector(p)), .49)
        self.assertAlmostEqual(np.trace(out @ out).real, 1., places=14)
        # The complete stored system remains reversible, despite either marginal mixing.
        back = evolve(h).conj().T @ out @ evolve(h)
        np.testing.assert_allclose(back, projector(psi), atol=2e-14)

    def test_12_adjacent_gate_mixture_improves_midpoint_error_so_grid_is_not_a_global_lower_bound(self):
        plus = projector(np.array([1., 1.])/np.sqrt(2))
        bell = projector(np.array([1., 0., 0., 1.])/np.sqrt(2))
        for k in (2, 4, 8, 16, 64):
            theta = np.pi/k
            target = phase_unitary(theta)
            u0 = I2; u1 = phase_unitary(2*theta)
            expected = np.sin(theta/2)**2
            for rho, rdim in ((plus, 1), (bell, 2)):
                mixture = .5*(apply_data_gate(u0, rho, rdim)+apply_data_gate(u1, rho, rdim))
                distance = trace_distance(mixture, apply_data_gate(target, rho, rdim))
                self.assertAlmostEqual(distance, expected, places=14)
            self.assertLess(expected, np.sin(theta/2))


if __name__ == '__main__':
    main(__name__, 'autonomous_program_resource_audit', report)
