"""Round 280: operational permissions, conserved observables and full accounting.

Uses standard adjoint-channel and amplitude-damping constructions. The charge,
interaction and supply of fresh ancillary systems are explicit model inputs.
"""
import unittest
import numpy as np
from growing_stream_audit import main

IDENTITY = np.eye(2, dtype=complex)
NUMBER = np.diag([0., 1.]).astype(complex)
X = np.array([[0., 1.], [1., 0.]], dtype=complex)
Y = np.array([[0., -1j], [1j, 0.]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)


def adjoint(kraus, observable):
    return sum(k.conj().T @ observable @ k for k in kraus)


def reset_kraus(state):
    return [np.outer(state, v.conj()) for v in IDENTITY]


def exchange(theta):
    """Basis 00, 01, 10, 11; excitation flows from A into a blank E."""
    c, s = np.cos(theta), np.sin(theta)
    u = np.eye(4, dtype=complex)
    u[1:3, 1:3] = [[c, s], [-s, c]]
    return u


def dilation(theta):
    return exchange(theta)[:, [0, 2]]


def local_channel(theta):
    v = dilation(theta).reshape(2, 2, 2)
    return [v[:, e, :] for e in range(2)]


def fixed_space_dimension(kraus):
    basis = [IDENTITY, X, Y, Z]
    matrix = np.stack([(adjoint(kraus, q)-q).reshape(-1) for q in basis], axis=1)
    real_matrix = np.vstack([matrix.real, matrix.imag])
    return 4-int(np.linalg.matrix_rank(real_matrix, tol=1e-11))


def two_site_gate(psi, unitary, first, second, sites):
    order = [first, second]+[j for j in range(sites) if j not in (first, second)]
    tensor = psi.reshape([2]*sites).transpose(order).reshape(4, -1)
    return (unitary @ tensor).reshape([2]*sites).transpose(np.argsort(order)).reshape(-1)


def fresh_environment(theta, collisions):
    """Explicit finite closed circuit: each of m environment qubits used once."""
    sites = collisions+1
    psi = np.zeros(2**sites, dtype=complex)
    psi[2**collisions] = 1.
    for e in range(1, sites):
        psi = two_site_gate(psi, exchange(theta), 0, e, sites)
    probs = np.abs(psi.reshape([2]*sites))**2
    charges = [float(probs.sum(axis=tuple(j for j in range(sites) if j != i))[1])
               for i in range(sites)]
    return psi, charges


def report():
    theta = np.pi/6
    rows = []
    for m in (1, 2, 3, 6, 8):
        _, charges = fresh_environment(theta, m)
        reused = np.linalg.matrix_power(exchange(theta), m)[:, 2]
        rows.append({'collisions': m, 'fresh_ancillas': m,
                     'visible_charge_fresh': charges[0],
                     'stored_environment_charge': sum(charges[1:]),
                     'full_charge': sum(charges),
                     'visible_charge_reused_single_ancilla': float(abs(reused[2])**2)})
    total = np.kron(NUMBER, IDENTITY)+np.kron(IDENTITY, NUMBER)
    flip = np.kron(X, IDENTITY)
    return {'round': 280,
            'scope': 'Operational conservation audit and supplied charge-preserving quantum circuit; no selected natural Hamiltonian, physical charge, energy unit or geometry.',
            'reset_fixed_hermitian_dimension': fixed_space_dimension(reset_kraus(np.array([1., 0.]))),
            'amplitude_damping_fixed_hermitian_dimension': fixed_space_dimension(local_channel(theta)),
            'charge_change_operator_norm_for_local_flip': float(np.linalg.norm(flip.conj().T @ total @ flip-total, 2)),
            'collision_angle': float(theta), 'collision_cases': rows,
            'interpretation': 'All allowed channels have only scalar common conserved observables; a chosen joint dynamics can conserve nontrivial additive charge while the observed subsystem relaxes.'}


class Audit(unittest.TestCase):
    def test_01_reset_has_only_scalar_fixed_observables(self):
        rng = np.random.default_rng(280)
        for state in (np.array([1., 0.]), np.array([1., 1j])/np.sqrt(2)):
            k = reset_kraus(state)
            self.assertEqual(fixed_space_dimension(k), 1)
            for _ in range(8):
                a = rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2))
                q = a+a.conj().T
                np.testing.assert_allclose(adjoint(k, q), np.vdot(state, q @ state)*IDENTITY, atol=1e-13)

    def test_02_trace_preservation_does_not_preserve_number(self):
        for theta in (0.1, 0.6, 1.2):
            k = local_channel(theta)
            np.testing.assert_allclose(adjoint(k, IDENTITY), IDENTITY, atol=1e-14)
            np.testing.assert_allclose(adjoint(k, NUMBER), np.cos(theta)**2*NUMBER, atol=1e-14)
            self.assertEqual(fixed_space_dimension(k), 1)

    def test_03_allowed_unitaries_do_not_select_conserved_charge(self):
        q = np.kron(NUMBER, IDENTITY)+np.kron(IDENTITY, NUMBER)
        for theta in (0.2, 0.8, 1.5):
            u = exchange(theta)
            np.testing.assert_allclose(u.conj().T @ q @ u, q, atol=1e-14)
        flip = np.kron(X, IDENTITY)
        self.assertAlmostEqual(np.linalg.norm(flip.conj().T @ q @ flip-q, 2), 1.)

    def test_04_dilation_intertwines_full_charge_and_is_isometric(self):
        q = np.kron(NUMBER, IDENTITY)+np.kron(IDENTITY, NUMBER)
        for theta in (0., 0.3, 1., np.pi/2):
            v = dilation(theta)
            np.testing.assert_allclose(v.conj().T @ v, IDENTITY, atol=1e-14)
            np.testing.assert_allclose(q @ v, v @ NUMBER, atol=1e-14)

    def test_05_environment_balance_is_an_operator_identity(self):
        for theta in (0.2, 0.7, 1.4):
            v = dilation(theta)
            environment = v.conj().T @ np.kron(IDENTITY, NUMBER) @ v
            np.testing.assert_allclose(environment+adjoint(local_channel(theta), NUMBER), NUMBER, atol=1e-14)
            np.testing.assert_allclose(environment, np.sin(theta)**2*NUMBER, atol=1e-14)

    def test_06_entangled_reference_and_reduced_channel(self):
        # Test the complete channel on a reference, not only classical populations.
        bell = np.array([1., 0., 0., 1.])/np.sqrt(2)
        rho = np.outer(bell, bell)
        for theta in (0.3, 0.9):
            out = np.kron(IDENTITY, dilation(theta)) @ bell  # order R,A,E
            tensor = out.reshape(2, 2, 2)
            ra = np.einsum('rae,sbe->rasb', tensor, tensor.conj()).reshape(4, 4)
            from_kraus = sum(np.kron(IDENTITY, k) @ rho @ np.kron(IDENTITY, k).conj().T for k in local_channel(theta))
            np.testing.assert_allclose(ra, from_kraus, atol=1e-14)
            reference = np.einsum('rae,sae->rs', tensor, tensor.conj())
            np.testing.assert_allclose(reference, IDENTITY/2, atol=1e-14)
            recovered = np.kron(IDENTITY, dilation(theta).conj().T) @ out
            np.testing.assert_allclose(recovered, bell, atol=1e-14)

    def test_07_fresh_environment_keeps_all_charge(self):
        theta = np.pi/6
        for m in range(1, 9):
            psi, charges = fresh_environment(theta, m)
            self.assertAlmostEqual(np.vdot(psi, psi).real, 1.)
            self.assertAlmostEqual(charges[0], np.cos(theta)**(2*m))
            self.assertAlmostEqual(sum(charges), 1.)
            for j in range(1, m+1):
                self.assertAlmostEqual(charges[j], np.sin(theta)**2*np.cos(theta)**(2*(j-1)))

    def test_08_finite_reused_environment_has_recurrence(self):
        theta = np.pi/6
        for m in range(1, 13):
            out = np.linalg.matrix_power(exchange(theta), m)[:, 2]
            self.assertAlmostEqual(abs(out[2])**2, np.cos(m*theta)**2)
            self.assertAlmostEqual(abs(out[1])**2+abs(out[2])**2, 1.)
        self.assertGreater(np.cos(6*theta)**2-np.cos(theta)**12, 0.8)


if __name__ == '__main__':
    main(__name__, 'operational_charge_audit', report)
