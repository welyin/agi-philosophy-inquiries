"""Round 269: supplied graph edges do not amplify quantum communication for free.

The bound uses Schmidt number across a fixed operational cut. Classical
communication is allowed free of charge; quantum transmissions and initial
shared entanglement are explicitly counted. No physical clock is assumed.
"""
import argparse
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np


def maximally_entangled(d):
    return np.eye(d, dtype=complex).reshape(-1)/np.sqrt(d)


def choi(kraus):
    size = kraus[0].shape[1]
    phi = maximally_entangled(size)
    vectors = [np.kron(np.eye(size), k)@phi for k in kraus]
    return sum(np.outer(v, v.conj()) for v in vectors)


def fidelity(matrix):
    d = int(round(np.sqrt(len(matrix))))
    phi = maximally_entangled(d)
    return float(np.real(phi.conj()@matrix@phi))


def block_protocol(d, q):
    """d-dimensional unknown input; q communicated qubits plus a free block flag."""
    width = min(2**q, d)
    if d % width:
        raise ValueError('Use equal blocks in this calibration.')
    encoders = [np.eye(d, dtype=complex)[start:start+width] for start in range(0, d, width)]
    return encoders, [a.conj().T@a for a in encoders]


def teleport_kraus(correct=True):
    eye = np.eye(2)
    x = np.array([[0, 1], [1, 0]])
    z = np.diag([1, -1])
    resource = eye/np.sqrt(2)
    result = []
    for pauli in (eye, x, z, x@z):
        bell = resource@pauli.T
        contracted = np.einsum('as,sb->ba', bell.conj(), resource)
        result.append(pauli.conj().T@contracted if correct else contracted)
    return result


def state_rank(psi, dimensions, alice_axes):
    alice_axes = tuple(alice_axes)
    bob_axes = tuple(i for i in range(len(dimensions)) if i not in alice_axes)
    matrix = psi.reshape(dimensions).transpose(alice_axes+bob_axes)
    a = int(np.prod([dimensions[i] for i in alice_axes], dtype=int))
    return int(np.linalg.matrix_rank(matrix.reshape(a, -1), tol=1e-10))


def report():
    rows = []
    for r in (1, 2, 3):
        for q in range(r+1):
            d = 2**r
            encoders, kraus = block_protocol(d, q)
            rows.append({'independent_input_qubits': r, 'quantum_transmissions': q,
                         'classical_flag_bits': r-q, 'initial_shared_ebits': 0,
                         'entanglement_fidelity': fidelity(choi(kraus)),
                         'proved_upper_bound': min(1, 2**(q-r)),
                         'classical_basis_inputs_preserved': True})
    return {'round': 269, 'block_measurement_protocols': rows,
            'two_inputs_one_quantum_transmission': {
                'best_entanglement_fidelity_with_free_classical_flag': fidelity(choi(block_protocol(4, 1)[1])),
                'two_quantum_transmissions_fidelity': fidelity(choi([np.eye(4)])),
                'one_transmission_plus_one_consumed_ebit_fidelity':
                fidelity(choi([np.kron(np.eye(2), k) for k in teleport_kraus()]))},
            'resource_bound': 'Schmidt_number(output) <= R0 * product(d_j); F_e <= min(1, R0*product(d_j)/D). For q qubits, e initial Bell pairs, r target qubits: r <= q+e for exact transfer. Returning the Bell resource unchanged and uncorrelated restores r <= q.',
            'scope': 'A conditional communication-cut budget for a fixed bipartition with local ancillas, LOCC and counted noiseless quantum transmissions. It distinguishes additional primitive channels, serial reuse and consumed entanglement; it does not derive spacetime, an energy cost, a physical rate or the necessity of dimension three.'}


class Audit(unittest.TestCase):
    def test_01_complete_channels_and_reference_marginals(self):
        for d in (2, 4, 8):
            for q in range(int(np.log2(d))+1):
                _, kraus = block_protocol(d, q)
                np.testing.assert_allclose(sum(k.conj().T@k for k in kraus), np.eye(d))
                j = choi(kraus)
                self.assertGreaterEqual(np.linalg.eigvalsh(j)[0], -1e-12)
                np.testing.assert_allclose(np.trace(j.reshape(d, d, d, d), axis1=1, axis2=3),
                                           np.eye(d)/d, atol=1e-13)

    def test_02_free_classical_flag_protocol_attains_bound(self):
        for row in report()['block_measurement_protocols']:
            self.assertAlmostEqual(row['entanglement_fidelity'], row['proved_upper_bound'], places=12)

    def test_03_basis_success_does_not_certify_unknown_states(self):
        _, kraus = block_protocol(4, 1)
        for column in np.eye(4):
            rho = np.outer(column, column)
            out = sum(k@rho@k.conj().T for k in kraus)
            np.testing.assert_allclose(out, rho)
        psi = np.array([1, 0, 1, 0])/np.sqrt(2)
        rho = np.outer(psi, psi)
        out = sum(k@rho@k.conj().T for k in kraus)
        self.assertAlmostEqual(float(np.real(psi@out@psi)), 0.5)

    def test_04_low_schmidt_rank_overlap_bound(self):
        rng = np.random.default_rng(269)
        for d in (2, 4, 8):
            phi = maximally_entangled(d)
            for rank in (1, d//2, d):
                for _ in range(16):
                    a = rng.normal(size=(d, rank))+1j*rng.normal(size=(d, rank))
                    b = rng.normal(size=(rank, d))+1j*rng.normal(size=(rank, d))
                    psi = (a@b).reshape(-1)
                    psi /= np.linalg.norm(psi)
                    self.assertLessEqual(abs(phi.conj()@psi)**2, rank/d+1e-12)

    def test_05_one_transmission_multiplies_rank_by_at_most_message_dimension(self):
        rng = np.random.default_rng(1269)
        for message_dim in (2, 3):
            for rank in (1, 2, 3):
                for _ in range(12):
                    a = rng.normal(size=(4*message_dim, rank))
                    b = rng.normal(size=(rank, 5))
                    psi = (a@b).reshape(-1)
                    before = state_rank(psi, (4, message_dim, 5), (0, 1))
                    after = state_rank(psi, (4, message_dim, 5), (0,))
                    self.assertLessEqual(after, message_dim*before)

    def test_06_serial_delivery_pays_two_uses(self):
        psi = maximally_entangled(4)
        # Axes: reference, first input qubit, second input qubit.
        self.assertEqual(state_rank(psi, (4, 2, 2), (0, 1, 2)), 1)
        self.assertEqual(state_rank(psi, (4, 2, 2), (0, 2)), 2)
        self.assertEqual(state_rank(psi, (4, 2, 2), (0,)), 4)
        self.assertAlmostEqual(fidelity(choi([np.eye(4)])), 1)

    def test_07_teleportation_requires_the_classical_correction(self):
        corrected, raw = teleport_kraus(), teleport_kraus(False)
        np.testing.assert_allclose(sum(k.conj().T@k for k in corrected), np.eye(2), atol=1e-14)
        self.assertAlmostEqual(fidelity(choi(corrected)), 1)
        np.testing.assert_allclose(choi(raw), np.eye(4)/4, atol=1e-14)

    def test_08_one_sent_qubit_and_one_consumed_pair_preserve_external_reference(self):
        channel = [np.kron(np.eye(2), k) for k in teleport_kraus()]
        np.testing.assert_allclose(choi(channel), choi([np.eye(4)]), atol=1e-14)
        rng = np.random.default_rng(2269)
        psi = rng.normal(size=20)+1j*rng.normal(size=20)
        psi /= np.linalg.norm(psi)
        outputs = [np.kron(np.eye(5), k)@psi for k in channel]
        rho = sum(np.outer(v, v.conj()) for v in outputs)
        np.testing.assert_allclose(rho, np.outer(psi, psi.conj()), atol=1e-14)

    def test_09_returned_resource_does_not_supply_extra_net_rank(self):
        for d, resource_rank in ((2, 2), (4, 2), (4, 4)):
            target = np.kron(maximally_entangled(d), maximally_entangled(resource_rank))
            rank = state_rank(target, (d, d, resource_rank, resource_rank), (0, 2))
            self.assertEqual(rank, d*resource_rank)
            self.assertGreater(rank, resource_rank*(d//2))

    def test_10_forgetting_block_flag_loses_even_basis_delivery(self):
        encoders, with_flag = block_protocol(4, 1)
        fixed_decoder = encoders[0].conj().T
        without_flag = [fixed_decoder@a for a in encoders]
        np.testing.assert_allclose(sum(k.conj().T@k for k in without_flag), np.eye(4))
        basis = np.array([0, 0, 1, 0], complex)
        rho = np.outer(basis, basis.conj())
        bad = sum(k@rho@k.conj().T for k in without_flag)
        self.assertAlmostEqual(float(np.real(basis.conj()@bad@basis)), 0)
        self.assertAlmostEqual(fidelity(choi(with_flag)), 0.5)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result['runtime'] = {'python': platform.python_version(), 'numpy': np.__version__}
    result['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    payload = json.dumps(result, ensure_ascii=False, indent=2)+'\n'
    target = Path(__file__).with_name('interface_channel_budget_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
