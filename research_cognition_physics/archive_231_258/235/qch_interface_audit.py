"""Round 235: audit the extra compatibility conditions in quantum causal histories.

QCH conditions are inputs from Hawkins et al. (hep-th/0302111), not consequences
of attaching arbitrary CP maps to a poset. NumPy checks certify small examples;
the note gives the dimension and reference-system arguments in general.
"""
import argparse
import json
from pathlib import Path
import platform
import unittest
import numpy as np


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)


def units(d):
    for i in range(d):
        for j in range(d):
            a = np.zeros((d, d), dtype=complex)
            a[i, j] = 1
            yield a


def ptrace(a, dims, keep):
    dims = list(dims)
    keep = set(keep)
    a = a.reshape(dims + dims)
    for k in reversed(range(len(dims))):
        if k not in keep:
            a = np.trace(a, axis1=k, axis2=k + len(dims))
            dims.pop(k)
    d = int(np.prod(dims))
    return a.reshape(d, d)


def bell(d):
    return np.eye(d).reshape(-1).astype(complex) / np.sqrt(d)


def state(v):
    return np.outer(v, v.conj())


def choi(channel, din):
    es = list(units(din))
    return sum(np.kron(e, channel(e)) for e in es)


def channel_checks(channel, din, dout):
    j = choi(channel, din)
    return {
        'choi_min_eigenvalue': float(np.linalg.eigvalsh(j).min()),
        'trace_preserving_error': float(np.linalg.norm(ptrace(j, [din, dout], [0]) - np.eye(din)))
    }


def monogamy_certificate(d):
    # Coordinates R,Y,Z; u_i=Bell_RY |i>_Z, v_i=Bell_RZ |i>_Y.
    u = np.zeros((d**3, d), dtype=complex)
    v = np.zeros_like(u)
    for i in range(d):
        for j in range(d):
            u[(j*d+j)*d+i, i] = 1/np.sqrt(d)
            v[(j*d+i)*d+j, i] = 1/np.sqrt(d)
    h = u@u.conj().T + v@v.conj().T
    expected = np.sort([0.]*(d**3-2*d) + [1-1/d]*d + [1+1/d]*d)
    return {
        'd': d, 'spectrum': np.linalg.eigvalsh(h).tolist(),
        'spectral_formula_error': float(np.max(np.abs(np.linalg.eigvalsh(h)-expected))),
        'overlap_error': float(np.linalg.norm(u.conj().T@v - np.eye(d)/d)),
        'maximum_sum_reference_fidelities': float(np.linalg.eigvalsh(h).max()),
        'two_identity_marginals_require': 2.,
    }


def random_unitary(d, rng):
    q, _ = np.linalg.qr(rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d)))
    return q


def fork_certificate(m, n, multiplicity=1):
    rng = np.random.default_rng(235000 + 100*m + 10*n + multiplicity)
    d = m*n*multiplicity
    u = random_unitary(d, rng)
    def left(a):
        return u.conj().T @ np.kron(np.kron(a, np.eye(n)), np.eye(multiplicity)) @ u
    def right(b):
        return u.conj().T @ np.kron(np.kron(np.eye(m), b), np.eye(multiplicity)) @ u
    aa, bb = list(units(m)), list(units(n))
    comm = max(np.linalg.norm(left(a)@right(b)-right(b)@left(a)) for a in aa for b in bb)
    hom = max(np.linalg.norm(left(a@b)-left(a)@left(b)) for a in aa for b in aa)
    hom = max(hom, max(np.linalg.norm(right(a@b)-right(a)@right(b)) for a in bb for b in bb))
    out = {'input_dimension': d, 'branch_dimensions': [m, n],
           'unobserved_multiplicity': multiplicity,
           'max_commutator_error': float(comm), 'max_homomorphism_error': float(hom),
           'complete_pair_in_this_construction': multiplicity == 1}
    if multiplicity == 1:
        out['joint_channel'] = channel_checks(lambda a:u@a@u.conj().T, d, d)
        entangled = bell(d)
        evolved = np.kron(np.eye(d), u) @ entangled
        recovered = np.kron(np.eye(d), u.conj().T) @ evolved
        out['reference_recovery_infidelity'] = float(abs(1-abs(np.vdot(entangled, recovered))**2))
    return out


def open_channel_certificate():
    depolarize = lambda a: np.trace(a)*I2/2
    swap = np.zeros((4, 4), dtype=complex)
    for i in range(2):
        for j in range(2):
            swap[j*2+i, i*2+j] = 1
    full = np.kron(np.eye(2), swap)
    initial = np.kron(state(bell(2)), I2/2)  # R,S,E
    evolved = full@initial@full.conj().T
    recovered_re = ptrace(evolved, [2, 2, 2], [0, 2])
    return {'depolarizing_channel': channel_checks(depolarize, 2, 2),
            'Heisenberg_multiplicativity_defect_on_X_squared': float(
                np.linalg.norm(depolarize(X@X)-depolarize(X)@depolarize(X))),
            'swap_dilation_channel_error': float(max(np.linalg.norm(
                ptrace(swap@np.kron(a, I2/2)@swap.conj().T, [2, 2], [0])-depolarize(a))
                for a in units(2))),
            'reference_preserved_in_environment_error': float(np.linalg.norm(recovered_re-state(bell(2))))}


def classical_copy_certificate():
    v = np.zeros((4, 2), dtype=complex)
    v[0, 0] = v[3, 1] = 1
    channel = lambda a:v@a@v.conj().T
    plus = state(np.array([1, 1])/np.sqrt(2))
    marginal = lambda a:ptrace(channel(a), [2, 2], [0])
    # The isometry copies orthogonal records, not the arbitrary quantum state.
    return {'channel': channel_checks(channel, 2, 4),
            'diagonal_record_error': float(max(np.linalg.norm(marginal(np.diag([p, 1-p]))-np.diag([p, 1-p]))
                                               for p in np.linspace(0, 1, 11))),
            'plus_input_marginal_trace_distance': float(np.abs(np.linalg.eigvalsh(marginal(plus)-plus)).sum()/2),
            'Heisenberg_homomorphism_defect': float(np.linalg.norm(
                v.conj().T@np.kron(X@X, I2)@v -
                (v.conj().T@np.kron(X, I2)@v)@(v.conj().T@np.kron(X, I2)@v)))}


def composition_certificate():
    direct = lambda a:X@a@X
    return {'each_pair_map_CPTP': channel_checks(direct, 2, 2),
            'max_direct_vs_two_identity_steps_error': float(max(np.linalg.norm(direct(a)-a) for a in units(2))),
            'consistent_identity_assignment_error': 0.}


def marginal_certificate():
    plus, minus = bell(2), np.array([1, 0, 0, -1])/np.sqrt(2)
    rp, rm = state(plus), state(minus)
    chp, chm = lambda a:np.trace(a)*rp, lambda a:np.trace(a)*rm
    return {'plus_preparation': channel_checks(chp, 2, 4),
            'minus_preparation': channel_checks(chm, 2, 4),
            'single_output_marginal_difference': float(max(np.linalg.norm(
                ptrace(rp, [2, 2], [i])-ptrace(rm, [2, 2], [i])) for i in [0, 1])),
            'joint_output_trace_distance': float(np.abs(np.linalg.eigvalsh(rp-rm)).sum()/2),
            'XX_expectations': [float(np.trace(r@np.kron(X, X)).real) for r in [rp, rm]]}


def report():
    return {'round': 235, 'date': '2026-09-22',
            'hypothesis': 'Are matrix algebras and individually CP edge maps sufficient for QCH?',
            'double_identity_fork': {'each_edge': channel_checks(lambda a:a, 2, 2),
                                    'cross_branch_X_Z_commutator_norm': float(np.linalg.norm(X@Z-Z@X))},
            'reference_monogamy': [monogamy_certificate(d) for d in [2, 3, 4]],
            'compatible_factor_forks': [fork_certificate(2, 2), fork_certificate(2, 3), fork_certificate(2, 2, 2)],
            'open_vs_complete_interface': open_channel_certificate(),
            'classical_record_copy': classical_copy_certificate(),
            'serial_composition': composition_certificate(),
            'marginals_do_not_determine_correlations': marginal_certificate(),
            'conclusion': 'Use QCH only after checking joint extension, commuting images and composition. Closed complete pairs conserve total Hilbert dimension; open channels need their environment represented.',
            'explicit_inputs': ['Chosen finite event order', 'Assignment of systems to events',
                                'Which interfaces are closed and complete', 'QCH extension axioms if invoked'],
            'not_claimed': ['New no-broadcasting theorem', 'Arbitrary causal sets satisfy QCH',
                            'Event count is spacetime volume', 'Quantum gravity has been derived'],
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__}}


class QCHTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = report()

    def check_channel(self, c):
        self.assertGreater(c['choi_min_eigenvalue'], -3e-13)
        self.assertLess(c['trace_preserving_error'], 3e-13)

    def test_local_cp_does_not_ensure_commuting_images(self):
        c = self.r['double_identity_fork']
        self.check_channel(c['each_edge'])
        self.assertAlmostEqual(c['cross_branch_X_Z_commutator_norm'], 2*np.sqrt(2))

    def test_reference_spectrum_for_all_three_dimensions(self):
        for c in self.r['reference_monogamy']:
            self.assertLess(c['spectral_formula_error'], 1e-13)
            self.assertLess(c['overlap_error'], 1e-13)
            self.assertAlmostEqual(c['maximum_sum_reference_fidelities'], 1+1/c['d'])
            self.assertLess(c['maximum_sum_reference_fidelities'], 2)

    def test_compatible_factors_and_complete_pairs(self):
        for c in self.r['compatible_factor_forks']:
            self.assertLess(c['max_commutator_error'], 1e-13)
            self.assertLess(c['max_homomorphism_error'], 1e-13)
            if c['complete_pair_in_this_construction']:
                self.check_channel(c['joint_channel'])
                self.assertLess(c['reference_recovery_infidelity'], 1e-13)
            else:
                self.assertEqual(c['unobserved_multiplicity'], 2)

    def test_open_channel_needs_environment_for_closed_description(self):
        c = self.r['open_vs_complete_interface']
        self.check_channel(c['depolarizing_channel'])
        self.assertAlmostEqual(c['Heisenberg_multiplicativity_defect_on_X_squared'], np.sqrt(2))
        self.assertLess(c['swap_dilation_channel_error'], 1e-13)
        self.assertLess(c['reference_preserved_in_environment_error'], 1e-13)

    def test_classical_records_not_arbitrary_quantum_states(self):
        c = self.r['classical_record_copy']
        self.check_channel(c['channel'])
        self.assertLess(c['diagonal_record_error'], 1e-13)
        self.assertAlmostEqual(c['plus_input_marginal_trace_distance'], .5)
        self.assertAlmostEqual(c['Heisenberg_homomorphism_defect'], np.sqrt(2))

    def test_composition_is_an_additional_compatibility_requirement(self):
        c = self.r['serial_composition']
        self.check_channel(c['each_pair_map_CPTP'])
        self.assertAlmostEqual(c['max_direct_vs_two_identity_steps_error'], np.sqrt(2))
        self.assertEqual(c['consistent_identity_assignment_error'], 0.)

    def test_joint_correlations_are_retained(self):
        c = self.r['marginals_do_not_determine_correlations']
        self.check_channel(c['plus_preparation'])
        self.check_channel(c['minus_preparation'])
        self.assertLess(c['single_output_marginal_difference'], 1e-13)
        self.assertAlmostEqual(c['joint_output_trace_distance'], 1.)
        np.testing.assert_allclose(c['XX_expectations'], [1, -1], atol=1e-13)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(QCHTests))
    if not result.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': result.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        target = Path(__file__).with_name('qch_interface_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != json.loads(json.dumps(data)):
            raise RuntimeError('Existing result differs; preserve and inspect before replacing.')
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))
