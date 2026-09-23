"""Round 398: tracial coupling strength vs prepared-sector communication.

A fixed-H controlled SWAP gives an explicit, dimension-dependent separation.
This tests an interpretation of round-397 weights, not the reconstruction
identity itself. No geometry, light-cone violation, or autonomous resource
supply is claimed. Control preparation is a real additional resource.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np

TARGET = Path(__file__).with_name('prepared_sector_geometry_audit_results.json')
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1, -1]).astype(complex)
PAULIS = (X, Y, Z)
SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], complex)
SINGLET = (np.eye(4)-SWAP)/2


def comm(a, b):
    return a@b-b@a


def hs(a):
    return float(np.linalg.norm(a)/np.sqrt(a.shape[0]))


def operators(m):
    d = 2**m
    enabled = np.zeros(d, complex)
    enabled[-1] = 1
    projector = np.outer(enabled, enabled.conj())
    h = np.kron(projector, SINGLET)
    # H is itself a projector: exponential at t=pi is I-2H.
    u = np.eye(4*d)-2*h
    embedding = np.kron(enabled[:, None], np.eye(4))
    return projector, h, u, embedding


def pair_strength(h, left, right):
    return float(np.sqrt(sum(hs(comm(comm(h, a), b))**2
                             for a in left for b in right)/64))


def direct_audit(m):
    projector, h, u, v = operators(m)
    d = projector.shape[0]
    left = [np.kron(np.eye(d), np.kron(p, I)) for p in PAULIS]
    right = [np.kron(np.eye(d), np.kron(I, p)) for p in PAULIS]
    c = pair_strength(h, left, right)
    values, vectors = np.linalg.eigh(h)
    exponential = (vectors*np.exp(-1j*np.pi*values))@vectors.conj().T
    witness = comm(u.conj().T@left[2]@u, right[0])
    return dict(control_qubits=m, control_dimension=d,
                measured_coupling_strength=c, predicted_coupling_strength=np.sqrt(3)/(4*np.sqrt(d)),
                hamiltonian_operator_norm=float(np.max(np.abs(values))),
                exponential_error=hs(exponential-u), unitary_error=hs(u.conj().T@u-np.eye(4*d)),
                code_intertwining_error=float(np.linalg.norm(u@v-v@SWAP)),
                prepared_control_isometry_error=float(np.linalg.norm(v.conj().T@v-np.eye(4))),
                communication_commutator_operator_norm=float(np.linalg.norm(witness, ord=2)),
                communication_commutator_hs_norm=hs(witness),
                predicted_commutator_hs_norm=2/np.sqrt(d))


def reference_audit():
    rng = np.random.default_rng(398)
    error = coherence_error = 0.
    for m in (1, 2, 3, 4):
        projector, h, u, v = operators(m)
        # Arbitrary unknown AB plus a three-dimensional internal reference R.
        state = rng.normal(size=(4, 3))+1j*rng.normal(size=(4, 3))
        state /= np.linalg.norm(state)
        actual = u@v@state
        expected = v@SWAP@state
        error = max(error, float(np.linalg.norm(actual-expected)))
        # Coherent control input is retained in the global unitary description.
        all_state = rng.normal(size=(h.shape[0], 3))+1j*rng.normal(size=(h.shape[0], 3))
        all_state /= np.linalg.norm(all_state)
        coherence_error = max(coherence_error, float(np.linalg.norm(u.conj().T@u@all_state-all_state)))
    return dict(full_unknown_ab_reference_error=error,
                coherent_control_and_reference_inverse_error=coherence_error,
                controls_return_unchanged_on_enabled_preparation=True,
                unknown_controls_assumed_automatically_enabled=False)


def unprepared_audit(m):
    projector, h, u, v = operators(m)
    d = 2**m
    maximum = 0.
    for i, j in itertools.product(range(4), repeat=2):
        matrix = np.zeros((4, 4), complex)
        matrix[i, j] = 1
        evolved = u@np.kron(np.eye(d)/d, matrix)@u.conj().T
        actual = np.einsum('cacb->ab', evolved.reshape(d, 4, d, 4))
        expected = SWAP@matrix@SWAP/d+(1-1/d)*matrix
        maximum = max(maximum, float(np.linalg.norm(actual-expected)))
    # Bell input RA, blank B. Channel A->B after discarding A,C.
    bell = np.array([1, 0, 0, 1], complex)/np.sqrt(2)
    bell_density = np.outer(bell, bell.conj())
    replacement = np.kron(I/2, np.diag([1., 0.]))
    output = bell_density/d+(1-1/d)*replacement
    fidelity = float(np.vdot(bell, output@bell).real)
    return dict(control_qubits=m, active_probability=1/d, complete_ab_channel_error=maximum,
                bell_transfer_fidelity=fidelity, predicted_bell_fidelity=(1+3/d)/4,
                exact_half_diamond_error_to_transfer=1-1/d,
                guarantee_for_all_control_states=False)


def spectator_audit():
    m = 2
    projector, h, u, v = operators(m)
    d = 2**m
    left = [np.kron(np.eye(d), np.kron(p, I)) for p in PAULIS]
    right = [np.kron(np.eye(d), np.kron(I, p)) for p in PAULIS]
    c = pair_strength(h, left, right)
    spectator = np.eye(4)
    padded = pair_strength(np.kron(h, spectator), [np.kron(p, spectator) for p in left],
                           [np.kron(p, spectator) for p in right])
    logical = pair_strength(v.conj().T@h@v, [np.kron(p, I) for p in PAULIS],
                            [np.kron(I, p) for p in PAULIS])
    return dict(ambient_strength=c, spectator_padded_strength=padded,
                prepared_code_strength=logical, expected_code_ratio=np.sqrt(d),
                pure_relabelling_or_idle_padding_causes_dilution=False,
                restriction_to_enabled_sector_is_additional_preparation=True)


def conversion_audit():
    # H_{AB-support} with two canonical supports; operator and tracial norms.
    h01 = .4*np.kron(np.kron(X, X), I)+.3*np.kron(np.kron(Y, Z), I)
    h012 = -.2*np.kron(np.kron(Z, X), Y)+.6*np.kron(np.kron(X, Y), Z)
    supports = [h01, h012]
    h = sum(supports)
    left = [np.kron(np.kron(p, I), I) for p in PAULIS]
    right = [np.kron(np.kron(I, p), I) for p in PAULIS]
    c = pair_strength(h, left, right)
    weights = [hs(s)**2 for s in supports]
    norm_sum = sum(float(np.linalg.norm(s, ord=2)) for s in supports)
    k, number = 3, len(supports)
    bound = 2**(k/2)*np.sqrt(number)*c
    return dict(maximum_support_size=k, co_support_count=number,
                coupling_strength=c, squared_strength_from_supports=sum(weights),
                individual_operator_norm_sum=norm_sum, conversion_bound=float(bound),
                pair_part_operator_norm=float(np.linalg.norm(h, ord=2)),
                support_and_overlap_bounds_are_extra_inputs=True)


@lru_cache(None)
def report():
    return dict(round=398,
                scope='A controlled-SWAP family refutes a dimension-independent conversion from round-397 tracial coupling strength to all-allowed-preparation communication suppression. Fixed positive bounded H; preparation, joint support and timing are inputs. Not a spacetime model or a violation of a Lieb-Robinson theorem.',
                direct_cases=[direct_audit(m) for m in range(5)],
                unknown_reference=reference_audit(),
                mixed_control_cases=[unprepared_audit(m) for m in (1, 2, 3)],
                padding_and_code=spectator_audit(), support_conversion=conversion_audit(),
                analytic_scale_table=[dict(control_qubits=m,
                    ambient_coupling_strength=float(np.sqrt(3)/(4*2**(m/2))),
                    prepared_transfer_time=float(np.pi), hamiltonian_norm=1.,
                    enabled_control_qubits_required=m,
                    prepared_transfer_fidelity=1.) for m in (4, 8, 16, 32)],
                round397_exact_support_identity_invalidated=False,
                background_independent_transfer_claimed=False,
                preparation_cost_ignored=False, full_position_generation_completed=False)


class Audit(unittest.TestCase):
    def test_ambient_strength_and_fixed_positive_hamiltonian(self):
        for case in report()['direct_cases']:
            self.assertAlmostEqual(case['measured_coupling_strength'], case['predicted_coupling_strength'])
            self.assertAlmostEqual(case['hamiltonian_operator_norm'], 1.)
            self.assertLess(case['exponential_error'], 2e-12)

    def test_exact_code_action_and_unknown_complex_reference(self):
        for case in report()['direct_cases']:
            for key in ('unitary_error', 'code_intertwining_error', 'prepared_control_isometry_error'):
                self.assertLess(case[key], 2e-12)
        for key in ('full_unknown_ab_reference_error', 'coherent_control_and_reference_inverse_error'):
            self.assertLess(report()['unknown_reference'][key], 2e-12)

    def test_uniform_and_tracial_communication_norms_separate(self):
        for case in report()['direct_cases']:
            self.assertAlmostEqual(case['communication_commutator_operator_norm'], 2.)
            self.assertAlmostEqual(case['communication_commutator_hs_norm'], case['predicted_commutator_hs_norm'])

    def test_mixed_controls_do_not_give_background_independent_transfer(self):
        for case in report()['mixed_control_cases']:
            self.assertLess(case['complete_ab_channel_error'], 2e-12)
            self.assertAlmostEqual(case['bell_transfer_fidelity'], case['predicted_bell_fidelity'])
            self.assertFalse(case['guarantee_for_all_control_states'])

    def test_idle_padding_invariant_but_prepared_restriction_changes_weight(self):
        case = report()['padding_and_code']
        self.assertAlmostEqual(case['spectator_padded_strength'], case['ambient_strength'])
        self.assertAlmostEqual(case['prepared_code_strength']/case['ambient_strength'], case['expected_code_ratio'])

    def test_extra_support_and_overlap_bounds_control_operator_budget(self):
        case = report()['support_conversion']
        self.assertAlmostEqual(case['coupling_strength']**2, case['squared_strength_from_supports'])
        self.assertLessEqual(case['pair_part_operator_norm'], case['individual_operator_norm_sum']+1e-12)
        self.assertLessEqual(case['individual_operator_norm_sum'], case['conversion_bound']+1e-12)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checked = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checked.wasSuccessful():
        raise SystemExit(1)
    result = dict(report())
    result['checks'] = dict(run=checked.testsRun, failures=len(checked.failures), errors=len(checked.errors))
    result['runtime'] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
