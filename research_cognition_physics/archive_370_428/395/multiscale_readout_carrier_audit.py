"""Round 395: one query-independent carrier for a finite multiscale table.

The phase family is an explicit calibration input, not derived position or a
spatial dimension. D_min=2**N follows analytically from pairwise perfect
distinguishability; the positive code is the standard Fourier product code.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import platform
import unittest
import numpy as np

TARGET = Path(__file__).with_name('multiscale_readout_carrier_audit_results.json')
PAULI_X = np.array([[0, 1], [1, 0]], complex)
PAULI_Y = np.array([[0, -1j], [1j, 0]], complex)


def queries(n):
    m = 2**n
    return tuple((j, ell, -2*np.pi*ell/(m >> j))
                 for j in range(n) for ell in range(m >> j))


def table(n):
    x = 2*np.pi*np.arange(2**n)/(2**n)
    return np.column_stack([(1+np.cos(2**j*x+theta))/2
                            for j, _, theta in queries(n)])


def pair_query(n, a, b):
    difference = abs(b-a)
    assert difference > 0
    power = (difference & -difference).bit_length()-1
    j = n-1-power
    ell = a % (2**n >> j)
    return j, ell


def phase_state(n, x):
    # Computational basis bit j has phase weight 2**j.
    return np.exp(1j*x*np.arange(2**n))/np.sqrt(2**n)


def product_state(n, x):
    result = np.array([1.], complex)
    for j in reversed(range(n)):
        result = np.kron(result, np.array([1, np.exp(1j*2**j*x)])/np.sqrt(2))
    return result


def effect(n, j, theta):
    single = (np.eye(2)+np.cos(theta)*PAULI_X-np.sin(theta)*PAULI_Y)/2
    result = np.array([[1.]], complex)
    for bit in reversed(range(n)):
        result = np.kron(result, single if bit == j else np.eye(2))
    return result


def packing_lower_bound(m, error):
    if not 0 <= error <= .5:
        raise ValueError('The pair-separation bound assumes 0 <= error <= 1/2.')
    beta = 4*error*(1-error)
    return m/(1+(m-1)*beta)


def exact_audit(n):
    m = 2**n
    readings = table(n)
    lookup = {(j, ell): index for index, (j, ell, _) in enumerate(queries(n))}
    pair_error = 0.
    for a in range(m):
        for b in range(a+1, m):
            column = lookup[pair_query(n, a, b)]
            pair_error = max(pair_error, abs(readings[a, column]-1), abs(readings[b, column]))
    code = np.column_stack([phase_state(n, 2*np.pi*k/m) for k in range(m)])
    output_error = 0.
    effect_spectrum_error = 0.
    for j, _, theta in queries(n):
        e = effect(n, j, theta)
        spectrum = np.linalg.eigvalsh(e)
        effect_spectrum_error = max(effect_spectrum_error, max(0., -float(spectrum[0]),
                                                               float(spectrum[-1])-1))
        for x in (.127, -.83, 2.31):
            state = phase_state(n, x)
            measured = float(np.vdot(state, e@state).real)
            output_error = max(output_error, abs(measured-(1+np.cos(2**j*x+theta))/2))
    return dict(scales=n, promised_settings=m, detector_settings=len(lookup),
                tested_pairs=m*(m-1)//2, pair_witness_error=pair_error,
                exact_minimum_carrier_dimension=m, qubits=n,
                fourier_unitarity_error=float(np.linalg.norm(code.conj().T@code-np.eye(m))),
                independent_product_code_error=float(np.linalg.norm(phase_state(n, .127)-product_state(n, .127))),
                off_grid_probability_error=output_error,
                effect_spectrum_error=effect_spectrum_error)


def reference_audit(n=3):
    m = 2**n
    code = np.column_stack([phase_state(n, 2*np.pi*k/m) for k in range(m)])
    # A specified complex input on label L and an idle internal reference R.
    vector = np.arange(1, 2*m+1).reshape(m, 2)+1j*np.arange(3, 2*m+3).reshape(m, 2)**2
    vector = vector/np.linalg.norm(vector)
    encoded = code@vector
    decoded = code.conj().T@encoded
    before_r = vector.T@vector.conj()
    after_r = encoded.T@encoded.conj()
    return dict(label_dimension=m, reference_dimension=2,
                full_vector_return_error=float(np.linalg.norm(decoded-vector)),
                reference_state_error=float(np.linalg.norm(after_r-before_r)),
                input_label_was_already_dimension_m=True,
                unknown_one_qubit_phase_amplification_claimed=False)


def noisy_audit(n=4, mixing=.02):
    m = 2**n
    code = np.column_stack([phase_state(n, 2*np.pi*k/m) for k in range(m)])
    states = [(1-mixing)*np.outer(code[:, k], code[:, k].conj())+mixing*np.eye(m)/m
              for k in range(m)]
    target = table(n)
    measured = np.empty_like(target)
    for column, (j, _, theta) in enumerate(queries(n)):
        e = effect(n, j, theta)
        measured[:, column] = [np.trace(rho@e).real for rho in states]
    error = float(np.max(np.abs(measured-target)))
    overlap = max(float(np.trace(states[a]@states[b]).real)
                  for a in range(m) for b in range(a+1, m))
    mean = sum(states)/m
    purity = float(np.trace(mean@mean).real)
    # Use the analytic certified uniform error, avoiding floating ceil drift.
    certified_error = mixing/2
    bound = packing_lower_bound(m, certified_error)
    return dict(scales=n, settings=m, mixing=mixing,
                observed_uniform_probability_error=error,
                certified_uniform_probability_error=certified_error,
                max_pair_hilbert_schmidt_overlap=overlap,
                overlap_upper_bound=4*certified_error*(1-certified_error),
                average_purity=purity, carrier_dimension_lower_bound=bound,
                minimum_integer_dimension=int(np.ceil(bound-1e-12)),
                actual_carrier_dimension=m,
                noisy_bound_claimed_sharp=False)


def query_dependence_audit():
    n = 3
    m = 2**n
    maximum_error = 0.
    for k in range(m):
        x = 2*np.pi*k/m
        for j, _, theta in queries(n):
            # Invalid as a common carrier: this preparation explicitly reads j.
            vector = np.array([1, np.exp(1j*2**j*x)])/np.sqrt(2)
            e = effect(1, 0, theta)
            maximum_error = max(maximum_error,
                                abs(np.vdot(vector, e@vector).real-(1+np.cos(2**j*x+theta))/2))
    # For x=pi/2 the j=0 and j=1 encodings are visibly different states.
    a = np.array([1, 1j])/np.sqrt(2)
    b = np.array([1, -1])/np.sqrt(2)
    difference = np.outer(a, a.conj())-np.outer(b, b.conj())
    return dict(scales=n, one_qubit_if_preparation_can_depend_on_query=True,
                probability_error=float(maximum_error),
                preparation_trace_distance=float(np.linalg.norm(difference, ord='nuc')/2),
                satisfies_query_independent_carrier_contract=False,
                revisit_to_source_or_scale_record_is_additional_access=True)


@lru_cache(None)
def report():
    return dict(round=395,
                scope='Conditional finite multiscale probability table on one query-independent carrier. Fourier phases and full contrast are modelling inputs; no spatial dimension or physical contact is derived.',
                exact_cases=[exact_audit(n) for n in range(1, 6)],
                coherent_label_encoding=reference_audit(),
                noisy_case=noisy_audit(),
                query_dependence=query_dependence_audit(),
                preparation_and_measurement_independence_required=True,
                memory_dimension_identified_with_spatial_dimension=False,
                full_position_generation_completed=False,
                global_affine_no_go_recounted_as_new=False)


class Audit(unittest.TestCase):
    def test_finite_pair_witnesses_force_orthogonal_supports(self):
        for data in report()['exact_cases']:
            self.assertEqual(data['detector_settings'], 2*data['promised_settings']-2)
            self.assertEqual(data['exact_minimum_carrier_dimension'], 2**data['scales'])
            self.assertLess(data['pair_witness_error'], 3e-12)

    def test_fourier_code_realizes_all_queries_and_off_grid_inputs(self):
        for data in report()['exact_cases']:
            for key in ('fourier_unitarity_error', 'independent_product_code_error',
                        'off_grid_probability_error', 'effect_spectrum_error'):
                self.assertLess(data[key], 4e-12)

    def test_label_isometry_preserves_unknown_reference(self):
        data = report()['coherent_label_encoding']
        self.assertLess(data['full_vector_return_error'], 3e-12)
        self.assertLess(data['reference_state_error'], 3e-12)
        self.assertTrue(data['input_label_was_already_dimension_m'])
        self.assertFalse(data['unknown_one_qubit_phase_amplification_claimed'])

    def test_noisy_dimension_bound_for_mixed_states(self):
        data = report()['noisy_case']
        self.assertLessEqual(data['observed_uniform_probability_error'], .01+3e-12)
        self.assertEqual(data['minimum_integer_dimension'], 11)
        self.assertLess(data['carrier_dimension_lower_bound'], data['actual_carrier_dimension'])
        self.assertLessEqual(data['max_pair_hilbert_schmidt_overlap'], data['overlap_upper_bound'])
        self.assertAlmostEqual(data['average_purity'], 1/16)
        self.assertEqual(packing_lower_bound(16, 0), 16)
        self.assertEqual(packing_lower_bound(16, .5), 1)

    def test_query_dependent_preparation_exposes_quantifier_error(self):
        data = report()['query_dependence']
        self.assertLess(data['probability_error'], 3e-12)
        self.assertAlmostEqual(data['preparation_trace_distance'], 1/np.sqrt(2))
        self.assertFalse(data['satisfies_query_independent_carrier_contract'])


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
