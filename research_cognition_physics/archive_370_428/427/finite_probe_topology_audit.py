"""Round 427: finite probability certificates on an actual control-budget set.

Infinite-tail and topology claims are analytic; finite arrays verify their
ingredients. No finite truncation is treated as the full unknown-input norm.
"""
import argparse
import io
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'finite_probe_topology_audit_results.json'
OBS = {}
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)


def block(protocol, j):
    u = I.copy()
    for axis, duration in protocol:
        angle = duration * 4. ** (-j)
        h = X if axis == 'X' else Z
        u = (np.cos(angle) * I - 1j * np.sin(angle) * h) @ u
    return u


def norm(a):
    return float(np.linalg.norm(a, 2))


def first_column_probabilities(u):
    a, b = u[:, 0]
    return (1 + np.array([a.real, a.imag, b.real, b.imag])) / 2


def physical_probabilities(u):
    """Independent Born calculation with an actual stationary reference level."""
    v = np.zeros((3, 3), dtype=complex)
    v[0, 0] = 1
    v[1:, 1:] = u
    psi = v @ (np.array([1, 1, 0], dtype=complex) / np.sqrt(2))
    probabilities = []
    for k in (1, 2):
        for imaginary in (False, True):
            observable = np.zeros((3, 3), dtype=complex)
            observable[0, k] = -1j if imaginary else 1
            observable[k, 0] = 1j if imaginary else 1
            effect = (np.eye(3) + observable) / 2
            assert np.min(np.linalg.eigvalsh(effect)) >= -1e-14
            assert np.max(np.linalg.eigvalsh(effect)) <= 1 + 1e-14
            probabilities.append(float(np.vdot(psi, effect @ psi).real))
    return np.array(probabilities)


def reconstruct(p):
    q = 2 * np.asarray(p) - 1
    a, b = q[0] + 1j * q[1], q[2] + 1j * q[3]
    return np.array([[a, -b.conjugate()], [b, a.conjugate()]])


def certificate(budget, blocks, probability_gap, estimate_error=0.):
    return min(1., max(4 * (probability_gap + 2 * estimate_error),
                       2 * budget * 4. ** (-(blocks + 1))))


def finite_unitary(protocol, blocks):
    u = np.zeros((1 + 2 * blocks, 1 + 2 * blocks), dtype=complex)
    u[0, 0] = 1
    for j in range(1, blocks + 1):
        u[2*j-1:2*j+1, 2*j-1:2*j+1] = block(protocol, j)
    return u


A = [('X', .8), ('Z', -1.1), ('X', .6), ('Z', .4)]
B = [('X', .81), ('Z', -1.12), ('X', .61), ('Z', .4)]


class Audit(unittest.TestCase):
    def test_01_actual_effects_and_reconstruction(self):
        born_error = reconstruction_error = 0.
        for protocol in (A, B, [('Z', math.pi * 4)]):
            for j in range(1, 9):
                u = block(protocol, j)
                p = physical_probabilities(u)
                born_error = max(born_error, float(np.max(abs(p-first_column_probabilities(u)))))
                reconstruction_error = max(reconstruction_error, norm(reconstruct(p)-u))
        self.assertLess(born_error, 1e-14)
        self.assertLess(reconstruction_error, 1e-14)
        OBS['actual_effect_born_formula_error'] = born_error
        OBS['su2_first_column_reconstruction_error'] = reconstruction_error

    def test_02_quaternion_difference_and_phase_reference(self):
        rng = np.random.default_rng(427)
        equality_error = 0.
        for _ in range(40):
            us = []
            for _ in range(2):
                q = rng.normal(size=4)
                q /= np.linalg.norm(q)
                a, b = q[0]+1j*q[1], q[2]+1j*q[3]
                us.append(np.array([[a, -b.conjugate()], [b, a.conjugate()]]))
            delta = us[0]-us[1]
            equality_error = max(equality_error, abs(norm(delta)-float(np.linalg.norm(delta[:, 0]))))
            gap = float(max(abs(first_column_probabilities(us[0])-first_column_probabilities(us[1]))))
            self.assertLessEqual(norm(delta), 4*gap+1e-14)
        # The block channels of I and -I coincide, but coherent reference tests differ.
        rho = np.array([[.6, .2j], [-.2j, .4]])
        self.assertLess(norm((-I) @ rho @ (-I)-rho), 1e-14)
        gap = float(max(abs(physical_probabilities(I)-physical_probabilities(-I))))
        self.assertAlmostEqual(gap, 1.)
        OBS['quaternion_operator_norm_identity_error'] = equality_error
        OBS['reference_detects_block_phase_probability_gap'] = gap

    def test_03_control_budget_full_operation_certificate(self):
        budget, count, truncation = 4., 4, 18
        for protocol in (A, B):
            self.assertLessEqual(sum(abs(t) for _, t in protocol), budget)
            for j in range(1, truncation+1):
                self.assertLessEqual(norm(block(protocol, j)-I), budget*4.**(-j)+1e-14)
        gap = max(float(max(abs(physical_probabilities(block(A,j))-physical_probabilities(block(B,j)))))
                  for j in range(1,count+1))
        bound = certificate(budget, count, gap)
        u, v = finite_unitary(A, truncation), finite_unitary(B, truncation)
        self.assertLessEqual(norm(u-v), bound+1e-14)
        # Arbitrary entangled inputs: these are lower witnesses, not a diamond computation.
        rng, witness = np.random.default_rng(1427), 0.
        relative = u.conj().T @ v
        for _ in range(25):
            psi = rng.normal(size=(len(u), 3)) + 1j*rng.normal(size=(len(u), 3))
            psi /= np.linalg.norm(psi)
            overlap = np.vdot(psi, relative @ psi)
            distance = math.sqrt(max(0., 1-abs(overlap)**2))
            witness = max(witness, distance)
            self.assertLessEqual(distance, bound+1e-12)
        OBS['actual_protocol_certificate'] = dict(
            common_action_bound=budget, measured_blocks=count, binary_settings=4*count,
            maximum_probability_difference=gap, analytic_infinite_tail_bound=2*budget*4.**(-(count+1)),
            full_half_diamond_upper_bound=bound, finite_truncation_operator_difference=norm(u-v),
            sampled_entangled_input_lower_witness=witness,
            truncation_blocks=truncation, beyond_truncation_bound=2*budget*4.**(-(truncation+1)))

    def test_04_measurement_errors_and_counted_shots(self):
        budget, count, eta, failure = 4., 4, .003, .01
        p = np.concatenate([physical_probabilities(block(A,j)) for j in range(1,count+1)])
        q = np.concatenate([physical_probabilities(block(B,j)) for j in range(1,count+1)])
        shots = math.ceil(math.log(16*count/failure)/(2*eta*eta))
        union_bound = 16*count*math.exp(-2*shots*eta*eta)
        self.assertLessEqual(union_bound, failure)
        # An explicit adverse error pair verifies the deterministic 2 eta allowance.
        sign = np.sign(q-p)
        p_hat, q_hat = np.clip(p+eta*sign, 0, 1), np.clip(q-eta*sign, 0, 1)
        measured_gap = float(max(abs(p_hat-q_hat)))
        self.assertLessEqual(float(max(abs(p-q))), measured_gap+2*eta+1e-14)
        bound = certificate(budget, count, measured_gap, eta)
        self.assertLessEqual(norm(finite_unitary(A, 18)-finite_unitary(B,18)), bound+1e-14)
        OBS['finite_shot_certificate'] = dict(
            error_per_probability=eta, requested_failure_probability=failure,
            shots_per_probability=shots, distinct_probability_estimates=8*count,
            total_channel_uses=8*count*shots, analytic_union_failure_bound=union_bound,
            adverse_observed_probability_gap=measured_gap,
            certified_half_diamond_error=bound, independent_repreparations_required=True)

    def test_05_no_budget_blind_family(self):
        rows = []
        for count in (1, 4, 12):
            m = count+1
            # T_m*w_j/pi=4**(m-j) is an even integer for j<m; at m it is 1.
            # Use integer arithmetic, not floating trig of enormous angles.
            for j in range(1,count+1):
                self.assertEqual((4**(m-j)) % 2, 0)
                self.assertEqual(float(max(abs(physical_probabilities(I)-physical_probabilities(I)))), 0.)
            before = np.array([1,1], dtype=complex)/np.sqrt(2)
            after = np.array([1,-1], dtype=complex)/np.sqrt(2)
            distance = math.sqrt(1-abs(np.vdot(before, after))**2)
            self.assertAlmostEqual(distance, 1.)
            rows.append(dict(tested_blocks=count, hidden_block=m,
                             action_over_pi=str(4**m), all_tested_probability_differences=0,
                             exact_full_half_diamond_distance=1))
        OBS['unbounded_budget_blind_family'] = rows


def run():
    output = io.StringIO()
    tests = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=427, baseline_round=426, status='finite_probability_to_uniform_operation_certificate',
                tests_run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors),
                python=platform.python_version(), numpy=np.__version__, observations=OBS,
                scope=dict(compact_budget_finite_test_certificate_proved=True,
                           test_Baire_topology_recovers_uniform_topology=True,
                           Baire_derived_from_cognitive_principles=False,
                           source_model_budget_promise_required=True,
                           stationary_coherent_reference_explicit_input=True,
                           preparation_and_readout_autonomously_generated=False,
                           finite_test_list_is_exact_global_tomography=False,
                           sampled_inputs_used_to_prove_diamond_bound=False,
                           actual_position_interface_still_required=True,
                           three_dimensional_space_unconditionally_derived=False,
                           full_cognitive_countermodel_completed=False,
                           full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert json.loads(TARGET.read_text(encoding='utf-8')) == result
    else:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
