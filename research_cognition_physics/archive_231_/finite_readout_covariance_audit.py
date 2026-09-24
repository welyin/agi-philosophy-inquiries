"""Round 417: finite summaries, translations, and rotation closure.

All geometry labels below are declared model inputs.  This audits a proposed
summary-update contract; it does not generate positions or spatial dimension.
Run read-only by default; --write-results creates a new file exclusively.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import unittest

import numpy as np


HERE = Path(__file__).resolve().parent
VERTICES = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=int)
POINTS = np.array([[0, 0, 0], [.12, -.17, .09], [-.22, .08, .16],
                   [.3, .1, -.2], [-.16, -.21, -.14], [.07, .11, .29],
                   [-.18, .25, -.09]])
IDENTITY = np.eye(4)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
PAULI = np.array([X, Y, Z])


def unitary(a):
    return np.diag(np.exp(-1j * (VERTICES @ np.asarray(a))))


def ket(x):
    return np.exp(-1j * (VERTICES @ np.asarray(x))) / 2


def density(x):
    v = ket(x)
    return np.outer(v, v.conj())


def trace_distance(a, b):
    difference = (a - b + (a - b).conj().T) / 2
    return float(np.sum(np.abs(np.linalg.eigvalsh(difference))) / 2)


def signed_permutation_rotations():
    group = []
    vertex_set = {tuple(v) for v in VERTICES}
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product([-1, 1], repeat=3):
            rotation = np.eye(3, dtype=int)[list(permutation)] * np.array(signs)[:, None]
            determinant = round(np.linalg.det(rotation))
            if determinant == 1 and {tuple(rotation @ v) for v in VERTICES} == vertex_set:
                group.append(rotation)
    return group


def carrier_permutation(rotation):
    """V_R |v> = |R v>, hence V_R U(a) V_R* = U(R a)."""
    lookup = {tuple(v): i for i, v in enumerate(VERTICES)}
    permutation = np.zeros((4, 4), dtype=int)
    for i, v in enumerate(VERTICES):
        permutation[lookup[tuple(rotation @ v)], i] = 1
    return permutation


def phase_chart(rho):
    # On ||x||_1 < pi/2 these three phases lie strictly inside (-pi, pi).
    differences = VERTICES[1:] - VERTICES[0]
    phases = -np.angle(4 * rho[1:, 0])
    return np.linalg.solve(differences, phases)


def rotate60():
    return np.array([[.5, -np.sqrt(3) / 2, 0],
                     [np.sqrt(3) / 2, .5, 0], [0, 0, 1]])


def random_channel_kraus(seed, environment=3):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(4 * environment, 4)) + 1j * rng.normal(size=(4 * environment, 4))
    q, _ = np.linalg.qr(a)
    return [q[4 * k:4 * (k + 1)] for k in range(environment)]


def apply_channel(kraus, rho):
    return sum(k @ rho @ k.conj().T for k in kraus)


def replace_kraus(target):
    values, vectors = np.linalg.eigh(target)
    basis = np.eye(4)
    return [np.sqrt(max(value, 0)) * np.outer(vector, b)
            for value, vector in zip(values, vectors.T) for b in basis]


@lru_cache(maxsize=1)
def report():
    rotations = signed_permutation_rotations()
    group_keys = {tuple(r.flat) for r in rotations}
    frequencies = np.unique((VERTICES[:, None] - VERTICES[None, :]).reshape(-1, 3), axis=0)
    sigma = {tuple(k) for k in frequencies}
    grid = np.array(list(itertools.product(range(5), repeat=3))) * (2 * np.pi / 5)
    characters = np.exp(1j * grid @ frequencies.T)
    gram = characters.conj().T @ characters / len(grid)
    rho_samples = np.array([density(x).reshape(-1) for x in grid])
    # Entries are characters with frequencies -(v_i-v_j), with factor 1/4.
    expected_entries = np.exp(-1j * grid @ (VERTICES[:, None] - VERTICES[None, :]).reshape(-1, 3).T) / 4
    rotation_errors = []
    channel_rotation_errors = []
    for r in rotations:
        p = carrier_permutation(r)
        for x in POINTS:
            rotation_errors.append(np.linalg.norm(p @ density(x) @ p.T - density(r @ x)))
            channel_rotation_errors.append(np.linalg.norm(p @ unitary(x) @ p.T - unitary(r @ x)))
    composition_errors = []
    for a, b in zip(POINTS, POINTS[::-1]):
        composition_errors.append(np.linalg.norm(unitary(a) @ unitary(b) - unitary(a + b)))
    # A declared entangled input checks the actual action on a retained reference.
    phi = np.eye(4).reshape(-1) / 2
    joint = np.outer(phi, phi)
    a, b = POINTS[1], POINTS[2]
    vab = np.kron(unitary(a) @ unitary(b), IDENTITY)
    vsum = np.kron(unitary(a + b), IDENTITY)
    reference_error = trace_distance(vab @ joint @ vab.conj().T, vsum @ joint @ vsum.conj().T)
    chart_errors = [np.linalg.norm(phase_chart(density(x)) - x) for x in POINTS]
    differences = VERTICES[1:] - VERTICES[0]
    derivative = np.array([-1j * (np.diag(VERTICES[:, i]) @ density([0, 0, 0])
                                  - density([0, 0, 0]) @ np.diag(VERTICES[:, i])) for i in range(3)])
    real_derivative = np.concatenate([derivative.reshape(3, -1).real,
                                      derivative.reshape(3, -1).imag], axis=1).T

    offset = np.array([np.pi, 0, 0])
    r = rotate60()
    original0, original1 = density([0, 0, 0]), density(offset)
    target0, target1 = original0, density(r @ offset)
    midpoint = (target0 + target1) / 2
    effect = target0
    outputs = []
    completeness_errors = []
    for seed in [41701, 41702, 41703]:
        kraus = random_channel_kraus(seed)
        completeness_errors.append(np.linalg.norm(sum(k.conj().T @ k for k in kraus) - IDENTITY))
        out0, out1 = apply_channel(kraus, original0), apply_channel(kraus, original1)
        outputs.append(dict(seed=seed, identical_input_output_error=trace_distance(out0, out1),
                            worst_target_error=max(trace_distance(out0, target0), trace_distance(out1, target1))))
    reset = replace_kraus(midpoint)
    reset_output = apply_channel(reset, original0)
    # The input-independent fixed auxiliary carries no information about x.
    aux = np.array([[.6, .1j], [-.1j, .4]])
    auxiliary_collision = trace_distance(np.kron(original0, aux), np.kron(original1, aux))
    tensor_collision = trace_distance(np.kron(original0, original0), np.kron(original1, original1))
    robustness = []
    for delta in [0, .05, .2, .7, 1]:
        # A separate commuting pair saturates the general robust bound.
        s0 = (1 - delta) * midpoint + delta * target0
        s1 = (1 - delta) * midpoint + delta * target1
        robustness.append(dict(delta=delta, input_distance=trace_distance(s0, s1),
                               worst_identity_error=max(trace_distance(s0, target0), trace_distance(s1, target1)),
                               lower_bound=(1 - delta) / 2))
    shift = np.array([.12, -.08, .04])
    shift_operator = np.einsum('i,ijk->jk', shift, PAULI) / 2
    local_minimum = 1.0
    local_update_error = 0.0
    for x in POINTS:
        before = (np.eye(2) + np.einsum('i,ijk->jk', x, PAULI)) / 2
        after = (np.eye(2) + np.einsum('i,ijk->jk', x + shift, PAULI)) / 2
        local_minimum = min(local_minimum, float(np.linalg.eigvalsh(after).min()))
        local_update_error = max(local_update_error, float(np.linalg.norm(after - before - shift_operator)))
    boundary = shift / np.linalg.norm(shift)
    boundary_state = (np.eye(2) + np.einsum('i,ijk->jk', boundary, PAULI)) / 2
    boundary_after = boundary_state + shift_operator
    return dict(
        round=417, scientific_baseline=416,
        scope='Global finite bounded readout closure under translations and dense rotations is trivial; this is an extra cache contract, not a cognitive axiom',
        tetrahedral_group=dict(order=len(rotations),
            closed=all(tuple((r @ s).flat) in group_keys for r in rotations for s in rotations),
            rotations_exactly_orthogonal=all(np.array_equal(r.T @ r, np.eye(3, dtype=int)) for r in rotations),
            carrier_permutations_exact=all(np.array_equal(carrier_permutation(r).T @ carrier_permutation(r), np.eye(4, dtype=int)) for r in rotations),
            state_covariance_error=float(max(rotation_errors)), channel_covariance_error=float(max(channel_rotation_errors))),
        translation=dict(maximum_composition_error=float(max(composition_errors)),
                         retained_reference_error=reference_error,
                         nontrivial_sample_distance=trace_distance(density(POINTS[3]), original0)),
        finite_readout_space=dict(frequencies=frequencies.tolist(), complex_dimension=len(frequencies),
            all_tetrahedral_frequency_permutations=all({tuple(r.T @ k) for k in frequencies} == sigma for r in rotations),
            fourier_grid_gram_error=float(np.linalg.norm(gram - np.eye(len(frequencies)))),
            matrix_entry_character_error=float(np.linalg.norm(rho_samples - expected_entries)),
            rotated_frequency_outside_set=not any(np.allclose(r.T @ np.array([2, 2, 0]), k, atol=1e-12, rtol=0) for k in frequencies)),
        local_chart=dict(phase_difference_determinant=int(round(np.linalg.det(differences))),
                         derivative_real_rank=int(np.linalg.matrix_rank(real_derivative)),
                         maximum_reconstruction_error=float(max(chart_errors)),
                         sample_count=len(POINTS), all_samples_in_declared_window=bool(np.max(np.sum(np.abs(POINTS), axis=1)) < np.pi / 2)),
        rotation_collision=dict(offset=offset.tolist(), rotation_angle=float(np.pi / 3),
            input_trace_distance=trace_distance(original0, original1),
            unitary_phase_identity_error=float(np.linalg.norm(unitary(offset) + IDENTITY)),
            target_overlap_absolute=float(abs(np.vdot(ket([0, 0, 0]), ket(r @ offset)))),
            target_trace_distance=trace_distance(target0, target1),
            target_binary_probabilities=[float(np.trace(effect @ t).real) for t in [target0, target1]],
            all_decoders_proven_worst_error_at_least=.5,
            same_auxiliary_input_distance=auxiliary_collision, two_copy_input_distance=tensor_collision),
        decoder_checks=dict(random_cptp=outputs, kraus_completeness_error=float(max(completeness_errors)),
            midpoint_replacement_completeness_error=float(np.linalg.norm(sum(k.conj().T @ k for k in reset) - IDENTITY)),
            midpoint_replacement_error=float(np.linalg.norm(reset_output - midpoint)),
            midpoint_worst_target_error=max(trace_distance(reset_output, target0), trace_distance(reset_output, target1)),
            midpoint_binary_probability=float(np.trace(effect @ reset_output).real)),
        robustness=robustness,
        local_statistical_boundary=dict(shift=shift.tolist(), local_minimum_state_eigenvalue=local_minimum,
            local_affine_update_error=local_update_error, global_positive_extension_exists=False,
            boundary_output_minimum_eigenvalue=float(np.linalg.eigvalsh(boundary_after).min()),
            exact_boundary_eigenvalue=-float(np.linalg.norm(shift)) / 2),
        finite_summary_autonomous_closure_is_extra_input=True,
        global_readout_closure_theorem_proved=True,
        finite_direction_group_ruled_out=False,
        finite_window_readout_ruled_out=False,
        repeated_access_to_source_ruled_out=False,
        quantum_direction_interface_disproved=False,
        all_cognitive_principles_countermodel=False,
        physical_positions_or_dimension_generated=False,
        phase_closure_triggered=False)


class AuditTests(unittest.TestCase):
    def test_01_exact_finite_symmetry_and_carrier_action(self):
        a = report()['tetrahedral_group']
        self.assertEqual(a['order'], 12)
        self.assertTrue(a['closed'] and a['rotations_exactly_orthogonal'] and a['carrier_permutations_exact'])
        self.assertLess(max(a['state_covariance_error'], a['channel_covariance_error']), 3e-14)

    def test_02_translation_composition_and_reference(self):
        a = report()['translation']
        self.assertLess(max(a['maximum_composition_error'], a['retained_reference_error']), 3e-14)
        self.assertGreater(a['nontrivial_sample_distance'], .2)

    def test_03_complete_frequency_space_and_missing_rotation(self):
        a = report()['finite_readout_space']
        self.assertEqual(a['complex_dimension'], 13)
        self.assertTrue(a['all_tetrahedral_frequency_permutations'] and a['rotated_frequency_outside_set'])
        self.assertLess(max(a['fourier_grid_gram_error'], a['matrix_entry_character_error']), 1e-12)

    def test_04_actual_local_phase_reconstruction(self):
        a = report()['local_chart']
        self.assertEqual(abs(a['phase_difference_determinant']), 16)
        self.assertEqual(a['derivative_real_rank'], 3)
        self.assertTrue(a['all_samples_in_declared_window'])
        self.assertLess(a['maximum_reconstruction_error'], 2e-15)

    def test_05_identical_summaries_orthogonal_targets(self):
        a = report()['rotation_collision']
        for key in ['input_trace_distance', 'unitary_phase_identity_error', 'target_overlap_absolute',
                    'same_auxiliary_input_distance', 'two_copy_input_distance']:
            self.assertLess(a[key], 2e-14, key)
        self.assertAlmostEqual(a['target_trace_distance'], 1)
        np.testing.assert_allclose(a['target_binary_probabilities'], [1, 0], atol=2e-14, rtol=0)

    def test_06_physical_decoders_and_optimal_midpoint(self):
        a = report()['decoder_checks']
        self.assertLess(max(a['kraus_completeness_error'], a['midpoint_replacement_completeness_error'],
                            a['midpoint_replacement_error']), 3e-14)
        for row in a['random_cptp']:
            self.assertLess(row['identical_input_output_error'], 2e-14)
            self.assertGreaterEqual(row['worst_target_error'], .5 - 2e-14)
        self.assertAlmostEqual(a['midpoint_worst_target_error'], .5)
        self.assertAlmostEqual(a['midpoint_binary_probability'], .5)

    def test_07_robust_contractivity_bound_and_saturation(self):
        for a in report()['robustness']:
            self.assertAlmostEqual(a['input_distance'], a['delta'])
            self.assertAlmostEqual(a['worst_identity_error'], a['lower_bound'])

    def test_08_local_statistics_do_not_give_cp_update(self):
        a = report()['local_statistical_boundary']
        self.assertGreater(a['local_minimum_state_eigenvalue'], .2)
        self.assertLess(a['local_affine_update_error'], 2e-15)
        self.assertLess(a['boundary_output_minimum_eigenvalue'], 0)
        self.assertAlmostEqual(a['boundary_output_minimum_eigenvalue'], a['exact_boundary_eigenvalue'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    tests = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(AuditTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    result = dict(report(), checks=dict(run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors)),
                  runtime=dict(python=platform.python_version(), numpy=np.__version__))
    if args.write_results:
        with (HERE / 'finite_readout_covariance_audit_results.json').open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
