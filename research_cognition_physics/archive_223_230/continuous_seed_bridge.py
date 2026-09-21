"""Round 229: one continuous reversible seed and dense controls give exact control.

The proof is analytic (see research_note_229.md). These NumPy certificates
check the finite conjugation/Jacobian construction, its stability and boundaries.
Sampled conjugators are mathematical witnesses, not a simulated implementation
of an unspecified actual gate set. No frozen result writer is called.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import platform
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from reversible_dynamics_bridge import dagger, hermitian_basis, sample_generator, unitary


def traceless(h):
    return h - np.trace(h) * np.eye(len(h)) / len(h)


def coordinates(matrices, basis):
    """Real Hilbert-Schmidt coordinates, one matrix per column."""
    values = np.einsum('aij,kji->ak', np.asarray(basis), np.asarray(matrices))
    if np.max(np.abs(values.imag)) > 1e-10:
        raise ValueError('Expected Hermitian matrices')
    return values.real


def random_unitary(n, rng):
    q, _ = np.linalg.qr(rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n)))
    return q


def seed(n):
    if n == 3:
        k = np.diag([1., 1., -2.])
    elif n >= 2 and n % 2 == 0:
        # A local two-level clock with a spectator: highly degenerate spectrum.
        k = np.kron(np.diag([1., -1.]), np.eye(n//2))
    else:
        raise ValueError('Certificate dimensions are 2, 3 and even n >= 4')
    return k / np.linalg.norm(k)


@lru_cache(maxsize=None)
def orbit_frame(n):
    """Choose n^2-1 independent conjugates out of a fixed seeded witness pool."""
    rng = np.random.default_rng(22900+n)
    m = n*n-1
    k = seed(n)
    basis = hermitian_basis(n)[1:]
    gates = [random_unitary(n, rng) for _ in range(4*m)]
    conjugates = [g@k@dagger(g) for g in gates]
    design = coordinates(conjugates, basis)
    remaining = design.copy()
    chosen = []
    for _ in range(m):
        j = int(np.argmax(np.sum(remaining*remaining, axis=0)))
        length = np.linalg.norm(remaining[:, j])
        if length < 1e-10:
            raise ArithmeticError('Witness pool failed to span the tangent space')
        direction = remaining[:, j] / length
        chosen.append(j)
        # Two projections keep accumulated roundoff below the rank threshold.
        for _ in range(2):
            remaining -= np.outer(direction, direction@remaining)
        remaining[:, chosen] = 0
    return k, [gates[j] for j in chosen], [conjugates[j] for j in chosen], design[:, chosen]


def pulse_product(generators, times, derivatives=False):
    n = len(generators[0])
    pulses = [unitary(k, t) for k, t in zip(generators, times)]
    prefix = [np.eye(n, dtype=complex)]
    for pulse in pulses:
        prefix.append(prefix[-1]@pulse)
    if not derivatives:
        return prefix[-1]
    suffix = [None]*(len(pulses)+1)
    suffix[-1] = np.eye(n, dtype=complex)
    for j in reversed(range(len(pulses))):
        suffix[j] = pulses[j]@suffix[j+1]
    jacobian = [prefix[j]@(-1j*k@pulses[j])@suffix[j+1]
                for j, k in enumerate(generators)]
    return prefix[-1], jacobian


def real_vector(matrix):
    return np.concatenate((matrix.real.ravel(), matrix.imag.ravel()))


def solve_near_identity(generators, target):
    """Finite numerical check of the local chart, not an exact symbolic solver."""
    times = np.zeros(len(generators))
    residuals = []
    for _ in range(25):
        value, derivatives = pulse_product(generators, times, True)
        residual = real_vector(value-target)
        error = float(np.linalg.norm(residual))
        residuals.append(error)
        if error < 2e-13:
            break
        jacobian = np.column_stack([real_vector(x) for x in derivatives])
        step = np.linalg.lstsq(jacobian, -residual, rcond=None)[0]
        for power in range(20):
            trial = times + step/(2**power)
            if np.linalg.norm(pulse_product(generators, trial)-target) < error:
                times = trial
                break
        else:
            raise ArithmeticError('Local inverse iteration did not descend')
    if residuals[-1] >= 2e-13:
        raise ArithmeticError('Local inverse iteration did not converge')
    return times, residuals


def rank_certificate(n):
    k, gates, generators, design = orbit_frame(n)
    singular = np.linalg.svd(design, compute_uv=False)
    m = n*n-1
    # ||gkg* - vkv*||_HS <= 2 ||k||_HS ||g-v||_op; ||k||_HS = 1.
    epsilon = float(singular[-1]/(8*np.sqrt(m)))
    rng = np.random.default_rng(229100+n)
    perturbations = []
    nearby_generators = []
    for g in gates:
        h = sample_generator(n, int(rng.integers(1, 1000000)))
        h /= np.linalg.norm(h, 2)
        v = unitary(h, epsilon)@g
        perturbations.append(float(np.linalg.norm(v-g, 2)))
        nearby_generators.append(v@k@dagger(v))
    nearby_design = coordinates(nearby_generators, hermitian_basis(n)[1:])
    nearby_singular = np.linalg.svd(nearby_design, compute_uv=False)
    eigenvalues, multiplicities = np.unique(np.diag(k).real, return_counts=True)
    return {
        'dimension': n, 'target_rank': m,
        'observed_rank': int(np.linalg.matrix_rank(design, tol=1e-10)),
        'seed_eigenvalues': eigenvalues.tolist(), 'seed_multiplicities': multiplicities.tolist(),
        'conjugates': len(generators), 'smallest_singular_value': float(singular[-1]),
        'largest_singular_value': float(singular[0]),
        'allowed_conjugator_operator_error': epsilon,
        'actual_sampled_operator_error': max(perturbations),
        'design_perturbation_norm': float(np.linalg.norm(nearby_design-design, 2)),
        'analytic_design_perturbation_bound': float(2*np.sqrt(m)*epsilon),
        'nearby_smallest_singular_value': float(nearby_singular[-1]),
        'analytic_smallest_singular_lower_bound': float(singular[-1]-2*np.sqrt(m)*epsilon),
    }


def inverse_certificate(n):
    _, _, generators, _ = orbit_frame(n)
    target_generator = sample_generator(n, 229500)
    target_generator /= np.linalg.norm(target_generator)
    target = unitary(target_generator, .002)
    times, errors = solve_near_identity(generators, target)
    value = pulse_product(generators, times)
    return {'dimension': n, 'pulse_factors': len(times), 'iterations': len(errors)-1,
            'residuals': errors, 'final_unitary_error': float(np.linalg.norm(value-target)),
            'largest_absolute_pulse_time': float(np.max(np.abs(times))),
            'unitarity_error': float(np.linalg.norm(dagger(value)@value-np.eye(n))),
            'determinant_error': float(abs(np.linalg.det(value)-1))}


def local_only_rank():
    rng = np.random.default_rng(229777)
    k = seed(4)
    gates = [np.kron(random_unitary(2, rng), random_unitary(2, rng)) for _ in range(60)]
    design = coordinates([g@k@dagger(g) for g in gates], hermitian_basis(4)[1:])
    return int(np.linalg.matrix_rank(design, tol=1e-10))


class ContinuousSeedTests(unittest.TestCase):
    def test_degenerate_seed_conjugates_span_all_traceless_directions(self):
        for n in (2, 3, 4, 6, 8):
            with self.subTest(n=n):
                c = rank_certificate(n)
                self.assertEqual(c['observed_rank'], n*n-1)
                self.assertGreater(c['smallest_singular_value'], .01)

    def test_finite_conjugator_errors_preserve_full_rank(self):
        for n in (2, 3, 4, 6, 8):
            with self.subTest(n=n):
                c = rank_certificate(n)
                self.assertLessEqual(c['actual_sampled_operator_error'], c['allowed_conjugator_operator_error']+1e-13)
                self.assertLessEqual(c['design_perturbation_norm'], c['analytic_design_perturbation_bound']+1e-12)
                self.assertGreaterEqual(c['nearby_smallest_singular_value'], c['analytic_smallest_singular_lower_bound']-1e-12)

    def test_conjugated_pulses_are_the_same_actual_one_parameter_group(self):
        k, gates, generators, _ = orbit_frame(4)
        for g, h in zip(gates, generators):
            np.testing.assert_allclose(unitary(h, .137), g@unitary(k, .137)@dagger(g), atol=5e-14)
            np.testing.assert_allclose(unitary(h, .4)@unitary(h, -.7), unitary(h, -.3), atol=5e-14)

    def test_product_jacobian_matches_independent_central_differences(self):
        _, _, generators, _ = orbit_frame(4)
        rng = np.random.default_rng(229888)
        times = rng.normal(size=len(generators))*.01
        _, analytic = pulse_product(generators, times, True)
        for j, derivative in enumerate(analytic):
            step = np.zeros(len(times)); step[j] = 1e-6
            numerical = (pulse_product(generators, times+step)-pulse_product(generators, times-step))/(2e-6)
            np.testing.assert_allclose(derivative, numerical, atol=4e-9)

    def test_local_inverse_reaches_independently_chosen_targets(self):
        for n in (2, 4):
            c = inverse_certificate(n)
            self.assertLess(c['final_unitary_error'], 2e-13)
            self.assertLess(c['unitarity_error'], 5e-13)
            self.assertLess(c['determinant_error'], 5e-13)

    def test_scalar_hamiltonian_has_no_projective_direction(self):
        k = 7*np.eye(4)
        np.testing.assert_allclose(traceless(k), 0, atol=0)
        rho = np.ones((4, 4))/4
        v = unitary(k, .27)
        np.testing.assert_allclose(v@rho@dagger(v), rho, atol=1e-14)

    def test_local_conjugators_cannot_replace_dense_global_controls(self):
        self.assertEqual(local_only_rank(), 3)
        self.assertLess(local_only_rank(), 4*4-1)

    def test_convex_continuity_does_not_supply_reversible_continuity(self):
        plus = np.ones((2, 2))/2
        z = np.diag([1., -1.])
        for p in (.1, .5, .9):
            output = (1-p)*plus + p*z@plus@z
            self.assertAlmostEqual(float(np.trace(output@output)), 1-2*p*(1-p))
            self.assertLess(float(np.trace(output@output)), 1)


def report():
    return {
        'round': 229, 'date': '2026-09-21',
        'hypothesis': 'In the full F+U+C+P theory, one nonconstant continuous reversible operation path on one object suffices for exact operation completeness on every existing object.',
        'analytic_result': 'For a nontrivial theory, existence of one reversible path, existence of one nontrivial continuous reversible one-parameter group, all unitary channels and all finite quantum instruments are equivalent under F+U+C+P.',
        'extra_input': 'Existence of that path (or Time with a nonscalar generator) is NOT inferred from bare C. A trivial Time group is insufficient for this proof.',
        'proof_inputs': ['round227 dense double-copy reversible group', 'Yamabe path-connected subgroup theorem', 'normal path component and simple su(N)', 'finite conjugate basis and inverse function theorem', 'round228 finite instrument dilation'],
        'rank_certificates': [rank_certificate(n) for n in (2, 3, 4, 6, 8)],
        'local_inverse_certificates': [inverse_certificate(n) for n in (2, 4)],
        'boundaries': {'scalar_projective_seed_rank': 0, 'two_qubit_local_only_seed_orbit_rank': local_only_rank(), 'two_qubit_full_rank': 15, 'convex_identity_phase_flip_at_half_pure_input_output_purity': .5},
        'numerical_scope': 'Seeded mathematical matrices validate finite algebra and local inversion. They do not construct an unknown actual dense subgroup or prove its continuity. No finite precision result is called exact.',
        'unresolved': 'Whether bare F+U+C+P already excludes the no-reversible-path branch. No complete countermodel is supplied.',
        'not_derived': ['a selected physical Hamiltonian or clock scale', 'efficient finite control cost', 'infinite-dimensional field theory or gravity'],
        'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ContinuousSeedTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        target = HERE/'continuous_seed_bridge_results.json'
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing result differs; review before replacing it.')
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
