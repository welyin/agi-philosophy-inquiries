"""Round 226: same-dimension full-ensemble steering fixes resource and POVMs.

Uses the coherent complex tensor representation established in round 225.
The permission conclusion concerns the steering auxiliary, not every target.
"""
import argparse
import json
from pathlib import Path
import platform
import unittest

import numpy as np

from reversible_dynamics_bridge import (
    dagger, hermitian_basis, real_superoperator, sample_generator, sample_state, unitary,
)
from tensor_process_bridge import choi, frame, maximally_entangled, partial_trace


def positive_power(x, power):
    vals, vecs = np.linalg.eigh(x)
    if vals.min() <= 0:
        raise ValueError('This certificate requires a positive definite matrix')
    return (vecs*(vals**power))@dagger(vecs)


def steer(w, e, target_dim, auxiliary_dim):
    return np.einsum('iajb,ba->ij', w.reshape(target_dim, auxiliary_dim, target_dim, auxiliary_dim), e)


def resource(v):
    vector = v.ravel()
    return np.outer(vector, vector.conj())


def random_povm(d, outcomes):
    rng = np.random.default_rng(22600+d+outcomes)
    mats = [rng.normal(size=(d, d))+1j*rng.normal(size=(d, d)) for _ in range(outcomes)]
    positives = [dagger(m)@m for m in mats]
    invsqrt = positive_power(sum(positives), -.5)
    return [invsqrt@x@invsqrt for x in positives]


def instance(d):
    rho = (sample_state(d, seed=226)+np.eye(d)/d)/2
    v = positive_power(rho, .5) @ unitary(sample_generator(d, seed=2261), .271)
    return rho, v, resource(v)


def resource_certificate(d):
    rho, v, w = instance(d)
    vinv = np.linalg.inv(v)
    effects = random_povm(d, d+2)
    ensemble = [steer(w, e, d, d) for e in effects]
    recovered = [(vinv@tau@dagger(vinv)).T for tau in ensemble]
    action = lambda e: steer(w, e, d, d)
    matrix = real_superoperator(action, hermitian_basis(d))
    wrong_w = choi(lambda x: v@x.T@dagger(v), d)
    wrong_eigs = np.linalg.eigvalsh(wrong_w)
    return {'dimension': d, 'target_min_eigenvalue': float(np.linalg.eigvalsh(rho).min()),
            'conditioning_rank': int(np.linalg.matrix_rank(matrix)),
            'conditioning_condition_number': float(np.linalg.cond(matrix)),
            'resource_rank': int(np.linalg.matrix_rank(w)),
            'resource_purity': float(np.trace(w@w).real),
            'target_marginal_error': float(np.linalg.norm(partial_trace(w, d, d, 'A')-rho)),
            'ensemble_sum_error': float(np.linalg.norm(sum(ensemble)-rho)),
            'povm_recovery_error': float(max(np.linalg.norm(e-f) for e, f in zip(effects, recovered))),
            'povm_sum_error': float(np.linalg.norm(sum(recovered)-np.eye(d))),
            'choi_reconstruction_error': float(np.linalg.norm(choi(lambda x: action(x.T), d)-w)),
            'wrong_order_branch_negative_eigenvalues': int(np.sum(wrong_eigs < -1e-12)),
            'expected_wrong_branch_negative_eigenvalues': d*(d-1)//2,
            'wrong_branch_min_eigenvalue': float(wrong_eigs.min())}


def noisy_resource_certificate(d, visibility):
    if not 0 < visibility <= 1:
        raise ValueError('Inverse formula is only valid for nonzero visibility')
    w = visibility*maximally_entangled(d)+(1-visibility)*np.eye(d*d)/(d*d)
    p = frame(d, seed=2262)[0]
    tau = p/d
    required_effect = (p.T-(1-visibility)*np.eye(d)/d)/visibility
    return {'dimension': d, 'visibility': visibility,
            'resource_purity': float(np.trace(w@w).real),
            'target_marginal_error': float(np.linalg.norm(partial_trace(w, d, d, 'A')-np.eye(d)/d)),
            'required_effect_min_eigenvalue': float(np.linalg.eigvalsh(required_effect).min()),
            'analytic_required_effect_min_eigenvalue': -(1-visibility)/(visibility*d),
            'inverse_equation_error': float(np.linalg.norm(steer(w, required_effect, d, d)-tau))}


def boundary_certificate():
    rho = np.diag([1., 0.])
    aux = np.diag([.3, .7])
    w = np.kron(rho, aux)
    probabilities = [.2, .3, .5]
    errors = [np.linalg.norm(steer(w, p*np.eye(2), 2, 2)-p*rho) for p in probabilities]
    action = lambda e: steer(w, e, 2, 2)
    return {'case': 'pure target, mixed same-dimension resource',
            'resource_rank': int(np.linalg.matrix_rank(w)),
            'resource_purity': float(np.trace(w@w).real),
            'conditioning_rank': int(np.linalg.matrix_rank(real_superoperator(action, hermitian_basis(2)))),
            'ensemble_error': float(max(errors)),
            'scope': 'All decompositions of a pure target are scalar multiples; full-rank premise is necessary.'}


def larger_auxiliary_certificate():
    # Target B and R0 share a Bell state; the auxiliary also has a mixed flag F.
    w = np.kron(maximally_entangled(2), np.eye(2)/2)
    effects = random_povm(2, 5)
    lifted = [np.kron(e, np.eye(2)) for e in effects]
    errors = [np.linalg.norm(steer(w, e, 2, 4)-f.T/2) for e, f in zip(lifted, effects)]
    image = np.column_stack([
        np.array([np.trace(b@steer(w, e, 2, 4)).real for b in hermitian_basis(2)])
        for e in hermitian_basis(4)])
    return {'target_dim': 2, 'auxiliary_dim': 4,
            'resource_rank': int(np.linalg.matrix_rank(w)),
            'resource_purity': float(np.trace(w@w).real),
            'conditioning_rank': int(np.linalg.matrix_rank(image)),
            'conditioning_nullity': int(image.shape[1]-np.linalg.matrix_rank(image)),
            'lifted_povm_sum_error': float(np.linalg.norm(sum(lifted)-np.eye(4))),
            'ensemble_error': float(max(errors)),
            'scope': 'A mixed resource still steers every ensemble if a redundant mixed flag enlarges the auxiliary.'}


class SteeringPermissionTests(unittest.TestCase):
    def test_full_rank_same_dimension_conditioning_is_invertible(self):
        for d in (2, 3, 4):
            c = resource_certificate(d)
            self.assertEqual(c['conditioning_rank'], d*d)
            self.assertLess(c['target_marginal_error'], 3e-13)

    def test_all_outcomes_of_povm_are_recovered_together(self):
        for d in (2, 3, 4):
            c = resource_certificate(d)
            self.assertLess(c['povm_recovery_error'], 3e-13)
            self.assertLess(c['povm_sum_error'], 3e-13)
            self.assertLess(c['ensemble_sum_error'], 3e-13)

    def test_sharp_projective_effects_are_also_unique_preimages(self):
        rho, v, w = instance(3)
        vinv = np.linalg.inv(v)
        for e in frame(3, seed=2263):
            tau = steer(w, e, 3, 3)
            recovered = (vinv@tau@dagger(vinv)).T
            np.testing.assert_allclose(recovered, e, atol=2e-14)
            np.testing.assert_allclose(recovered@recovered, recovered, atol=2e-14)

    def test_positive_joint_state_selects_the_transpose_conditioning_branch(self):
        for d in (2, 3, 4):
            c = resource_certificate(d)
            self.assertEqual(c['resource_rank'], 1)
            self.assertAlmostEqual(c['resource_purity'], 1)
            self.assertLess(c['choi_reconstruction_error'], 3e-13)
            self.assertEqual(c['wrong_order_branch_negative_eigenvalues'], d*(d-1)//2)

    def test_maximally_mixed_target_requires_maximal_entanglement(self):
        for d in (2, 3, 4):
            v = unitary(sample_generator(d, seed=2264), .339)/np.sqrt(d)
            w = resource(v)
            np.testing.assert_allclose(partial_trace(w, d, d, 'A'), np.eye(d)/d, atol=2e-14)
            np.testing.assert_allclose(partial_trace(w, d, d, 'B'), np.eye(d)/d, atol=2e-14)
            self.assertAlmostEqual(np.trace(w@w).real, 1)

    def test_every_nonzero_isotropic_noise_blocks_an_exact_pure_branch(self):
        for d in (2, 3, 4):
            for visibility in (.5, .9, .99, 1.):
                c = noisy_resource_certificate(d, visibility)
                self.assertAlmostEqual(c['required_effect_min_eigenvalue'], c['analytic_required_effect_min_eigenvalue'])
                self.assertLess(c['inverse_equation_error'], 3e-13)
                if visibility < 1:
                    self.assertLess(c['required_effect_min_eigenvalue'], -1e-4)

    def test_pure_target_is_an_explicit_boundary_exception(self):
        c = boundary_certificate()
        self.assertEqual(c['resource_rank'], 2)
        self.assertEqual(c['conditioning_rank'], 1)
        self.assertAlmostEqual(c['resource_purity'], .58)
        self.assertLess(c['ensemble_error'], 2e-14)

    def test_mixed_flag_is_an_explicit_unequal_dimension_exception(self):
        c = larger_auxiliary_certificate()
        self.assertEqual(c['resource_rank'], 2)
        self.assertAlmostEqual(c['resource_purity'], .5)
        self.assertEqual(c['conditioning_rank'], 4)
        self.assertEqual(c['conditioning_nullity'], 12)
        self.assertLess(c['ensemble_error'], 2e-14)
        self.assertLess(c['lifted_povm_sum_error'], 2e-14)


def report():
    return {'round': 226, 'date': '2026-09-20',
            'analytic_result': 'In the coherent quantum representation, a same-Hilbert-dimension state steering every finite ensemble of a full-rank target is pure with full Schmidt rank, and the steering auxiliary admits every finite POVM. Conversely, such a pure resource and full auxiliary POVMs implement all ensembles.',
            'inherited_input': 'U with its fixed auxiliary type, equal real linear dimension, and full finite-ensemble quantifiers; round225 coherent quantum tensors.',
            'additional_physical_axioms': [],
            'resource_certificates': [resource_certificate(d) for d in (2, 3, 4)],
            'noisy_resource_certificates': [noisy_resource_certificate(d, p) for d in (2, 3, 4) for p in (.5, .9, .99, 1.)],
            'boundary_examples': [boundary_certificate(), larger_auxiliary_certificate()],
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
            'permission_scope': 'All finite POVMs on the designated steering auxiliary. Does not yet establish all measurements on every target type or all quantum instruments.',
            'not_derived': ['deterministic purification of an unknown state', 'uniqueness of purification under actual allowed reversible auxiliary operations', 'availability of every CPTP map', 'a particular Hamiltonian or spacetime'],
            'next': 'Audit whether C and finite composition supply enough reversible operations to transfer permissions and implement all channels; state any needed topological closure explicitly.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SteeringPermissionTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        target = Path(__file__).with_name('steering_permission_bridge_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing result differs; review before replacing it.')
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
