"""Round 217: embedding, joining and sharing are different integration contracts.

The general Bell-sharing bound is proved analytically in the note and agrees
with known isotropic-state extendibility. Finite matrices verify the proof's
identities and attainability; no numerical search stands in for the theorem.
"""
import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import unittest

import numpy as np


def validate(d, k):
    if type(d) is not int or type(k) is not int or d < 2 or k < 1:
        raise ValueError("Integer local dimension d>=2 and partner count k>=1 required.")


def bell_vector(d):
    return np.eye(d).reshape(-1)/np.sqrt(d)


def bell_isometry(d, k, partner):
    """Insert Phi on A,B_partner; other B registers are spectator inputs."""
    validate(d, k)
    if type(partner) is not int or not 1 <= partner <= k:
        raise ValueError("Partner index is between 1 and k.")
    rest = [j for j in range(1, k+1) if j != partner]
    result = np.zeros((d**(k+1), d**(k-1)))
    for column, digits in enumerate(product(range(d), repeat=k-1)):
        for a in range(d):
            labels = [0]*(k+1)
            labels[0] = labels[partner] = a
            for j, value in zip(rest, digits):
                labels[j] = value
            result[np.ravel_multi_index(tuple(labels), (d,)*(k+1)), column] = 1/np.sqrt(d)
    return result


def bell_sum(d, k):
    operators = [bell_isometry(d, k, j) for j in range(1, k+1)]
    return sum(v@v.T for v in operators)


def keep(density, d, sites, total):
    sites = tuple(sites)
    if len(set(sites)) != len(sites) or any(i < 0 or i >= total for i in sites):
        raise ValueError("Distinct in-range retained sites required.")
    if density.shape != (d**total, d**total):
        raise ValueError("Matrix shape does not match subsystem count.")
    rest = tuple(i for i in range(total) if i not in sites)
    order = sites+rest+tuple(total+i for i in sites+rest)
    grouped = density.reshape((d,)*(2*total)).transpose(order)
    return np.einsum('abcb->ac', grouped.reshape(d**len(sites), d**len(rest),
                                               d**len(sites), d**len(rest)))


def sharp_values(d, k):
    validate(d, k)
    return {'maximum_sum_overlap': F(d+k-1, d),
            'maximum_minimum_overlap': F(d+k-1, d*k),
            'minimum_maximum_trace_error': F((d-1)*(k-1), d*k),
            'isotropic_visibility_at_boundary': F(d+k, k*(d+1))}


def isotropic_state(d, overlap):
    overlap = float(overlap)
    if not 0 <= overlap <= 1:
        raise ValueError("Overlap must be in [0,1].")
    vector = bell_vector(d)
    phi = np.outer(vector, vector)
    return overlap*phi+(1-overlap)*(np.eye(d*d)-phi)/(d*d-1)


def top_state(d, k):
    """Numerical representative of the analytically specified top projector."""
    h = bell_sum(d, k)
    values, vectors = np.linalg.eigh(h)
    optimum = float(sharp_values(d, k)['maximum_sum_overlap'])
    selected = vectors[:, np.abs(values-optimum) < 1e-10]
    if not selected.shape[1]:
        raise ArithmeticError("Analytic top eigenspace absent.")
    return selected@selected.T/selected.shape[1]


def trace_distance(first, second):
    return float(np.abs(np.linalg.eigvalsh(first-second)).sum()/2)


def enlarged_center_state(d, k):
    """Order A_1,...,A_k,B_1,...,B_k: retain distinct logical identities."""
    validate(d, k)
    vector = np.eye(d**k).reshape(-1)/np.sqrt(d**k)
    return np.outer(vector, vector)


def classical_star_join(tables):
    """Exact gluing of any finite p(A,B_j) with the same marginal p(A)."""
    tables = [tuple(tuple(F(x) for x in row) for row in table) for table in tables]
    if not tables or not tables[0] or not tables[0][0]:
        raise ValueError("Nonempty probability tables required.")
    dimension_a = len(tables[0])
    marginal = tuple(sum(row) for row in tables[0])
    if sum(marginal) != 1:
        raise ValueError("Each table must be normalized.")
    for table in tables:
        if len(table) != dimension_a or not table[0] or any(len(row) != len(table[0]) for row in table):
            raise ValueError("Each table must be rectangular with the same A.")
        if any(x < 0 for row in table for x in row) or tuple(sum(row) for row in table) != marginal:
            raise ValueError("Nonnegative tables with a common A marginal required.")
    result = {}
    for a, mass in enumerate(marginal):
        for labels in product(*(range(len(t[0])) for t in tables)):
            result[(a,)+labels] = (mass*np.prod([t[a][b]/mass for t, b in zip(tables, labels)])
                                    if mass else F(0))
    return result


class IntegrationExtensionContractTests(unittest.TestCase):
    def test_isometries_and_pair_overlap_identity(self):
        for d, k in ((2, 2), (2, 4), (3, 3)):
            v = [bell_isometry(d, k, j) for j in range(1, k+1)]
            for a in v:
                np.testing.assert_allclose(a.T@a, np.eye(d**(k-1)), atol=8e-16)
            for i in range(k):
                for j in range(i):
                    block = v[i].T@v[j]
                    np.testing.assert_allclose(block@block.T, np.eye(d**(k-1))/d**2, atol=8e-16)

    def test_exact_top_eigenvalue_and_constructive_eigenvector(self):
        for d, k in ((2, 1), (2, 2), (2, 5), (3, 2), (3, 3), (4, 2)):
            h = bell_sum(d, k)
            optimum = float(sharp_values(d, k)['maximum_sum_overlap'])
            vector = sum(bell_isometry(d, k, j)[:, 0] for j in range(1, k+1))
            vector /= np.linalg.norm(vector)
            np.testing.assert_allclose(h@vector, optimum*vector, atol=2e-15)
            self.assertAlmostEqual(np.linalg.eigvalsh(h)[-1], optimum, places=12)

    def test_real_state_attains_sharp_trace_errors_on_all_pairs(self):
        for d, k in ((2, 1), (2, 2), (2, 3), (2, 5), (3, 2), (3, 3)):
            state = top_state(d, k)
            values = sharp_values(d, k)
            target = isotropic_state(d, 1)
            self.assertGreaterEqual(np.linalg.eigvalsh(state).min(), -3e-15)
            self.assertAlmostEqual(np.trace(state), 1., places=13)
            for j in range(1, k+1):
                reduced = keep(state, d, (0, j), k+1)
                np.testing.assert_allclose(reduced, isotropic_state(d, values['maximum_minimum_overlap']), atol=3e-15)
                self.assertAlmostEqual(trace_distance(reduced, target), float(values['minimum_maximum_trace_error']), places=12)

    def test_two_qubit_exact_polynomial_certificate(self):
        # 2H is integer: certify its eigenvalues are among 0,1,3 without eigensolving.
        h = bell_sum(2, 2)
        integer = np.rint(2*h).astype(np.int64)
        np.testing.assert_allclose(2*h, integer, atol=5e-16)
        eye = np.eye(8, dtype=np.int64)
        np.testing.assert_array_equal(integer@(integer-eye)@(integer-3*eye), np.zeros((8, 8), dtype=np.int64))
        numerator = integer@integer-integer
        np.testing.assert_array_equal(numerator@numerator, 6*numerator)
        self.assertEqual(np.trace(numerator), 12)
        state = numerator/12
        np.testing.assert_allclose(state, top_state(2, 2), atol=8e-16)
        self.assertEqual(sharp_values(2, 2)['minimum_maximum_trace_error'], F(1, 4))

    def test_overlap_agreement_does_not_make_bell_demands_joinable(self):
        phi = isotropic_state(2, 1)
        np.testing.assert_allclose(keep(phi, 2, (0,), 2), np.eye(2)/2, atol=3e-16)
        self.assertLess(sharp_values(2, 2)['maximum_sum_overlap'], 2)
        product_extension = np.kron(phi, np.eye(2)/2)
        np.testing.assert_allclose(keep(product_extension, 2, (0, 1), 3), phi, atol=3e-16)
        self.assertAlmostEqual(trace_distance(keep(product_extension, 2, (0, 2), 3), phi), .75)

    def test_more_environment_does_not_change_pair_witness_bound(self):
        h = bell_sum(2, 2)
        enlarged = np.kron(h, np.eye(3))
        self.assertAlmostEqual(np.linalg.eigvalsh(enlarged)[-1], 1.5, places=13)

    def test_full_complex_inputs_also_obey_real_witness(self):
        rng = np.random.default_rng(217)
        h, target = bell_sum(2, 3), isotropic_state(2, 1)
        for _ in range(5):
            raw = rng.normal(size=(16, 4))+1j*rng.normal(size=(16, 4))
            state = raw@raw.conj().T
            state /= np.trace(state)
            self.assertLessEqual(np.trace(h@state).real, 2+1e-14)
            errors = [trace_distance(keep(state, 2, (0, j), 4), target) for j in (1, 2, 3)]
            self.assertGreaterEqual(max(errors), 1/3-1e-14)
            for j in (1, 2, 3):
                self.assertLessEqual(trace_distance(keep(state.real, 2, (0, j), 4), target), errors[j-1]+1e-14)

    def test_classical_star_gluing_with_unequal_alphabets_and_zero_mass(self):
        tables = [((F(1, 6), F(1, 6)), (F(1, 2), F(1, 6)), (0, 0)),
                  ((F(1, 3), 0, 0), (F(1, 6), F(1, 6), F(1, 3)), (0, 0, 0)),
                  ((0, F(1, 3)), (F(1, 3), F(1, 3)), (0, 0))]
        joint = classical_star_join(tables)
        self.assertEqual(sum(joint.values()), 1)
        for j, table in enumerate(tables, 1):
            for a, row in enumerate(table):
                for b, probability in enumerate(row):
                    self.assertEqual(sum(p for x, p in joint.items() if x[0] == a and x[j] == b), probability)

    def test_sharp_error_grows_with_partner_count_but_not_total_impossibility(self):
        for d in (2, 3, 7):
            errors = [sharp_values(d, k)['minimum_maximum_trace_error'] for k in range(1, 21)]
            self.assertEqual(errors[0], 0)
            self.assertTrue(all(a < b for a, b in zip(errors, errors[1:])))
            self.assertLess(errors[-1], F(d-1, d))

    def test_product_composition_preserves_unknown_external_correlations(self):
        rng = np.random.default_rng(218)
        raw = rng.normal(size=(6, 4))
        joint = raw@raw.T
        joint /= np.trace(joint)
        partner = np.diag([.2, .8])
        k = np.array([[.7, .1], [.1, .5]])
        old_branch = np.kron(k, np.eye(3))@joint@np.kron(k.T, np.eye(3))
        extended = np.kron(np.kron(k, np.eye(3)), np.eye(2))
        np.testing.assert_allclose(extended@np.kron(joint, partner)@extended.T,
                                   np.kron(old_branch, partner), atol=3e-16)

    def test_enlarging_center_preserves_distinct_old_bell_relations(self):
        for d, k in ((2, 2), (2, 3), (3, 2)):
            state = enlarged_center_state(d, k)
            phi = isotropic_state(d, 1)
            for j in range(k):
                np.testing.assert_allclose(keep(state, d, (j, k+j), 2*k), phi, atol=5e-16)
            center = keep(state, d, tuple(range(k)), 2*k)
            np.testing.assert_allclose(center, np.eye(d**k)/d**k, atol=5e-16)
            self.assertAlmostEqual(np.trace(state@state), 1., places=13)

    def test_invalid_parameters_or_inconsistent_classical_overlap_rejected(self):
        for d, k in ((1, 2), (2, 0), (True, 2), (2, 1.5)):
            with self.assertRaises(ValueError):
                sharp_values(d, k)
        with self.assertRaises(ValueError):
            classical_star_join([((1, 0), (0, 0)), ((0, 0), (0, 1))])


def report():
    rows = []
    for d, k in ((2, 1), (2, 2), (2, 3), (2, 5), (2, 10), (3, 2), (3, 3)):
        rows.append({'d': d, 'partners': k, **{name: str(value) for name, value in sharp_values(d, k).items()}})
    return {'round': 217, 'date': '2026-09-20',
            'question': 'Which integration requirement constrains state spaces without demanding incompatible old relations?',
            'cognitive_motivation': 'Preserve declared old capabilities while integrating into a larger cognitive whole',
            'additional_model_inputs': ['finite matrix state spaces', 'fixed identities of the shared A and partner subsystems',
                                        'same-time pair marginals', 'full Bell-pair effect and trace distance'],
            'old_rounds_reused': [80, 81, 84, 88, 126, 180, 189, 207, 216],
            'independent_composition_excludes_classical_or_real_models': False,
            'matching_single_party_overlap_suffices_for_quantum_joining': False,
            'sharp_bell_pair_sharing_error': '(d-1)*(k-1)/(d*k)',
            'sharp_bound_applies_to_all_complex_global_states': True,
            'attained_by_real_global_states': True,
            'unrestricted_extra_environment_removes_fixed_subsystem_obstruction': False,
            'distinct_logical_center_factors_allow_exact_real_bell_relations': True,
            'center_dimension_for_k_independent_d_level_bell_partners': 'd**k under the declared distinct-factor contract',
            'classical_star_with_common_A_marginal_always_has_extension': True,
            'all_classical_loopy_marginals_are_joinable': False,
            'rows': rows,
            'literature': {'url': 'https://arxiv.org/abs/1305.1342',
                           'checked': 'Definitions II.1/II.4, theorem III.8 and normalization F=Phi_plus/d',
                           'known_sharing_bound_claimed_original': False},
            'physical_protocol_for_merging_arbitrary_unknown_inputs_built': False,
            'quantum_necessity_or_full_SoCA_derived': False,
            'next': 'Define operational compatibility from actual preparation/history and allowed joint tests, before asking whether every compatible interface embeds functorially; do not infer local tomography or arbitrary relation-copying.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(IntegrationExtensionContractTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['automated_checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        target = Path(__file__).with_name('integration_extension_contract_results.json')
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
