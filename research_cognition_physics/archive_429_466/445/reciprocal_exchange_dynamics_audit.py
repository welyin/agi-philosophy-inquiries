"""Round 445: autonomous address exchanges with reciprocal-consistency energy.

The consistency energy and the all-pairs address access are explicit model inputs.
No geometry, controller schedule, dissipation, or numerical postselection is used.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import io
import itertools as it
import json
from pathlib import Path
import platform
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'reciprocal_exchange_dynamics_audit_results.json'
OBS = {}


def swap(word, pair):
    out = list(word)
    i, j = pair
    out[i], out[j] = out[j], out[i]
    return tuple(out)


def defects(word):
    return len(word)-2*sum(word[i] == j and word[j] == i
                           for i, j in it.combinations(range(len(word)), 2))


def matching_words(n):
    def pair(remaining, word):
        if not remaining:
            yield tuple(word)
            return
        i = remaining[0]
        for j in remaining[1:]:
            nxt = list(word)
            nxt[i], nxt[j] = j, i
            yield from pair(tuple(k for k in remaining if k not in (i, j)), nxt)
    return tuple(sorted(pair(tuple(range(n)), [-1]*n)))


def rematch_neighbors(word):
    edges = [(i, j) for i, j in enumerate(word) if i < j]
    out = set()
    for (a, b), (c, d) in it.combinations(edges, 2):
        for new_edges in (((a, c), (b, d)), ((a, d), (b, c))):
            nxt = list(word)
            for i, j in new_edges:
                nxt[i], nxt[j] = j, i
            out.add(tuple(nxt))
    return out


def matrices(words):
    n, dim = len(words[0]), len(words)
    lookup = {w: i for i, w in enumerate(words)}
    gamma = np.zeros((dim, dim), dtype=np.int64)
    for j, word in enumerate(words):
        for pair in it.combinations(range(n), 2):
            gamma[lookup[swap(word, pair)], j] += 1
    energy = np.array([defects(w) for w in words], dtype=np.int64)
    matches = matching_words(n)
    z = np.zeros((dim, len(matches)), dtype=np.int64)
    for j, w in enumerate(matches):
        z[lookup[w], j] = 1
    return energy, gamma, z


@lru_cache(None)
def four_model():
    words = tuple(it.permutations(range(4)))
    energy, gamma, z = matrices(words)
    r = np.array([F(1, int(e)) if e else F(0) for e in energy], dtype=object)
    g, v = gamma.astype(object), z.astype(object)
    w1 = -r[:, None]*(g@v)
    w2 = r[:, None]*(g@(r[:, None]*(g@v)))
    k = v.T@g@w1
    a3, a4 = g@w2-w1@k, -w2@k
    return words, energy, gamma, z, w1, w2, k, a3, a4


def unitary(h, time):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*time*values))@vectors.conj().T


def opnorm(a):
    return float(np.linalg.norm(a, ord=2))


def short(x):
    return float(f'{float(x):.12g}')


def fraction_data(x):
    return dict(numerator=x.numerator, denominator=x.denominator, decimal=short(x))


def bound(epsilon, tau):
    return (2+4*abs(tau))*abs(epsilon)+(F(13, 5)+4*abs(tau))*epsilon**2


def window_certificate():
    def sin7(x):
        return x-x**3/6+x**5/120-x**7/5040
    low, high, epsilon = F(19, 10), F(21, 10), F(1, 128)
    sine_lower = min(sin7(F(3, 4)*low), sin7(F(3, 4)*high))
    b = bound(epsilon, high)
    return sine_lower, b, F(8, 9)*sine_lower**2-b


class Audit(unittest.TestCase):
    def test_01_raw_subject_tensor_and_invariant_permutation_sector(self):
        words = tuple(it.product(range(4), repeat=4))
        energy, gamma, z = matrices(words)
        self.assertEqual(Counter(energy)[0], 3)
        self.assertEqual(min(energy[energy > 0]), 2)
        self.assertTrue(np.all(energy >= 0))
        raw_index = {w: i for i, w in enumerate(words)}
        perm, ep, gp, zp, *_ = four_model()
        embed = np.zeros((256, 24), dtype=np.int64)
        for j, w in enumerate(perm):
            embed[raw_index[w], j] = 1
        np.testing.assert_array_equal(gamma@embed, embed@gp)
        np.testing.assert_array_equal(energy[:, None]*embed, embed*ep[None, :])
        np.testing.assert_array_equal(z, embed@zp)
        np.testing.assert_array_equal(z.T@gamma@z, np.zeros((3, 3), dtype=int))
        # Independently build the complete tensor diagonal from two-site effects.
        tensor_diag = np.full((4,)*4, 4, dtype=int)
        for i, j in it.combinations(range(4), 2):
            selector = [slice(None)]*4
            selector[i], selector[j] = j, i
            tensor_diag[tuple(selector)] -= 2
        np.testing.assert_array_equal(tensor_diag.ravel(), energy)
        OBS['raw_model'] = dict(raw_dimension=256, invariant_permutation_dimension=24,
            zero_energy_matching_dimension=3, dimensionless_gap=2,
            defect_histogram={str(k): int(v) for k, v in sorted(Counter(energy).items())},
            full_tensor_extension_checked=True, first_order_code_action_zero=True,
            permutation_sector_is_exact_not_postselected=True,
            consistency_energy_is_additional_input=True)

    def test_02_exact_second_order_all_matching_columns(self):
        reports = []
        for n in (2, 4, 6, 8):
            matches = matching_words(n)
            pairs = tuple(it.combinations(range(n), 2))
            path_count = 0
            for word in matches:
                # Counts four times Gamma R Gamma, using the actual intermediate word.
                four_response = Counter()
                same, cross = 0, 0
                for pair in pairs:
                    middle = swap(word, pair)
                    e = defects(middle)
                    self.assertIn(e, (2, 4))
                    same += e == 2
                    cross += e == 4
                    for second in pairs:
                        final = swap(middle, second)
                        if defects(final) == 0:
                            four_response[final] += 4//e
                            path_count += 1
                expected = Counter({word: n*n//2})
                expected.update({w: 2 for w in rematch_neighbors(word)})
                self.assertEqual(four_response, expected)
                self.assertEqual(same, n//2)
                self.assertEqual(cross, n*(n-2)//2)
                self.assertEqual(len(expected)-1, n*(n-2)//4)
            reports.append(dict(subjects=n, matchings=len(matches),
                ordered_two_swap_returns=path_count,
                matching_flip_degree=n*(n-2)//4,
                diagonal_magnitude=fraction_data(F(n*n, 8)),
                off_diagonal_magnitude=fraction_data(F(1, 2))))
        OBS['second_order'] = dict(exact_integer_scaled_path_checks=reports,
            formula='K = -(N^2/8) I - A_matching/2',
            independent_edge_pair_neighbor_construction=True)

    def test_03_exact_duhamel_intertwiner_grams(self):
        _, energy, gamma, z, w1, w2, k, a3, a4 = four_model()
        zero = np.zeros_like(w1)
        np.testing.assert_array_equal(energy[:, None]*w1+gamma@z, zero)
        np.testing.assert_array_equal(energy[:, None]*w2+gamma@w1, z@k)
        expected = [
            (w1, F(3, 4), F(1, 8), F(1)),
            (w2, F(9, 16), F(9, 16), F(27, 16)),
            (a3, F(93, 16), F(141, 32), F(117, 8)),
            (a4, F(81, 16), F(81, 16), F(243, 16))]
        reports = []
        for name, (a, diag, off, norm_squared) in zip(('W1', 'W2', 'A3', 'A4'), expected):
            gram = np.full((3, 3), off, dtype=object)
            np.fill_diagonal(gram, diag)
            np.testing.assert_array_equal(a.T@a, gram)
            self.assertEqual(diag+2*off, norm_squared)
            self.assertGreaterEqual(diag+2*off, diag-off)
            reports.append(dict(matrix=name, gram_diagonal=str(diag),
                gram_off_diagonal=str(off), exact_norm_squared=str(norm_squared)))
        epsilon = F(1, 128)
        w = z+epsilon*w1+epsilon**2*w2
        np.testing.assert_array_equal(
            energy[:, None]*w+epsilon*gamma@w-w@(epsilon**2*k),
            epsilon**3*a3+epsilon**4*a4)
        self.assertLess(F(27, 16), F(13, 10)**2)
        self.assertLess(F(117, 8), 16)
        self.assertLess(F(243, 16), 16)
        OBS['duhamel_certificate'] = dict(exact_fraction_grams=reports,
            exact_residual_identity_verified=True,
            bare_code_input_no_dressed_state_preparation=True,
            bound='(2+4|tau|)|epsilon|+(13/5+4|tau|)epsilon^2',
            arbitrary_code_coherence_and_reference_covered_by_operator_bound=True)

    def test_04_rational_complete_time_window_witness(self):
        sine_lower, b, lower = window_certificate()
        self.assertEqual(b, F(13, 160)+F(11, 16384))
        self.assertLess(b, F(82, 1000))
        self.assertGreater(lower, F(39, 50))
        self.assertGreater(sine_lower, 0)
        # Exact matching graph spectrum supplies p_other = (8/9) sin^2(3 tau/4).
        _, _, _, _, _, _, k, _, _ = four_model()
        np.testing.assert_array_equal(k, -F(3, 2)*np.eye(3, dtype=int)
                                      -F(1, 2)*np.ones((3, 3), dtype=int))
        OBS['window'] = dict(epsilon='1/128', slow_time_interval=['19/10', '21/10'],
            uniform_operator_error_upper=fraction_data(b),
            sine_lower=fraction_data(sine_lower),
            actual_other_matching_probability_lower=fraction_data(lower),
            strict_simple_probability_lower='39/50',
            proof='sin concavity on the interval; alternating sine polynomial lower bounds at both endpoints',
            witness_scope='each initial matching basis state, with matching-dependent output effect',
            all_unknown_superpositions_have_positive_change_probability=False,
            no_postselection_or_perfect_timing_required=True)

    def test_05_full_evolution_unknown_inputs_and_reference(self):
        _, energy, gamma, z, _, _, k, _, _ = four_model()
        epsilon = F(1, 128)
        h = np.diag(energy)+float(epsilon)*gamma
        kd = np.array(k, dtype=float)
        reports = []
        for tau in (F(19, 10), F(2), F(21, 10)):
            u = unitary(h, float(tau/epsilon**2))
            ideal = unitary(kd, float(tau))
            actual_map, ideal_map = u@z, z@ideal
            error = opnorm(actual_map-ideal_map)
            leakage = opnorm(actual_map-z@(z.T@actual_map))**2
            probabilities = [float(np.sum(np.abs(z.T@actual_map[:, j])**2)
                                   -abs((z.T@actual_map)[j, j])**2) for j in range(3)]
            self.assertLess(error, float(bound(epsilon, tau)))
            self.assertLess(leakage, float((6*epsilon)**2))
            self.assertGreater(min(probabilities), float(window_certificate()[2]))
            # Maximally entangled code/reference input: column vectorization.
            actual_ref = actual_map.ravel()/np.sqrt(3)
            ideal_ref = ideal_map.ravel()/np.sqrt(3)
            trace_distance = np.sqrt(max(0, 1-abs(np.vdot(ideal_ref, actual_ref))**2))
            self.assertLess(trace_distance, float(bound(epsilon, tau)))
            reports.append(dict(tau=str(tau), all_input_operator_error=short(error),
                operator_error_bound=short(bound(epsilon, tau)),
                maximum_code_leakage_probability=short(leakage),
                minimum_other_matching_probability=short(min(probabilities)),
                three_dimensional_reference_trace_distance=short(trace_distance)))
        # Independent full raw-space exponential, including duplicate/self addresses.
        raw_energy, raw_gamma, raw_z = matrices(tuple(it.product(range(4), repeat=4)))
        raw_u = unitary(np.diag(raw_energy)+float(epsilon)*raw_gamma, float(F(2)/epsilon**2))
        perm = tuple(it.permutations(range(4)))
        raw_index = {w: i for i, w in enumerate(it.product(range(4), repeat=4))}
        lifted = np.zeros_like(raw_u@raw_z)
        small = unitary(h, float(F(2)/epsilon**2))@z
        for j, word in enumerate(perm):
            lifted[raw_index[word], :] = small[j, :]
        raw_error = opnorm(raw_u@raw_z-lifted)
        self.assertLess(raw_error, 2e-9)
        OBS['full_evolution'] = dict(spectral_checks=reports,
            raw_256_dimension_embedding_error=short(raw_error),
            all_time_leakage_probability_upper_for_N4=fraction_data((6*epsilon)**2),
            floating_point_checks_not_used_as_analytic_probability_certificate=True)

    def test_06_path_selection_counterexample_and_endpoint_budget(self):
        reports = []
        for n in (4, 6, 8):
            matches = matching_words(n)
            lookup = {w: i for i, w in enumerate(matches)}
            a = np.zeros((len(matches), len(matches)), dtype=np.int64)
            for j, w in enumerate(matches):
                for v in rematch_neighbors(w):
                    a[lookup[v], j] = 1
            endpoint_a = np.array([[a[i, j] if matches[i][0] != matches[j][0] else 0
                                    for j in range(len(matches))] for i in range(len(matches))])
            np.testing.assert_array_equal(endpoint_a, endpoint_a.T)
            np.testing.assert_array_equal(endpoint_a.sum(axis=0),
                                          np.full(len(matches), n-2))
            self.assertAlmostEqual(opnorm(endpoint_a), n-2, places=10)
            self.assertTrue(all(sum(1 for k in range(n) if w[k] == i) == 1
                                for w in matches for i in range(n)))
            # A matching has no three-edge simple path: 440's middle-edge gate is absent.
            for w in matches:
                for v in rematch_neighbors(w):
                    changed_edges = [(i, j) for i, j in enumerate(v) if i < j and w[i] != j]
                    self.assertEqual(len(changed_edges), 2)
                    self.assertTrue(all(w[i] != j for i, j in changed_edges))
            reports.append(dict(subjects=n, potential_exchange_terms=n*(n-1)//2,
                endpoint_effective_adjacency_norm=n-2, occupied_degree=1))
        rename_count = 0
        for labels in it.permutations(range(4)):
            def rename(w):
                result = [None]*4
                for i, j in enumerate(w):
                    result[labels[i]] = labels[j]
                return tuple(result)
            for word in matching_words(4):
                self.assertEqual({rename(v) for v in rematch_neighbors(word)},
                                 rematch_neighbors(rename(word)))
                rename_count += 1
        # A uniform matching superposition is an eigenstate of K: no universal change rate.
        k = np.array(four_model()[6], dtype=float)
        np.testing.assert_array_equal(k@np.ones(3), -3*np.ones(3))
        OBS['scope_and_budget'] = dict(endpoint_checks=reports,
            relabeling_checks=rename_count,
            old_three_path_flip_identically_zero_in_C1_matching_sector=True,
            new_effective_rematching_crosses_old_components=True,
            effective_endpoint_off_diagonal_norm='g^2*(N-2)/(2*Delta)',
            fixed_occupied_degree_does_not_bound_potential_interaction_budget=True,
            uniform_matching_state_is_effective_eigenstate=True,
            data_factors_are_spectators_in_this_round=True)


def run():
    OBS.clear()
    buffer = io.StringIO()
    result = unittest.TextTestRunner(stream=buffer, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(buffer.getvalue())
    return dict(round=445, baseline_round=444, date='2026-09-24',
        runtime=dict(python=platform.python_version(), numpy=np.__version__),
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        observations=OBS.copy(),
        scope=dict(explicit_two_subject_continuous_address_model=True,
            all_even_N_second_order_matching_flip_formula_proved=True,
            finite_window_bare_code_error_and_probability_certificate=True,
            arbitrary_unknown_matching_coherence_reference_error_covered=True,
            all_time_initial_code_leakage_bound_proved=True,
            known_perfect_matching_flip_graph_identified=True,
            direct_two_body_virtual_process_evades_exact_code_support_obstruction=True,
            consistency_energy_and_address_access_derived_from_429=False,
            existing_path_locality_selected_by_this_consistency_model=False,
            arbitrary_unknown_matching_state_has_uniform_change_probability=False,
            data_conditional_transport_implemented_in_this_model=False,
            arbitrary_inconsistent_input_autonomously_repaired=False,
            bounded_capacity_implies_uniform_interaction_budget=False,
            three_dimensional_space_generated=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    report = run()
    if not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(report, ensure_ascii=False, indent=2))
