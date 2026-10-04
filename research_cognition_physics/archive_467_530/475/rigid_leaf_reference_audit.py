"""Round 475: an all-time S4 selection rule for equal-time leaf distances.

Scientific baseline 472. No import of 473 or 474. Preparation, fixed six-tree
sector and later data-only distance instruments remain explicit model inputs.
"""
import argparse
from fractions import Fraction
from functools import lru_cache, reduce
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import branching_tree_distance_audit as trees_core

TARGET = Path(__file__).with_name('rigid_leaf_reference_audit_results.json')
OBS = {}
I = np.eye(2, dtype=complex)
P = (np.array([[0, 1], [1, 0]], dtype=complex),
     np.array([[0, -1j], [1j, 0]], dtype=complex),
     np.diag([1, -1]).astype(complex))
LEAVES = (0, 3, 4, 5)


def tensor(xs):
    return reduce(np.kron, xs)


def pauli_at(v, a):
    return tensor([P[a] if k == v else I for k in range(6)])


def sign(p):
    return (-1)**sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p)))


def dot(i, j):
    return sum(pauli_at(i, a) @ pauli_at(j, a) for a in range(3))


def chirality(i, j, k):
    return sum(sign(p)*pauli_at(i, p[0]) @ pauli_at(j, p[1]) @ pauli_at(k, p[2])
               for p in itertools.permutations(range(3)))


def cycles(p):
    seen = set()
    sizes = []
    for i in range(len(p)):
        if i in seen:
            continue
        j = i
        size = 0
        while j not in seen:
            seen.add(j)
            size += 1
            j = p[j]
        sizes.append(size)
    return tuple(sorted(sizes))


def char22(p):
    return {(1, 1, 1, 1): 2, (1, 1, 2): 0, (2, 2): 2,
            (1, 3): -1, (4,): 0}[cycles(p)]


@lru_cache(None)
def system():
    trees, _, h, f, _ = trees_core.six_vertex_sector()
    return trees, np.asarray(h, dtype=np.int64), np.asarray(f, dtype=np.int64)


@lru_cache(None)
def permutations():
    trees, _, _ = system()
    out = []
    for p in itertools.permutations(range(4)):
        full = list(range(6))
        for j in range(4):
            full[LEAVES[j]] = LEAVES[p[j]]
        dm = []
        for x in range(64):
            bits = [(x >> (5-v)) & 1 for v in range(6)]
            y = sum(bits[v] << (5-full[v]) for v in range(6))
            dm.append(y)
        gm = []
        for tree in trees:
            target = frozenset(tuple(sorted((full[a], full[b]))) for a, b in tree)
            gm.append(trees.index(target))
        out.append((p, np.array(dm), np.array(gm), char22(p)))
    return out


def source_numerators():
    return [tensor([I, b, I, I+P[0], I+P[1], I+P[2]]) for b in (I,)+P]


def twirled_numerators():
    # Each equals 384 times the relevant density or affine coefficient.
    return [6*np.eye(64)+chirality(3, 4, 5),
            2*dot(1, 3)+chirality(1, 4, 5),
            2*dot(1, 4)-chirality(1, 3, 5),
            2*dot(1, 5)+chirality(1, 3, 4)]


def eta(r, nonorthogonal=False):
    factors = [I, I+sum(r[a]*P[a] for a in range(3)), I,
               I+P[0], I+P[0 if nonorthogonal else 1], I+P[2]]
    return tensor(factors)/64


def distance(a, b):
    return np.array([len(trees_core.path_between(t, 6, a, b))-1 for t in system()[0]])


def gauss_integer(a):
    return bool(np.array_equal(a.real, np.rint(a.real)) and
                np.array_equal(a.imag, np.rint(a.imag)))


def project_graph_invariant(a):
    return sum(a[np.ix_(gm, gm)] for _, _, gm, _ in permutations())/24


class Audit(unittest.TestCase):
    def test_01_exact_symmetries_and_characters(self):
        _, h, f = system()
        fg = np.kron(np.eye(64, dtype=np.int64), f)
        v = h-fg
        for _, dm, gm, _ in permutations():
            total = np.array([6*d+g for d in dm for g in gm])
            self.assertTrue(np.array_equal(fg[np.ix_(total, total)], fg))
            self.assertTrue(np.array_equal(v[np.ix_(total, total)], v))
        for a in range(3):
            collective = np.kron(sum(pauli_at(v, a) for v in range(6)), np.eye(6))
            self.assertTrue(np.array_equal(h@collective, collective@h))
        records = []
        inner_products = [0, 0, 0]
        for p, _, _, character in permutations():
            fixed = lambda q: sum(q[j] == j for j in range(4))
            sq = tuple(p[p[j]] for j in range(4))
            cube = tuple(p[sq[j]] for j in range(4))
            chi = fixed(p)
            exterior2 = (chi**2-fixed(sq))//2
            exterior3 = (chi**3-3*chi*fixed(sq)+2*fixed(cube))//6
            cs = (chi, exterior2, exterior3)
            for j in range(3):
                inner_products[j] += character*cs[j]
            records.append(dict(cycle=list(cycles(p)), character22=character,
                                natural_and_exterior_characters=list(cs)))
        self.assertEqual(inner_products, [0, 0, 0])
        self.assertEqual(sum(x[3]**2 for x in permutations()), 24)
        # Group-algebra convolution proves E_22**2=E_22 without assuming
        # an external character table theorem for this finite certificate.
        for q, _, _, cq in permutations():
            coefficient = 0
            for p, _, _, cp in permutations():
                inverse = tuple(p.index(j) for j in range(4))
                relative = tuple(inverse[q[j]] for j in range(4))
                coefficient += cp*char22(relative)
            self.assertEqual(coefficient, 12*cq)
        self.assertEqual(sum(x[3] for x in permutations()), 0)
        OBS['exact_symmetries'] = dict(permutations=24, dimension=384,
            flip_and_exchange_covariant_separately=True,
            collective_SU2_generators_commute_exactly=True,
            character_inner_products_24_times=inner_products,
            character_records=records, central_projection_convolution_exact=True,
            arbitrary_real_kappa_and_J=True)

    def test_02_independent_exact_twirl(self):
        collective = [sum(pauli_at(v, a) for v in range(6)) for a in range(3)]
        def four_casimir(x):
            return sum(a@(a@x-x@a)-(a@x-x@a)@a for a in collective)
        # BEFORE computing: row/column absolute sums of A_mu are <= 6.
        # Each four-Casimir component costs <= 3*(2*6)**2 = 432.
        # Raw Gaussian components <= 2; 128 also covers all intermediate
        # products, sums and final factor-six cross multiplication.
        prior_bound = 128*math.prod(4*ell*(ell+1)+432 for ell in range(1, 5))
        self.assertLess(prior_bound, 2**53)
        max_intermediate = 0
        for raw, predicted in zip(source_numerators(), twirled_numerators()):
            value = raw.copy()
            denominator = 1
            # The support has at most four Pauli vectors: only spins l=0..4.
            for ell in range(1, 5):
                factor = 4*ell*(ell+1)
                value = factor*value-four_casimir(value)
                denominator *= factor
                self.assertTrue(gauss_integer(value))
                max_intermediate = max(max_intermediate,
                    int(max(abs(value.real).max(), abs(value.imag).max())))
            self.assertTrue(np.array_equal(6*value, denominator*predicted))
            for a in collective:
                self.assertTrue(np.array_equal(a@predicted, predicted@a))
        self.assertLess(max_intermediate, 2**53)
        OBS['exact_twirl'] = dict(affine_source_coefficients=4,
            density_denominator=384, casimir_denominator=denominator,
            max_gaussian_integer_component=max_intermediate,
            a_priori_integer_component_bound=prior_bound,
            support_spin_cutoff=4, predicted_pair_and_chirality_formula_exact=True,
            twirl_used_as_expectation_identity_not_a_free_physical_preparation=True)

    def test_03_selection_rule_certificates(self):
        sizes = []
        for n in twirled_numerators():
            weighted = sum(c*n[np.ix_(dm, dm)] for _, dm, _, c in permutations())
            self.assertTrue(np.array_equal(weighted, np.zeros((64, 64))))
            sizes.append(int(np.count_nonzero(n)))
        qvectors = []
        for a, b in itertools.combinations(LEAVES, 2):
            q = 3*distance(a, b)-8
            weighted = sum(c*q[gm] for _, _, gm, c in permutations())
            self.assertTrue(np.array_equal(weighted, 12*q))
            qvectors.append(q)
        self.assertEqual(np.linalg.matrix_rank(qvectors), 2)
        # Exact trivial projection of the three matching indicators.
        matching = [np.diag((distance(0, b) == 2).astype(int)) for b in (3, 4, 5)]
        self.assertTrue(np.array_equal(sum(matching), np.eye(6)))
        for m in matching:
            self.assertTrue(np.array_equal(24*project_graph_invariant(m), 8*np.eye(6)))
        OBS['selection_rule'] = dict(affine_coefficients_central_projection_zero=4,
            centered_leaf_distances_in_22=6, centered_module_dimension=2,
            unnormalized_central_projector_eigenvalue=12,
            twirl_numerator_nonzero_entries=sizes,
            graph_S4_invariant_initial_state_may_be_mixed=True,
            graph_reference_must_be_jointly_invariant_not_only_its_marginal=True,
            all_time_proof_is_symmetry_not_Taylor_extrapolation=True)

    def test_04_full_unitary_and_unknown_reference(self):
        _, h, f = system()
        fg = np.kron(np.eye(64), f)
        vv = h-fg
        sg = np.ones((6, 6))/6
        rows = []
        worst = 0.0
        for kap, coupling in ((1., 1.), (-.6, 1.4), (1.25, -.3)):
            eig, vec = np.linalg.eigh(kap*fg+coupling*vv)
            for t in (.17, .63):
                u = (vec*np.exp(-1j*t*eig))@vec.conj().T
                for r in ((0., 0., 0.), (.2, -.3, .4), (1., 0., 0.)):
                    state = np.kron(eta(r), sg)
                    evolved = u@state@u.conj().T
                    population = np.diag(evolved).real.reshape(64, 6).sum(axis=0)
                    errors = [abs(population@distance(a, b)-8/3)
                              for a, b in itertools.combinations(LEAVES, 2)]
                    worst = max(worst, *errors)
                    rows.append(dict(kappa=kap, J=coupling, time=t,
                                     source=list(r), max_mean_error=max(errors)))
        self.assertLess(worst, 3e-12)
        eig, vec = np.linalg.eigh(h)
        u = (vec*np.exp(-1j*.37*eig))@vec.T
        # Independent Bell source-reference input, mixed leaf0/internal2,
        # pure +X/+Y/+Z reference leaves, graph |s>. No copying of source.
        basis = np.eye(2)
        plus = (basis[:, 0]+basis[:, 1])/np.sqrt(2)
        plus_y = (basis[:, 0]+1j*basis[:, 1])/np.sqrt(2)
        means_r = [np.zeros((2, 2), dtype=complex) for _ in range(3)]
        masks = [(distance(0, b) == 2) for b in (3, 4, 5)]
        for a, b in itertools.product(range(2), repeat=2):
            columns = [np.kron(tensor([basis[:, a], basis[:, q], basis[:, b],
                          plus, plus_y, basis[:, 0]]), np.ones(6)/np.sqrt(6))/np.sqrt(2)
                       for q in range(2)]
            evolved = (u@np.column_stack(columns)).reshape(64, 6, 2)
            for j, mask in enumerate(masks):
                w = evolved[:, mask, :].reshape(-1, 2)
                means_r[j] += w.T@w.conj()/4
        reference_error = max(np.max(abs(m-np.eye(2)/6)) for m in means_r)
        self.assertLess(reference_error, 3e-12)
        # Complete graph-reference extension: invariant projection of the
        # effective graph effect equals I/3, tested independently of source R.
        fixed = eta((.2, -.3, .4))
        graph_effect_errors = []
        for mask in masks:
            effect = u.conj().T@np.kron(np.eye(64), np.diag(mask.astype(int)))@u
            eg = np.einsum('ba,aibj->ij', fixed, effect.reshape(64, 6, 64, 6))
            graph_effect_errors.append(float(np.max(abs(project_graph_invariant(eg)-np.eye(6)/3))))
        self.assertLess(max(graph_effect_errors), 3e-12)
        # A uniform classical graph label stored in R has invariant graph
        # marginal, but already at t=0 its reference distance moment differs.
        marginal_only_error = max(abs(Fraction(int(d), 6)-Fraction(4, 9))
                                  for d in distance(0, 3))
        self.assertEqual(marginal_only_error, Fraction(1, 9))
        OBS['numerical_supplement'] = dict(full_unitary_samples=rows,
            marginal_only_graph_invariance_t0_reference_counterexample=str(marginal_only_error),
            maximum_leaf_mean_error=worst,
            Bell_source_reference_matching_moment_error=float(reference_error),
            graph_reference_effect_errors=graph_effect_errors,
            numerical_witnesses_do_not_replace_exact_selection_proof=True)

    def test_05_nonorthogonal_preparation_counterexample(self):
        _, h, _ = system()
        raw = tensor([I, I, I, I+P[0], I+P[0], I+P[2]]).real.astype(np.int64)
        initial_numerator = np.kron(raw, np.ones((6, 6), dtype=np.int64))
        op = np.kron(np.eye(64, dtype=np.int64), np.diag(distance(3, 4)))
        derivative_values = []
        for k in range(1, 5):
            op = h@op-op@h
            numerator = int(np.einsum('ij,ji->', initial_numerator, op))
            if k % 2:
                self.assertEqual(numerator, 0)
            derivative_values.append(str(Fraction(numerator, 384)*((-1)**(k//2))))
        self.assertEqual(derivative_values[:3], ['0', '0', '0'])
        self.assertEqual(derivative_values[3], '-44/3')
        # Bound all integer products in this four-step exact trace computation.
        conservative_integer_bound = 384**2*4*3*18**4
        self.assertLess(conservative_integer_bound, 2**63)
        self.assertFalse(np.array_equal(h@np.kron(np.eye(64), np.diag(distance(3, 4))),
                                       np.kron(np.eye(64), np.diag(distance(3, 4)))@h))
        OBS['preparation_boundary'] = dict(replacement_reference_axes=['X', 'X', 'Z'],
            source_maximally_mixed=True, kappa=1, J=1,
            d34_first_four_derivatives=derivative_values,
            fourth_order_coefficient='-11/18', remainder_order=6,
            real_state_H_observable_give_even_scalar_expectation=True,
            exact_integer_trace_bound=conservative_integer_bound,
            distance_not_a_conserved_operator=True,
            no_universal_rigid_reference_claim=True)

    def test_06_joint_shape_distribution_and_geometry_boundary(self):
        ds = np.array([distance(0, b) for b in (3, 4, 5)])
        mean = Fraction(8, 3)
        self.assertTrue(np.array_equal(ds.sum(axis=0), np.full(6, 8)))
        covariance = []
        for a in ds:
            row = []
            for b in ds:
                value = sum((Fraction(int(x))-mean)*(Fraction(int(y))-mean)
                            for x, y in zip(a, b))/6
                row.append(str(value))
            covariance.append(row)
        self.assertEqual(covariance, [['2/9', '-1/9', '-1/9'],
                                     ['-1/9', '2/9', '-1/9'],
                                     ['-1/9', '-1/9', '2/9']])
        # If these means themselves are asserted to be Euclidean lengths,
        # a regular four-leaf simplex side 8/3 has squared circumradius 8/3.
        radius_squared = Fraction(3, 8)*mean**2
        source_radius_squared = Fraction(3, 2)**2
        self.assertEqual(radius_squared, Fraction(8, 3))
        self.assertGreater(radius_squared, source_radius_squared)
        OBS['scope_and_readout'] = dict(matching_class_probabilities=['1/3']*3,
            every_leaf_distance_distribution={'2': '1/3', '3': '2/3'},
            variance='2/9', representative_distance_covariances=covariance,
            finite_joint_equal_time_shape_distribution_fixed=True,
            single_run_shape_rigidity_or_multitime_stationarity_claimed=False,
            graph_coherences_not_reconstructed_by_six_distance_means=True,
            initial_source_leaf_mean='3/2', leaf_mean='8/3',
            Euclidean_leaf_circumradius_squared=str(radius_squared),
            proposed_equal_source_radius_squared=str(source_radius_squared),
            mean_lengths_do_not_form_one_Euclidean_anchor_configuration=True,
            data_only_estimation_available_conditionally_from_round=472,
            independent_copies_fresh_probes_storage_timing_reference_all_counted=True,
            no_sharp_free_graph_measurement_claim=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=475, baseline_round=472, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(all_time_all_real_couplings_conditional_selection_rule=True,
            no_dependency_on_round473_or_round474=True,
            source_reference_extension_and_jointly_invariant_graph_reference_extension=True,
            preparation_symmetry_is_explicit_model_input=True,
            no_lossless_new_reference_preparation_claim=True,
            classical_rigid_frame_derived=False, displacement_group_derived=False,
            three_dimensional_space_derived=False, full_GR_goal_completed=False,
            phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert result == json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
