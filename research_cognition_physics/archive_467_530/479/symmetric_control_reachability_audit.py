"""Round 479: a two-coefficient symmetry restriction on fixed-source readouts.

Scientific baseline 475. The six-tree sector, specified reference preparation,
fixed covariant CPTP maps, and data-only readout permissions remain inputs.
This does not classify physical dimension or all possible internal controls.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import rigid_leaf_reference_audit as old
import sequential_tree_distance_readout_audit as reader

TARGET = Path(__file__).with_name('symmetric_control_reachability_audit_results.json')
OBS = {}
K = np.array([[0, 1, -1], [-1, 0, 1], [1, -1, 0]], dtype=np.int64)
ROOT = np.column_stack([-np.ones(3, dtype=np.int64), np.eye(3, dtype=np.int64)])
I = np.eye(2, dtype=complex)


def short(x):
    return float(f'{float(x):.11g}')


def exact_rank(matrix):
    a = [[F(int(v)) for v in row] for row in matrix]
    row = 0
    for col in range(len(a[0])):
        pivot = next((j for j in range(row, len(a)) if a[j][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        scale = a[row][col]
        a[row] = [v/scale for v in a[row]]
        for j in range(row+1, len(a)):
            scale = a[j][col]
            if scale:
                a[j] = [v-scale*w for v, w in zip(a[j], a[row])]
        row += 1
        if row == len(a):
            break
    return row


def wedge_matrix(p, degree):
    basis = list(itertools.combinations(range(4), degree))
    out = np.zeros((len(basis), len(basis)), dtype=np.int64)
    for j, support in enumerate(basis):
        image = [p[a] for a in support]
        parity = (-1)**sum(image[a] > image[b] for a in range(degree)
                          for b in range(a+1, degree))
        out[basis.index(tuple(sorted(image))), j] = parity
    return out


def qgraph():
    return np.array([old.distance(1, j)-old.distance(1, 0)
                     for j in (3, 4, 5)], dtype=np.int64)


@lru_cache(None)
def representations():
    q = qgraph()
    out = []
    for p, dm, gm, _ in old.permutations():
        numerator = q[:, gm]@q.T@(4*np.eye(3, dtype=np.int64)-np.ones((3, 3), dtype=np.int64))
        assert np.all(numerator % 8 == 0)
        rq = numerator//8
        assert np.array_equal(rq@q, q[:, gm])
        out.append((p, dm, gm, rq))
    return out


def intertwiner_constraints(degree):
    dim = math.comb(4, degree)
    rows = []
    for p, _, _, rq in representations():
        v = wedge_matrix(p, degree)
        for i in range(3):
            for j in range(dim):
                row = np.zeros((3, dim), dtype=np.int64)
                row[i, :] += v[:, j]
                row[:, j] -= rq[i, :]
                rows.append(row.ravel())
    return rows


def incidence():
    out = np.zeros((4, 6), dtype=np.int64)
    for j, (a, b) in enumerate(itertools.combinations(range(4), 2)):
        out[a, j], out[b, j] = 1, -1
    return out


def unitary(h):
    v, w = np.linalg.eigh(h)
    return (w*np.exp(-1j*v))@w.conj().T


@lru_cache(None)
def s12():
    mapping = []
    for x in range(64):
        bits = [(x >> (5-j)) & 1 for j in range(6)]
        bits[1], bits[2] = bits[2], bits[1]
        mapping.append(sum(bits[j] << (5-j) for j in range(6)))
    return np.kron(np.eye(64)[mapping], np.eye(6))


def schedule(items):
    h = old.system()[1]
    w = np.eye(384, dtype=complex)
    for duration, area in items:
        w = unitary(duration*h+area*s12())@w
    return w


def response(w, graph=None):
    if graph is None:
        graph = np.ones((6, 6))/6
    rows = []
    for q in qgraph():
        effect = w.conj().T@(np.tile(q, 64)[:, None]*w)
        data_effect = np.einsum('ji,aibj->ab', graph, effect.reshape(64, 6, 64, 6))
        rows.append([np.trace(data_effect@n).real/64 for n in old.source_numerators()])
    rows = np.asarray(rows)
    return rows[:, 0], rows[:, 1:]


def fit(a):
    alpha = float(np.trace(a)/3)
    beta = float(np.sum(a*K)/6)
    return alpha, beta, float(np.linalg.norm(a-alpha*np.eye(3)-beta*K))


def output_graph(w, r):
    state = w@np.kron(old.eta(r), np.ones((6, 6))/6)@w.conj().T
    return state.reshape(64, 6, 64, 6).trace(axis1=0, axis2=2)


def exact_channel_coefficients():
    h = old.system()[1]
    raw = [np.kron(n, np.ones((6, 6))) for n in old.source_numerators()]
    assert 4*384**2*18**3 < 2**53
    out = []
    for degree in range(4):
        rows = []
        for q in qgraph():
            op = np.diag(np.tile(q, 64))
            for _ in range(degree):
                op = h@op-op@h
            row = []
            for n in raw:
                v = np.einsum('ij,ji->', n, op)*(1j**degree)
                assert v.imag == 0 and v.real == round(v.real)
                row.append(F(int(v.real), 384*math.factorial(degree)))
            rows.append(row)
        out.append(rows)
    return out


class Audit(unittest.TestCase):
    def close(self, a, b, tol=3e-11):
        self.assertLess(float(np.linalg.norm(np.asarray(a)-np.asarray(b))), tol)

    def test_01_complete_finite_intertwiner_certificate(self):
        ranks, dimensions = [], []
        maps = (None, ROOT, ROOT@incidence(), None)
        for degree in range(4):
            dim = math.comb(4, degree)
            rank = exact_rank(intertwiner_constraints(degree))
            ranks.append(rank)
            dimensions.append(3*dim-rank)
            if maps[degree] is not None:
                for p, _, _, rq in representations():
                    self.assertTrue(np.array_equal(maps[degree]@wedge_matrix(p, degree),
                                                   rq@maps[degree]))
        self.assertEqual(ranks, [3, 11, 17, 12])
        self.assertEqual(dimensions, [0, 1, 1, 0])
        for p, _, _, rq in representations():
            self.assertEqual(int(np.trace(rq)), sum(p[j] == j for j in range(4))-1)
        OBS['finite_intertwiner_certificate'] = dict(
            permutation_count=24, input_modules=['trivial', 'V', 'wedge2(V)', 'wedge3(V)'],
            linear_equation_ranks=ranks, solution_dimensions=dimensions,
            rational_elimination_not_floating_rank=True,
            explicit_nonzero_maps=['ROOT', 'ROOT times incidence'],
            output_is_Q_difference_representation=True,
            all_intertwiners_classified_not_only_selected_Hamiltonians=True)

    def test_02_source_decomposition_and_plane(self):
        # 475: 384 T eta = N0 + rx Nx + ry Ny + rz Nz.
        pair_embedding = np.vstack([np.zeros((1, 3), dtype=int), np.eye(3, dtype=int)])
        chirality_embedding = np.zeros((6, 3), dtype=np.int64)
        pairs = list(itertools.combinations(range(4), 2))
        chirality_embedding[pairs.index((2, 3)), 0] = 1
        chirality_embedding[pairs.index((1, 3)), 1] = -1
        chirality_embedding[pairs.index((1, 2)), 2] = 1
        self.assertTrue(np.array_equal(ROOT@pair_embedding, np.eye(3)))
        self.assertTrue(np.array_equal(ROOT@incidence()@chirality_embedding, -K))
        self.assertTrue(np.array_equal(K.T, -K))
        self.assertTrue(np.array_equal(K.T@K, 3*np.eye(3)-np.ones((3, 3))))
        self.assertTrue(np.array_equal(K@np.ones(3), np.zeros(3)))
        coeff = exact_channel_coefficients()
        zero = [[F(0)]*4 for _ in range(3)]
        self.assertEqual(coeff[0], zero)
        self.assertEqual(coeff[1], zero)
        self.assertEqual(coeff[2], [[F(0)]+[F(-2, 3) if i == j else F(0)
                                                     for j in range(3)] for i in range(3)])
        self.assertEqual(coeff[3], [[F(0)]+[F(int(K[i, j]), 9)
                                                     for j in range(3)] for i in range(3)])
        OBS['exact_response_restriction'] = dict(
            channel_form='Q(r)=a*r+b*K*r', K=K.tolist(), constant_offset_zero=True,
            free_H_a_leading='-2*t^2/3', free_H_b_leading='t^3/9',
            no_hidden_third_intertwiner=True,
            fixed_source_ex_identity='Q_y+Q_z=0',
            fixed_nonzero_source_plane='span(r,K*r)',
            source_parallel_111_at_most_one_dimension=True, zero_source_always_zero=True,
            full_input_response_may_still_have_rank_three=True,
            plane_dimension_bound_not_reachability_of_the_whole_plane=True)

    def test_03_actual_symmetric_controls_and_jacobian(self):
        trials = [[(.17, 0.)], [(.13, .27), (.21, -.19), (.08, .35)],
                  [(.31, -.4)], [(.07, .6), (.22, .1)]]
        rows = []
        for items in trials:
            w = schedule(items)
            bias, a = response(w)
            alpha, beta, residual = fit(a)
            self.close(bias, 0)
            self.assertLess(residual, 4e-13)
            values = a@np.array([1., 0., 0.])
            self.assertLess(abs(values[1]+values[2]), 4e-13)
            rows.append(dict(schedule=[list(x) for x in items], a=short(alpha), b=short(beta),
                             matrix_residual=short(residual), fixed_source_Q=[short(x) for x in values]))
        h = old.system()[1]
        t, u, theta = .13, .21, .27
        left, middle, right = unitary(u*h), unitary(theta*s12()), unitary(t*h)
        w = left@middle@right
        derivatives = [left@middle@(-1j*h@right),
                       (-1j*h@left)@middle@right,
                       left@(-1j*s12()@middle)@right]
        rho = np.kron(old.eta((1., 0., 0.)), np.ones((6, 6))/6)
        jacobian = np.array([[2*np.trace((np.tile(q, 64)[:, None]*dw)@rho@w.conj().T).real
                              for dw in derivatives] for q in qgraph()])
        singular = np.linalg.svd(jacobian, compute_uv=False)
        self.assertLess(float(singular[-1]), 4e-13)
        self.assertGreater(float(singular[1]), 1e-5)
        self.close(jacobian[1]+jacobian[2], 0, 5e-13)
        OBS['actual_control_witnesses'] = dict(samples=rows,
            derivative_parameters=['first_wait', 'second_wait', 'middle_exchange_area'],
            analytic_unitary_derivatives_not_finite_differences=True,
            Jacobian_singular_values=[short(x) for x in singular],
            numerical_rank_two_example_not_global_controllability_proof=True,
            internally_selected_pulses_are_extra_control_permissions=True)

    def test_04_general_covariant_channel_and_broken_symmetry(self):
        h = old.system()[1]
        # A graph-only potential is an explicit extra, symmetry-breaking input.
        trees, _, f = old.system()
        n0 = np.diag([int((0, 1) in tree) for tree in trees]).astype(np.int64)
        g = f+n0
        sg = np.ones((6, 6))/6
        wg = unitary(.2*g)
        p = np.diag(wg@sg@wg.conj().T).real
        values = qgraph()@p
        self.assertGreater(abs(values[1]+values[2]), .01)
        for q in qgraph():
            op = np.diag(q)
            nested = g@(g@op-op@g)-(g@op-op@g)@g
            self.assertEqual(F(-int(nested.sum()), 6), F(8, 3))
        # Average this CPTP map over the 24 conjugated leaf permutations.
        # This restores covariance and yields a=0=b (data are ignored).
        covariant_q = np.zeros(3)
        for _, _, gm, _ in representations():
            inverse = np.argsort(gm)
            w = wg[np.ix_(inverse, inverse)]
            covariant_q += qgraph()@np.diag(w@sg@w.conj().T).real/24
        self.close(covariant_q, 0)
        OBS['scope_of_covariance'] = dict(
            asymmetric_graph_generator='F+n_10, an extra graph potential',
            exact_Q_second_derivative=[str(F(8, 3))]*3,
            asymmetric_Q=[short(x) for x in values],
            fixed_ex_plane_violated_when_joint_leaf_covariance_dropped=True,
            collective_SU2_still_present_in_asymmetric_counterexample=True,
            S4_averaged_CPTP_example_Q=[short(x) for x in covariant_q],
            no_assertion_that_this_new_graph_control_is_available=True,
            no_no_go_for_all_reference_preserving_controls=True)

    def test_05_two_separate_unknown_reference_contracts(self):
        w = schedule([(.13, .27), (.21, -.19), (.08, .35)])
        _, a = response(w)
        alpha, beta, _ = fit(a)
        basis = np.eye(2)
        plus = (basis[:, 0]+basis[:, 1])/math.sqrt(2)
        plus_y = (basis[:, 0]+1j*basis[:, 1])/math.sqrt(2)
        moments = [np.zeros((2, 2), complex) for _ in range(3)]
        for left, right in itertools.product(range(2), repeat=2):
            columns = [np.kron(old.tensor([basis[:, left], basis[:, q], basis[:, right],
                        plus, plus_y, basis[:, 0]]), np.ones(6)/math.sqrt(6))/math.sqrt(2)
                       for q in range(2)]
            out = w@np.column_stack(columns)
            for i, q in enumerate(qgraph()):
                moments[i] += out.T@(np.tile(q, 64)[:, None]*out.conj())/4
        expected = [sum(a[i, j]*old.P[j].T/2 for j in range(3)) for i in range(3)]
        for got, want in zip(moments, expected):
            self.close(got, want)
        # Jointly S4-invariant G-R input can have branch-dependent a and b.
        _, amixed = response(w, np.eye(6)/6)
        am, bm, residual = fit(amixed)
        self.assertLess(residual, 4e-13)
        self.assertGreater(abs(alpha-am)+abs(beta-bm), 1e-5)
        correlated = (np.kron(old.eta((1., 0., 0.)), np.ones((6, 6))/6)
                      +np.kron(old.eta((-1., 0., 0.)), np.eye(6)/6))/2
        evolved = w@correlated@w.conj().T
        pg = np.diag(evolved).real.reshape(64, 6).sum(axis=0)
        correlated_q = qgraph()@pg
        self.close(correlated_q, (a[:, 0]-amixed[:, 0])/2)
        self.assertGreater(float(np.linalg.norm(correlated_q)), 1e-3)
        OBS['unknown_reference_contracts'] = dict(
            independent_graph_arbitrary_source_R_formula='Q_R=a*s_R+b*K*s_R; s_mu=Tr_source(sigma_mu rho_source_R)',
            Bell_source_R_largest_error=short(max(np.linalg.norm(g-e) for g, e in zip(moments, expected))),
            graph_R_joint_invariance_required=True,
            independent_source_joint_graph_R_formula='Q_R(r)=A_R*r+B_R*K*r',
            graph_R_branch_coefficients=[[short(alpha), short(beta)], [short(am), short(bm)]],
            graph_R_coefficients_are_reference_operators_not_one_universal_scalar_pair=True,
            arbitrary_initial_source_graph_reference_correlations_not_included=True,
            correlated_source_with_zero_marginal_Q=[short(x) for x in correlated_q],
            product_input_condition_cannot_be_replaced_by_its_two_marginals=True,
            CPTP_bound_does_not_claim_every_channel_preserves_unknown_quantum_information=True)

    def test_06_actual_data_readout_of_fixed_source_plane(self):
        w = schedule([(.13, .27), (.21, -.19), (.08, .35)])
        rho = output_graph(w, (1., 0., 0.))
        true = qgraph()@np.diag(rho).real
        probe = F(1, 65536)
        means = []
        for leaf in (0, 3, 4, 5):
            minus, plus, _ = reader.instrument((1, leaf), probe)
            signed = reader.apply_map((plus-minus)/float(probe), rho, g=6, r=1)
            means.append(float(np.trace(signed).real))
        estimated = means[0]-np.asarray(means[1:])
        bias = (math.expm1(50*float(probe))-50*float(probe))/float(probe)
        self.close(true[1]+true[2], 0)
        self.assertLess(abs(estimated[1]+estimated[2]), 4*bias)
        self.assertLess(float(np.max(abs(estimated-true))), 2*bias)
        OBS['actual_readout_and_goal_scope'] = dict(
            true_Q=[short(x) for x in true], actual_instrument_Q=[short(x) for x in estimated],
            probe_wait=str(probe), single_edge_bias_bound=short(bias),
            fixed_source_plane_record_error_bound='4 times single-edge error for Q_y+Q_z',
            actual_CP_instruments_from_round472=True,
            fixed_preparation_and_readout_endpoint=1,
            old_correlations_isolated_before_fresh_independent_probes=True,
            repeats_controls_clocks_storage_and_joint_instrument_errors_remain_inputs=True,
            no_new_copy_budget_optimization_or_single_shot_exact_coordinate_claim=True,
            legal_three_coordinate_quotient_not_disproved=True,
            only_fixed_preparation_and_covariant_control_family_restricted=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=479, baseline_round=475, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(exact_two_coefficient_covariant_response_classification=True,
            fixed_source_control_readouts_lie_in_at_most_a_plane=True,
            whole_plane_reachability_not_proved=True,
            source_preparation_and_joint_leaf_covariance_are_restrictive_inputs=True,
            asymmetric_controls_changed_preparations_or_correlations_not_excluded=True,
            arbitrary_controls_or_cognitive_axioms_not_refuted=True,
            full_GR_goal_completed=False, phase_closure_triggered=False))


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

