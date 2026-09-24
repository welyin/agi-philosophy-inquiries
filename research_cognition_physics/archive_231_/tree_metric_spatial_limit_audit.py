"""Round 467: tree metric limits and the sharp fluctuation cost of a square mean.

All algebraic certificates use Q(sqrt(2)) with Fraction coefficients. NumPy is
used only for the separate coherent graph/reference example. No new dynamics
or measured spatial metric is asserted. Default results are exclusive-create;
--dry-run is read-only and --check compares with the frozen result.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

TARGET = Path(__file__).with_name("tree_metric_spatial_limit_audit_results.json")
OBS = {}


@dataclass(frozen=True)
class Q2:
    p: F = F(0)
    q: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "p", F(self.p))
        object.__setattr__(self, "q", F(self.q))

    @staticmethod
    def cast(x):
        return x if isinstance(x, Q2) else Q2(x)

    def __add__(self, other):
        z = self.cast(other)
        return Q2(self.p + z.p, self.q + z.q)

    __radd__ = __add__

    def __neg__(self):
        return Q2(-self.p, -self.q)

    def __sub__(self, other):
        return self + (-self.cast(other))

    def __rsub__(self, other):
        return self.cast(other) - self

    def __mul__(self, other):
        z = self.cast(other)
        return Q2(self.p*z.p + 2*self.q*z.q, self.p*z.q + self.q*z.p)

    __rmul__ = __mul__

    def __truediv__(self, other):
        z = self.cast(other)
        den = z.p*z.p - 2*z.q*z.q
        if not den:
            raise ZeroDivisionError
        return self * Q2(z.p/den, -z.q/den)

    def sign(self):
        # Exact ordering; never infer an algebraic sign from floating point.
        p, q = self.p, self.q
        if not q:
            return (p > 0) - (p < 0)
        if not p:
            return (q > 0) - (q < 0)
        if p*q > 0:
            return (p > 0) - (p < 0)
        diff = p*p - 2*q*q
        return ((p > 0) - (p < 0))*((diff > 0) - (diff < 0))

    def __abs__(self):
        return self if self.sign() >= 0 else -self

    def __float__(self):
        return float(self.p) + math.sqrt(2)*float(self.q)

    def record(self):
        return {"rational": str(self.p), "sqrt2_coefficient": str(self.q),
                "approx": float(format(float(self), '.13g'))}


Z, ONE, R2 = Q2(), Q2(1), Q2(0, 1)
PAIRS = list(itertools.combinations(range(4), 2))
SIDES = ((0, 1), (0, 3), (1, 2), (2, 3))
SPLITS = (((0, 1), (2, 3)), ((0, 3), (1, 2)), ((0, 2), (1, 3)))


def maximum(items):
    best = items[0]
    for x in items[1:]:
        if (x-best).sign() > 0:
            best = x
    return best


def leaf_distances(edges):
    adj = {}
    for u, v, weight in edges:
        w = Q2.cast(weight)
        adj.setdefault(u, []).append((v, w))
        adj.setdefault(v, []).append((u, w))
    assert len(edges) == len(adj)-1
    result = {}
    for source in range(4):
        seen = {source: Z}
        todo = [source]
        while todo:
            u = todo.pop()
            for v, w in adj[u]:
                if v not in seen:
                    seen[v] = seen[u] + w
                    todo.append(v)
        assert len(seen) == len(adj)
        for target in range(source+1, 4):
            result[source, target] = seen[target]
    return result, adj


def quartet(split, leaf_weights, central):
    edges = [(4, 5, Q2.cast(central))]
    for center, leaves in zip((4, 5), split):
        edges += [(center, leaf, Q2.cast(leaf_weights[leaf])) for leaf in leaves]
    return leaf_distances(edges)[0]


def sums(d):
    return (d[0, 1]+d[2, 3], d[0, 3]+d[1, 2], d[0, 2]+d[1, 3])


def square(side=ONE):
    return {p: side if p in SIDES else R2*side for p in PAIRS}


def mean(ds, weights):
    return {p: sum((w*d[p] for w, d in zip(weights, ds)), Z) for p in PAIRS}


def avg_square_pair():
    a, b = (2-R2)/2, 2*(R2-1)
    return [quartet(split, [a]*4, b) for split in SPLITS[:2]], a, b


def unit_subdivision(split, m, k):
    edges, nxt = [], 6

    def path(u, v, length):
        nonlocal nxt
        for _ in range(length-1):
            edges.append((u, nxt, 1))
            u, nxt = nxt, nxt+1
        edges.append((u, v, 1))

    path(4, 5, k)
    for center, leaves in zip((4, 5), split):
        for leaf in leaves:
            path(center, leaf, m)
    return leaf_distances(edges)


class Audit(unittest.TestCase):
    def test_01_exact_tree_four_point_certificates(self):
        lo, hi = F(1414213562373095, 10**15), F(1414213562373096, 10**15)
        self.assertLess(lo*lo, 2)
        self.assertGreater(hi*hi, 2)
        count = 0
        for leaves in itertools.permutations([F(1, 7), F(2, 5), F(3, 4), F(7, 3)]):
            for split in SPLITS:
                for b in (Z, Q2(F(1, 3)), R2):
                    d = quartet(split, leaves, b)
                    ss = sums(d)
                    top = maximum(list(ss))
                    self.assertGreaterEqual(sum(x == top for x in ss), 2)
                    # General pointwise inequality used by the variance proof.
                    self.assertGreaterEqual((ss[0]+ss[1]+abs(ss[0]-ss[1])-2*ss[2]).sign(), 0)
                    count += 1
        OBS['four_point'] = dict(exact_quartets_checked=count,
            root2_bracket=[str(lo), str(hi)],
            mature_source='Buneman 1974, A Note on the Metric Properties of Trees',
            finite_checks_not_the_general_proof=True)

    def test_02_sharp_additive_and_relative_obstructions(self):
        target, delta = square(), R2-1
        eps, eta = delta/2, 3-2*R2
        additive = quartet(SPLITS[0], [Q2(F(1, 2))]*4, eps)
        relative = quartet(SPLITS[0], [Q2(F(1, 2))]*4, eta)
        self.assertEqual(maximum([abs(additive[p]-target[p]) for p in PAIRS]), eps)
        self.assertEqual(maximum([abs(relative[p]-target[p])/target[p] for p in PAIRS]), eta)
        self.assertEqual(2*R2-2*eps, 2+2*eps)
        self.assertEqual(2*R2*(1-eta), 2*(1+eta))
        self.assertGreater(eps.sign(), 0)
        self.assertGreater(eta.sign(), 0)
        OBS['sharp_square_approximation'] = dict(
            side_length=1, additive_error_minimum=eps.record(),
            pairwise_relative_error_minimum=eta.record(),
            both_attained_by_strict_positive_degree_three_tree=True,
            subdivision_rescaling_and_arbitrary_tree_size_do_not_remove_four_point_condition=True,
            excludes_near_isometric_riemannian_length_dimension_at_least_two=True,
            does_not_identify_Hausdorff_or_spectral_dimension=True)

    def test_03_square_mean_and_sharp_variance(self):
        ds, a, b = avg_square_pair()
        avg = mean(ds, [F(1, 2)]*2)
        self.assertEqual(avg, square())
        variances = {p: sum(((d[p]-avg[p])*(d[p]-avg[p])/2 for d in ds), Z) for p in PAIRS}
        delta = R2-1
        side_total = sum((variances[p] for p in SIDES), Z)
        self.assertEqual(side_total, 4*delta*delta)
        for p in SIDES:
            self.assertEqual(variances[p], delta*delta)
        for p in set(PAIRS)-set(SIDES):
            self.assertEqual(variances[p], Z)
        diff = [sums(d)[0]-sums(d)[1] for d in ds]
        self.assertEqual(sum(diff, Z), Z)
        self.assertEqual(sum((x*x/2 for x in diff), Z), 16*delta*delta)
        self.assertEqual(sum((abs(x)/2 for x in diff), Z), 4*delta)
        # An independent rational mixture also obeys the general mean inequality.
        other = [quartet(split, [Q2(F(i+1, 3)) for i in range(4)], Q2(2)) for split in SPLITS]
        ws = [F(1, 2), F(1, 3), F(1, 6)]
        av = sums(mean(other, ws))
        absolute = sum((w*abs(sums(d)[0]-sums(d)[1]) for w, d in zip(ws, other)), Z)
        self.assertGreaterEqual((absolute-(2*av[2]-av[0]-av[1])).sign(), 0)
        OBS['sharp_fluctuation_cost'] = dict(
            leaf_length=a.record(), central_length=b.record(),
            mean_is_exact_Euclidean_square=True,
            side_variance_sum_minimum=side_total.record(),
            side_rms_standard_deviation_minimum=delta.record(),
            opposite_side_sum_difference_rms_minimum=(4*delta).record(),
            attained_by_equal_mixture_of_two_trees=True,
            mean_metric_not_identified_with_communication_metric=True)

    def test_04_unit_edge_subdivision_approaches_mean_square(self):
        target = square(2+R2)
        cases = []
        for m in (1, 2, 7, 31, 64):
            k = math.isqrt(8*m*m)
            if 32*m*m > (2*k+1)**2:
                k += 1
            self.assertLessEqual((2*k-1)**2, 32*m*m)
            self.assertLessEqual(32*m*m, (2*k+1)**2)
            ds = []
            for split in SPLITS[:2]:
                d, adj = unit_subdivision(split, m, k)
                self.assertEqual(len(adj), 4*m+k+1)
                self.assertEqual(max(map(len, adj.values())), 3)
                ds.append({p: x/m for p, x in d.items()})
            av = mean(ds, [F(1, 2)]*2)
            error = maximum([abs(av[p]-target[p]) for p in PAIRS])
            self.assertGreaterEqual((Q2(F(1, 2*m))-error).sign(), 0)
            side_variance_sum = sum(((d[p]-av[p])*(d[p]-av[p])/2
                                     for d in ds for p in SIDES), Z)
            robust_rms_lower = (R2-1)*(2+R2)-3*error
            if robust_rms_lower.sign() > 0:
                self.assertGreaterEqual((side_variance_sum-4*robust_rms_lower*robust_rms_lower).sign(), 0)
            cases.append(dict(m=m, central_unit_edges=k, vertices=4*m+k+1,
                mean_error=error.record(), error_upper=str(F(1, 2*m)),
                robust_side_rms_lower=robust_rms_lower.record()))
        OBS['unit_edge_mean_limit'] = dict(cases=cases,
            scaling_factor='1/m', target_square_side='2+sqrt(2)',
            only_four_selected_leaf_distances_claimed=True,
            distribution_prepared_by_current_H_claimed=False,
            whole_network_three_dimensional_limit_claimed=False)

    def test_05_coherent_graph_and_reference_interface(self):
        ds, _, _ = avg_square_pair()
        metrics = {p: np.diag([float(d[p]) for d in ds]) for p in PAIRS}
        max_residual = 0.
        graph_coherence = 0.
        for phi in (0., .37, 1.11):
            psi = np.array([1/math.sqrt(2), 0, np.exp(1j*phi)/2, np.exp(1j*phi)/2])
            rho = np.outer(psi, psi.conj())
            graph_rho = psi.reshape(2, 2) @ psi.reshape(2, 2).conj().T
            graph_coherence = max(graph_coherence, abs(graph_rho[0, 1]))
            side_var = 0.
            for p, diag in metrics.items():
                op = np.kron(diag, np.eye(2))
                mu = float(np.trace(rho@op).real)
                var = float(np.trace(rho@op@op).real)-mu*mu
                max_residual = max(max_residual, abs(mu-float(square()[p])))
                if p in SIDES:
                    side_var += var
            max_residual = max(max_residual, abs(side_var-float(4*(R2-1)*(R2-1))))
        self.assertLess(max_residual, 2e-14)
        self.assertGreater(graph_coherence, .35)
        # Operator inequality on all three quartet sectors, not just two states.
        all_ds = [quartet(split, [ONE]*4, R2) for split in SPLITS]
        exact_slack = []
        for d in all_ds:
            x, y, z = sums(d)
            slack = (x+y+abs(x-y))/2-z
            self.assertGreaterEqual(slack.sign(), 0)
            exact_slack.append(slack.record())
        OBS['quantum_interface'] = dict(
            graph_reference_joint_dimension=4,
            nonzero_graph_offdiagonal=float(format(graph_coherence, '.13g')),
            maximum_float_residual=float(format(max_residual, '.13g')),
            exact_operator_slack_on_three_splits=exact_slack,
            arbitrary_joint_state_extension_from_commuting_diagonal_operators=True,
            graph_measurement_dephasing_or_postselection_required=False,
            actual_reader_or_preparation_device_built=False)

    def test_06_mean_distance_does_not_replace_round465_weight(self):
        # Unit-edge instance m=1,k=3: a side pair has distances 2 and 5.
        distances = [unit_subdivision(s, 1, 3)[0][0, 1].p for s in SPLITS[:2]]
        self.assertEqual(sorted(distances), [F(2), F(5)])
        weight_mean = sum((1/(d+1)**2 for d in distances), F(0))/2
        distance_mean = sum(distances, F(0))/2
        weight_at_mean = 1/(distance_mean+1)**2
        self.assertEqual(weight_mean, F(5, 72))
        self.assertEqual(weight_at_mean, F(4, 81))
        self.assertEqual(weight_mean-weight_at_mean, F(13, 648))
        OBS['existing_propagation_interface_boundary'] = dict(
            unit_edge_distances=[str(x) for x in distances],
            mean_distance=str(distance_mean),
            mean_inverse_square_weight=str(weight_mean),
            inverse_square_weight_of_mean=str(weight_at_mean),
            exact_gap=str(weight_mean-weight_at_mean),
            Jensen_inequality_is_mature_boundary_not_new_dynamics=True,
            mean_distance_alone_not_substituted_into_round465_bound=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=467, baseline_round=465, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(independent_of_round466=True, tree_four_point_theorem_is_mature=True,
            tree_path_metric_equals_physical_length_is_tested_extra_hypothesis=True,
            exact_and_near_isometric_length_limits_only=True,
            excludes_dimension_three_under_that_hypothesis=True,
            all_relational_or_averaged_geometries_excluded=False,
            Euclidean_square_mean_requires_nonvanishing_instantaneous_fluctuations=True,
            quantum_coherence_does_not_remove_commuting_distance_inequality=True,
            temporal_averaging_ergodicity_or_self_averaging_derived=False,
            mean_metric_generated_by_existing_H=False,
            spatial_dimension_derived=False, full_GR_goal_completed=False,
            new_physics_axiom_adopted=False, historical_files_changed=False,
            phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--dry-run', action='store_true')
    modes.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert result == json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as output:
            output.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
