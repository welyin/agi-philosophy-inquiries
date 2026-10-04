"""Round 426: closed finite-budget control groups and a strict closure boundary.

General topology and infinite-block statements are proved in the note.
Finite tests certify formulas and examples, not infinite controllability.
"""
import argparse
from fractions import Fraction as F
import io
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'closed_control_lie_audit_results.json'
OBS = {}
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)


def norm(a):
    return float(np.linalg.norm(a, 2))


def exp_h(h, time=1.):
    values, vectors = np.linalg.eigh(h)
    return (vectors * np.exp(-1j * time * values)) @ vectors.conj().T


def bracket(a, b):
    return 1j * (a @ b - b @ a)


def fraction_rank(rows):
    a = [[F(x) for x in row] for row in rows]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((k for k in range(rank, len(a)) if a[k][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [x / scale for x in a[rank]]
        for k in range(rank + 1, len(a)):
            scale = a[k][col]
            if scale:
                a[k] = [x - scale * y for x, y in zip(a[k], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def integer_vector(a):
    values = np.concatenate([a.real.ravel(), a.imag.ravel()])
    assert np.max(np.abs(values - np.rint(values))) < 1e-10
    return [int(x) for x in np.rint(values)]


def integer_lie_rank(generators):
    basis = []
    rows = []
    queue = [g.copy() for g in generators]
    while queue:
        candidate = queue.pop(0)
        row = integer_vector(candidate)
        divisor = math.gcd(*row)
        if not divisor:
            continue
        candidate /= divisor
        row = [x // divisor for x in row]
        if fraction_rank(rows + [row]) == len(rows):
            continue
        basis.append(candidate)
        rows.append(row)
        queue.extend(bracket(g, candidate) for g in generators)
    return len(basis)


def odd_sqrt(x):
    x = np.asarray(x)
    return np.sign(x) * np.sqrt(np.abs(x))


def bernstein_odd_sqrt(x, degree):
    """Bernstein polynomial after mapping [-1,1] to [0,1]."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    nodes = np.arange(degree + 1)
    coeff = odd_sqrt(2 * nodes / degree - 1)
    log_binom = np.array([math.lgamma(degree + 1) - math.lgamma(k + 1) -
                         math.lgamma(degree - k + 1) for k in nodes])
    for index, value in np.ndenumerate(x):
        y = (value + 1) / 2
        if y <= 0:
            out[index] = -1.
        elif y >= 1:
            out[index] = 1.
        else:
            weights = np.exp(log_binom + nodes * np.log(y) + (degree - nodes) * np.log1p(-y))
            out[index] = weights @ coeff
    return out


def protocol_endpoint(weight, protocol):
    out = I.copy()
    for duration, ux, uz in protocol:
        assert abs(ux) + abs(uz) <= 1 + 1e-14
        out = exp_h(weight * (ux * X + uz * Z), duration) @ out
    return out


class Audit(unittest.TestCase):
    def test_01_banach_doubling_and_escape(self):
        rng = np.random.default_rng(42601)
        minimum_slack = float('inf')
        delta = 0.5
        for _ in range(100):
            a = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
            a *= rng.uniform(0.01, delta) / norm(a)
            t = np.eye(3) + a
            slack = norm(t @ t - np.eye(3)) - (2 - delta) * norm(a)
            minimum_slack = min(minimum_slack, slack)
        self.assertGreaterEqual(minimum_slack, -1e-13)
        unitary = exp_h(Z, 0.01)
        base = norm(unitary - I)
        predicted = math.floor(math.log(delta / base) / math.log(2 - delta)) + 1
        first_exit = next(k for k in range(predicted + 1) if norm(np.linalg.matrix_power(unitary, 2 ** k) - I) >= delta)
        self.assertLessEqual(first_exit, predicted)
        OBS['banach_doubling_minimum_sampled_slack'] = minimum_slack
        OBS['unitary_first_exit_doublings'] = first_exit
        OBS['general_bound_doublings'] = predicted

    def test_02_closed_control_groups_do_not_select_three(self):
        qutrit_x = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], complex)
        qutrit_z = np.diag([-4., -1., 5.]).astype(complex)
        ranks = [integer_lie_rank([Z]), integer_lie_rank([X, Z]),
                 integer_lie_rank([qutrit_x, qutrit_z])]
        self.assertEqual(ranks, [1, 3, 8])
        OBS['exact_rational_lie_dimensions'] = ranks
        OBS['dimension_examples_are_operational_groups_not_physical_spaces'] = True

    def test_03_nested_brackets_and_exact_vandermonde(self):
        error = 0.
        for j in range(1, 7):
            w = 4. ** (-j)
            current = w * Z
            for k in range(5):
                error = max(error, norm(current - w ** (2 * k + 1) * Z))
                current = -0.25 * bracket(w * X, bracket(w * X, current))
        self.assertLess(error, 1e-16)
        certificates = []
        for size in range(1, 5):
            weights = [F(1, 4 ** j) for j in range(1, size + 1)]
            matrix = [[w ** (2 * k + 1) for k in range(size)] for w in weights]
            determinant = math.prod(weights) * math.prod(weights[j] ** 2 - weights[i] ** 2
                                                         for i in range(size) for j in range(i + 1, size))
            self.assertNotEqual(determinant, 0)
            self.assertEqual(fraction_rank(matrix), size)
            certificates.append(dict(blocks=size, determinant=str(determinant), rank=size))
        OBS['nested_bracket_profile_error'] = error
        OBS['finite_block_exact_vandermonde_certificates'] = certificates

    def test_04_uniform_odd_polynomial_approximation(self):
        points = np.linspace(-1, 1, 1201)
        rows = []
        previous = float('inf')
        for degree in (15, 63, 255):
            values = bernstein_odd_sqrt(points, degree)
            sampled = float(np.max(np.abs(values - odd_sqrt(points))))
            odd_error = float(np.max(np.abs(values + values[::-1])))
            proven = np.sqrt(2) * degree ** (-0.25)
            self.assertLess(sampled, proven)
            self.assertLess(sampled, previous)
            self.assertLess(odd_error, 3e-13)
            rows.append(dict(degree=degree, sampled_error=sampled,
                             analytic_uniform_error_bound=float(proven), oddness_error=odd_error))
            previous = sampled
        OBS['bernstein_approximation_certificates'] = rows

    def test_05_target_tail_and_infinite_cost_certificates(self):
        certificates = []
        for j in (1, 4, 10, 20):
            theta, w = F(1, 2 ** j), F(1, 4 ** j)
            lower = (theta - theta ** 3 / 6) / w
            self.assertEqual(lower, 2 ** j - F(1, 6 * 2 ** j))
            target = exp_h(Z, float(theta))
            distance = norm(target - I)
            self.assertAlmostEqual(distance, 2 * np.sin(float(theta) / 2), places=14)
            self.assertLessEqual(distance, float(theta) + 1e-15)
            certificates.append(dict(block=j, exact_required_cost_lower_bound=str(lower),
                                     unitary_distance_to_identity=distance))
        OBS['weak_block_exact_cost_certificates'] = certificates
        OBS['target_tail_operator_norm_upper_bound_after_block_N'] = '2^(-(N+1))'

    def test_06_all_protocol_speed_bound_witnesses(self):
        rng = np.random.default_rng(42606)
        protocol = []
        for _ in range(24):
            controls = rng.normal(size=2)
            controls /= np.sum(np.abs(controls))
            protocol.append((float(rng.uniform(0.01, 0.2)), *controls))
        budget = sum(t * (abs(x) + abs(z)) for t, x, z in protocol)
        maximum_speed_ratio = 0.
        maximum_bound_violation = 0.
        plus = np.array([1, 1], complex) / np.sqrt(2)
        for j in range(1, 13):
            w, theta = 4. ** (-j), 2. ** (-j)
            u = protocol_endpoint(w, protocol)
            target = exp_h(Z, theta)
            maximum_speed_ratio = max(maximum_speed_ratio, norm(u - I) / (w * budget))
            # Pure output-state trace distance is a legitimate channel witness.
            a, b = u @ plus, target @ plus
            witness = norm(np.outer(a, a.conj()) - np.outer(b, b.conj()))
            analytic_lower = max(0., np.sin(theta) - budget * w)
            maximum_bound_violation = max(maximum_bound_violation, analytic_lower - witness)
        self.assertLessEqual(maximum_speed_ratio, 1 + 1e-12)
        self.assertLessEqual(maximum_bound_violation, 1e-12)
        OBS['sample_protocol_budget'] = float(budget)
        OBS['sample_protocol_maximum_speed_ratio'] = float(maximum_speed_ratio)
        OBS['channel_witness_bound_maximum_violation'] = float(maximum_bound_violation)

    def test_07_uniform_finite_error_resource_certificate(self):
        certificates = []
        for budget in (1, 2, 3, 8, 31, 64):
            j = 1
            while F(1, 2 ** j) > F(1, 2 * budget):
                j += 1
            theta = F(1, 2 ** j)
            lower = theta - theta ** 3 / 6 - budget * theta ** 2
            universal = F(1, 12 * budget)
            self.assertGreater(lower, universal)
            certificates.append(dict(budget=budget, chosen_block=j,
                                     exact_error_lower_bound=str(lower),
                                     conservative_error_bound=str(universal)))
        OBS['finite_error_resource_certificates'] = certificates

    def test_08_commutator_approximation_counts_resource(self):
        rows = []
        for repetitions in (4, 16, 64, 256):
            step = 1 / np.sqrt(repetitions)
            error = 0.
            for w in np.linspace(0, 0.25, 31):
                a, b = exp_h(w * X, step), exp_h(w * Z, step)
                cycle = a @ b @ a.conj().T @ b.conj().T
                actual = np.linalg.matrix_power(cycle, repetitions)
                # [-iwX,-iwZ]=2iw^2 Y.
                target = exp_h(-2 * w * w * Y)
                error = max(error, norm(actual - target))
            rows.append(dict(repetitions=repetitions, sampled_uniform_error=error,
                             protocol_action=4 * repetitions * step))
        self.assertTrue(all(rows[k + 1]['sampled_uniform_error'] < rows[k]['sampled_uniform_error'] for k in range(3)))
        self.assertGreater(rows[-1]['protocol_action'], rows[0]['protocol_action'])
        OBS['commutator_approximation_with_action'] = rows


def run():
    output = io.StringIO()
    tests = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=426, baseline_round=425, status='verified_conditional_control_to_Lie_bridge',
                tests_run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors),
                python=platform.python_version(), numpy=np.__version__, observations=OBS,
                scope=dict(frozen_round415_compactness_reused=True,
                           baire_exact_reachable_control_group_is_Lie=True,
                           closure_assumption_derived_from_cognition=False,
                           half_operation_needed_for_this_bridge=False,
                           exact_budget_continuity_given_as_input=False,
                           complete_operation_topology_still_input=True,
                           infinitely_weak_block_target_in_norm_closure=True,
                           weak_block_target_has_finite_exact_budget=False,
                           approximate_implementation_excluded=False,
                           finite_error_lower_bound_claimed_optimal=False,
                           actual_position_interface_still_input=True,
                           operational_group_dimension_identified_as_space=False,
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
            stream.write(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
