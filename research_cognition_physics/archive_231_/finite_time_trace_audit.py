"""Round 390: finite-time coalescence under integrable linear generators.

This is a conditional dynamics audit, not a definition of actual history.
NumPy only. Default read-only; --write-results creates its output exclusively.
"""
import argparse
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np


HERE = Path(__file__).resolve().parent
TARGET = HERE / 'finite_time_trace_audit_results.json'
TOL = 3e-11
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)
P0 = np.diag([1., 0.]).astype(complex)
P1 = I - P0
LOWER = np.array([[0, 1], [0, 0]], dtype=complex)
BELL = np.outer(np.array([1, 0, 0, 1]), np.array([1, 0, 0, 1])) / 2


def tn(a):
    return float(np.abs(np.linalg.eigvalsh((a + a.conj().T) / 2)).sum())


def distance(a, b):
    return tn(a - b) / 2


def random_state(rng, d):
    a = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
    a = a @ a.conj().T
    return a / np.trace(a)


def rotation(pauli, angle):
    return math.cos(angle/2)*I - 1j*math.sin(angle/2)*pauli


SEGMENTS = ((.7, .4, I),
            (1.1, .3, rotation(Y, .83)),
            (.8, .5, rotation(Z, .57) @ rotation(Y, 1.19)))


def kraus_channel(operators, a):
    return sum((k @ a @ k.conj().T for k in operators), np.zeros_like(a))


def amplitude_kraus(q, r):
    k0 = np.diag([1., math.sqrt(q)]).astype(complex)
    k1 = math.sqrt(1-q)*LOWER
    return [r @ k @ r.conj().T for k in (k0, k1)]


def generator(a, h, jumps):
    value = -1j*(h @ a - a @ h)
    for v in jumps:
        w = v.conj().T @ v
        value += v @ a @ v.conj().T - (w @ a + a @ w)/2
    return value


def supermatrix(channel, d=2):
    basis = np.eye(d*d, dtype=complex)
    return np.column_stack([channel(v.reshape((d, d), order='F')).reshape(-1, order='F')
                            for v in basis])


def choi(channel, d=2):
    result = np.zeros((d*d, d*d), dtype=complex)
    for a in range(d):
        for b in range(d):
            e = np.zeros((d, d), dtype=complex)
            e[a, b] = 1
            result += np.kron(e, channel(e))
    return result


def matrices():
    answer = []
    for gamma, duration, r in SEGMENTS:
        k = amplitude_kraus(math.exp(-gamma*duration), r)
        v = math.sqrt(gamma)*r @ LOWER @ r.conj().T
        answer.append((supermatrix(lambda a: kraus_channel(k, a)),
                       supermatrix(lambda a: generator(a, 0*I, [v]))))
    return answer


def composed(a, reference=1):
    for gamma, duration, r in SEGMENTS:
        operators = [np.kron(k, np.eye(reference))
                     for k in amplitude_kraus(math.exp(-gamma*duration), r)]
        a = kraus_channel(operators, a)
    return a


def rk4(g, value, duration, steps=96):
    dt = duration/steps
    for _ in range(steps):
        k1 = g @ value
        k2 = g @ (value + dt*k1/2)
        k3 = g @ (value + dt*k2/2)
        k4 = g @ (value + dt*k3)
        value = value + dt*(k1 + 2*k2 + 2*k3 + k4)/6
    return value


def depolarize(a, lam):
    d = len(a)
    return lam*a + (1-lam)*np.trace(a)*np.eye(d)/d


def system_reset(a, reference):
    reduced = np.trace(a.reshape(2, reference, 2, reference), axis1=0, axis2=2)
    return np.kron(I/2, reduced)


def partial_trace(a, keep):
    axes = (1, 3) if keep == 0 else (0, 2)
    return np.trace(a.reshape(2, 2, 2, 2), axis1=axes[0], axis2=axes[1])


def experiment():
    pairs = matrices()
    f = pairs[2][0] @ pairs[1][0] @ pairs[0][0]
    ode = np.eye(4, dtype=complex)
    for (_, duration, _), (_, g) in zip(SEGMENTS, pairs):
        ode = rk4(g, ode, duration)
    budget = 2*sum(gamma*duration for gamma, duration, _ in SEGMENTS)
    rng = np.random.default_rng(390)
    ratios = []
    for reference in (1, 2, 3):
        local = []
        for _ in range(32):
            a, b = [random_state(rng, 2*reference) for _ in range(2)]
            local.append(distance(composed(a, reference), composed(b, reference))/distance(a, b))
        ratios.append(dict(reference_dimension=reference, pairs=len(local),
                           min_ratio=min(local), max_ratio=max(local)))
    lam = .2
    inverse_p0 = (P0-(1-lam)*I/2)/lam
    return dict(
        round=390, scientific_base_through_round=389,
        hypothesis='Dissipation alone allows distinct complete states to merge exactly in finite time.',
        scope=dict(
            ontology='Only cognition; reference and environment here are parts inside the declared whole.',
            additional_inputs=['Fixed finite-dimensional complete state space',
                               'Same state-independent linear time-local generator for both preparations',
                               'Measurable generator with finite integrated norm including the endpoint',
                               'Chosen continuous parameter and physical CPTP propagation'],
            conclusion='These inputs forbid coalescence of distinct states even without unitarity.',
            actual_event_trace_definition_derived=False,
            same_trajectory_recurrence_excluded=False,
            recurrence_assumed_allowed_by_cognitive_principles=False,
            overwrite_without_rollback_is_new_candidate=True,
            forward_recurrence_accepted_in_user_clarification=True,
            no_recurrence_axiom_imposed=False,
            global_dynamics_derived_from_ontology=False,
            three_dimensional_space_derived=False,
            budget_is_thermodynamic_or_energy_cost=False),
        time_ordered_case=dict(
            jump_rates=[g for g, _, _ in SEGMENTS], durations=[dt for _, dt, _ in SEGMENTS],
            reference_uniform_majorant=budget, certified_ratio_lower_bound=math.exp(-budget),
            generator_commutator_frobenius=float(np.linalg.norm(pairs[0][1]@pairs[1][1]-pairs[1][1]@pairs[0][1])),
            wrong_order_frobenius=float(np.linalg.norm(f-pairs[0][0]@pairs[1][0]@pairs[2][0])),
            rk4_vs_exact_max_error=float(np.max(np.abs(ode-f))),
            superoperator_smallest_singular_value=float(np.linalg.svd(f, compute_uv=False)[-1]),
            sampled_contractions=ratios),
        sharp_qubit_example=dict(rate=1., target_pair_ratio=1e-6,
                                 required_duration=math.log(1e6),
                                 inverse_at_lambda_point_two_eigenvalues=np.linalg.eigvalsh(inverse_p0).tolist(),
                                 unamplified_generator_norm=1.,
                                 bell_input_amplified_norm=tn(system_reset(BELL, 2)-BELL),
                                 arbitrary_reference_safe_majorant=1.5),
        finite_endpoint_example=dict(T=1., lambda_rule='max(1-t/T,0)^2',
                                     differentiability='C1 at T', channel_rank_at_T=1,
                                     generator_rate_before_T='2/(T-t)',
                                     cumulative_rate_at_endpoint='infinite'),
        bounded_internal_swap=dict(g=1., global_hamiltonian_norm=1.,
                                   reset_time=math.pi/2,
                                   reduced_lambda='cos(g*t)^2',
                                   reduced_rate_before_T='2*g*tan(g*t)',
                                   complete_state_distance_at_T=1., local_state_distance_at_T=0.,
                                   information_transferred_to_internal_partner=True),
        sources=['https://arxiv.org/pdf/0801.4100', 'https://arxiv.org/pdf/1710.06771'])


class Checks(unittest.TestCase):
    def test_01_exact_piecewise_channels_are_cptp(self):
        for gamma, duration, r in SEGMENTS:
            k = amplitude_kraus(math.exp(-gamma*duration), r)
            np.testing.assert_allclose(sum(v.conj().T@v for v in k), I, atol=TOL)
            j = choi(lambda a: kraus_channel(k, a))
            self.assertGreaterEqual(np.linalg.eigvalsh(j)[0], -TOL)
        self.assertGreaterEqual(np.linalg.eigvalsh(choi(composed))[0], -TOL)

    def test_02_noncommuting_time_order_matches_independent_ode(self):
        m = experiment()['time_ordered_case']
        self.assertGreater(m['generator_commutator_frobenius'], .05)
        self.assertGreater(m['wrong_order_frobenius'], .01)
        self.assertLess(m['rk4_vs_exact_max_error'], 3e-11)
        self.assertGreater(m['superoperator_smallest_singular_value'], .1)

    def test_03_bound_holds_with_entangled_internal_references(self):
        m = experiment()['time_ordered_case']
        for r in m['sampled_contractions']:
            self.assertGreaterEqual(r['min_ratio'], m['certified_ratio_lower_bound']-TOL)
            self.assertLessEqual(r['max_ratio'], 1+TOL)
        self.assertAlmostEqual(m['reference_uniform_majorant'], 2.02)

    def test_04_linear_inverse_is_not_a_physical_recovery(self):
        lam = .2
        inv = (P0-(1-lam)*I/2)/lam
        np.testing.assert_allclose(depolarize(inv, lam), P0, atol=TOL)
        np.testing.assert_allclose(np.linalg.eigvalsh(inv), [-2, 3], atol=TOL)
        f = supermatrix(lambda a: depolarize(a, lam))
        self.assertAlmostEqual(float(np.linalg.det(f).real), lam**3)

    def test_05_sharp_pair_contraction_and_finite_accuracy_time(self):
        rng = np.random.default_rng(391)
        for t in (.1, 1., 4., math.log(1e6)):
            lam = math.exp(-t)
            for _ in range(6):
                a, b = [random_state(rng, 2) for _ in range(2)]
                self.assertAlmostEqual(distance(depolarize(a, lam), depolarize(b, lam)), lam*distance(a, b))
        t = math.log(1e6)
        self.assertAlmostEqual(math.exp(-t), 1e-6)
        self.assertGreater(math.exp(-(t-.01)), 1e-6)

    def test_06_unamplified_norm_cannot_be_used_for_arbitrary_reference(self):
        rng = np.random.default_rng(392)
        for _ in range(64):
            a = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
            a += a.conj().T
            self.assertLessEqual(tn(depolarize(a, 0)-a), tn(a)+TOL)
        self.assertAlmostEqual(tn(depolarize(P0, 0)-P0), 1.)
        self.assertAlmostEqual(tn(system_reset(BELL, 2)-BELL), 1.5)
        jumps = [np.kron(p/2, I) for p in (X, Y, Z)]
        np.testing.assert_allclose(generator(BELL, np.zeros((4, 4)), jumps),
                                   system_reset(BELL, 2)-BELL, atol=TOL)

    def test_07_c1_channel_can_reset_at_finite_singular_endpoint(self):
        for t in (0., .2, .8, .99, 1., 1.2):
            lam = max(1-t, 0)**2
            f = supermatrix(lambda a: depolarize(a, lam))
            self.assertGreaterEqual(np.linalg.eigvalsh(choi(lambda a: depolarize(a, lam)))[0], -TOL)
            if t < 1:
                self.assertAlmostEqual(lam, math.exp(2*math.log(1-t)))
                self.assertAlmostEqual(-2*(1-t), -(2/(1-t))*lam)
            else:
                self.assertEqual(np.linalg.matrix_rank(f), 1)
                np.testing.assert_allclose(depolarize(P0, lam), depolarize(P1, lam), atol=TOL)
        # The left derivative at T tends to the zero right derivative.
        for h in (1e-2, 1e-3, 1e-4):
            quotient = (depolarize(P0, h*h)-depolarize(P0, 0))/h
            np.testing.assert_allclose(quotient, h*(P0-I/2), atol=2e-12)

    def test_08_singular_reduced_rate_does_not_require_unbounded_total_h(self):
        swap = np.array([[1, 0, 0, 0], [0, 0, 1, 0],
                         [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex)
        self.assertAlmostEqual(np.linalg.norm(swap, 2), 1.)
        for t in (.1, .6, math.pi/2-1e-5, math.pi/2):
            u = math.cos(t)*np.eye(4) - 1j*math.sin(t)*swap
            outputs = [u @ np.kron(a, I/2) @ u.conj().T for a in (P0, P1)]
            for a, whole in zip((P0, P1), outputs):
                np.testing.assert_allclose(partial_trace(whole, 0), depolarize(a, math.cos(t)**2), atol=TOL)
            self.assertAlmostEqual(distance(*outputs), 1.)
            self.assertAlmostEqual(distance(*(partial_trace(a, 0) for a in outputs)), math.cos(t)**2)
        for a, whole in zip((P0, P1), outputs):
            np.testing.assert_allclose(partial_trace(whole, 1), a, atol=TOL)
        self.assertGreater(2*math.tan(math.pi/2-1e-5), 1e5)

    def test_09_gksl_majorant_with_hamiltonian_and_arbitrary_hermitian_input(self):
        h = .3*X + .2*Z
        jumps = [math.sqrt(.4)*LOWER, .2*Y]
        b = 2*np.linalg.norm(h, 2) + 2*sum(np.linalg.norm(v, 2)**2 for v in jumps)
        rng = np.random.default_rng(393)
        for r in (1, 2, 3):
            for _ in range(20):
                a = rng.normal(size=(2*r, 2*r)) + 1j*rng.normal(size=(2*r, 2*r))
                a += a.conj().T
                out = generator(a, np.kron(h, np.eye(r)), [np.kron(v, np.eye(r)) for v in jumps])
                self.assertLessEqual(tn(out), b*tn(a)+TOL)


def report():
    return experiment()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    run = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Checks))
    if not run.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result['checks'] = dict(run=run.testsRun, failures=len(run.failures), errors=len(run.errors))
    result['runtime'] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
