"""Round 401: one fixed half-line history clock and persistent finite archives.

The Feynman clock and Bessel kernel are established tools. The audited interface
is an unconditional, reference-preserving bound for every fixed archived prefix.
An infinite internal workspace and a frozen-output contract are explicit inputs.
"""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
from strict_causal_limit_audit import bessel_real_quadrature

TARGET = Path(__file__).with_name('persistent_prefix_history_audit_results.json')
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)


def bessel_sequence(count, z, nodes=16384):
    theta = 2*np.pi*np.arange(nodes)/nodes
    # J_n(z) = mean exp(i*n*theta-i*z*sin(theta)).
    return np.fft.ifft(np.exp(-1j*z*np.sin(theta)))[:count]


def halfline_amplitudes(count, time):
    if time == 0:
        return np.eye(1, count, dtype=complex)[0]
    n = np.arange(count)
    bessel = bessel_sequence(count+2, 2*time)
    return (-1j)**n*(bessel[:count]+bessel[2:count+2])


def sine_integral(count, time, nodes=32768):
    theta = np.pi*(np.arange(nodes)+.5)/nodes
    phase = np.exp(-2j*time*np.cos(theta))*np.sin(theta)
    return np.array([2*np.mean(np.sin((n+1)*theta)*phase) for n in range(count)])


def prefix_bound(count, time):
    if count == 0:
        return 0.
    if time == 0:
        return 1.
    return min(1., count*(count+1)*(2*count+1)/(6*time*time))


def completion_time(count, epsilon):
    return math.sqrt(count*(count+1)*(2*count+1)/(6*epsilon))


def path_amplitudes(length, time):
    theta = np.pi*np.arange(1, length+1)/(length+1)
    basis = np.sqrt(2/(length+1))*np.sin(np.outer(np.arange(1, length+1), theta))
    return basis @ (np.exp(-2j*time*np.cos(theta))*basis[0])


def truncation_bound(length, time):
    """Norm bound after removing the scalar 2I: 2 sum_{k>=L}(2t)^k/k!."""
    if time == 0:
        return 0.
    x = 2*abs(time)
    if x >= length+1:
        return 2.
    log_value = (math.log(2)+length*math.log(x)-math.lgamma(length+1)
                 -math.log1p(-x/(length+1)))
    return min(2., math.exp(log_value))


def swap_gate(a, b, qubits):
    size = 2**qubits
    gate = np.zeros((size, size), complex)
    for j in range(size):
        bit_a = (j >> (qubits-1-a)) & 1
        bit_b = (j >> (qubits-1-b)) & 1
        k = j if bit_a == bit_b else j ^ (1 << (qubits-1-a)) ^ (1 << (qubits-1-b))
        gate[k, j] = 1.
    return gate


def on_site(gate, site, qubits):
    result = np.ones((1, 1), complex)
    for j in range(qubits):
        result = np.kron(result, gate if site == j else I)
    return result


def unitary(h, time):
    vals, vecs = np.linalg.eigh(h)
    return (vecs*np.exp(-1j*time*vals)) @ vecs.conj().T


def history_matrices(gates):
    d = gates[0].shape[0]
    length = len(gates)+1
    h = 2*np.eye(length*d, dtype=complex)
    prefixes = [np.eye(d, dtype=complex)]
    for j, gate in enumerate(gates):
        h[(j+1)*d:(j+2)*d, j*d:(j+1)*d] = gate
        h[j*d:(j+1)*d, (j+1)*d:(j+2)*d] = gate.conj().T
        prefixes.append(gate @ prefixes[-1])
    w = np.zeros_like(h)
    for j, v in enumerate(prefixes):
        w[j*d:(j+1)*d, j*d:(j+1)*d] = v
    return h, w, prefixes


def output_density(pure_data_reference, input_size, archive_size, ref_size):
    pure = pure_data_reference.reshape(input_size, archive_size*ref_size)
    return pure.T @ pure.conj()


def trace_distance(a, b):
    return float(np.abs(np.linalg.eigvalsh((a-b+(a-b).conj().T)/2)).sum()/2)


def initial_payload(qubits, reference_size, rng=None):
    size = 2**qubits
    if rng is None:
        assert reference_size == size
        data_reference = np.eye(size, dtype=complex)/np.sqrt(size)
    else:
        data_reference = rng.normal(size=(size, reference_size))+1j*rng.normal(size=(size, reference_size))
        data_reference /= np.linalg.norm(data_reference)
    initial = np.zeros((size, size, reference_size), complex)
    initial[:, 0, :] = data_reference
    return initial.reshape(size*size, reference_size)


@lru_cache(None)
def report():
    kernel = []
    for time in (.25, 2., 7., 19.):
        a = halfline_amplitudes(80, time)
        sine = sine_integral(80, time)
        bessel = bessel_sequence(82, 2*time)
        recurrence = (-1j)**np.arange(80)*(np.arange(80)+1)/time*bessel[1:81]
        kernel.append(dict(time=time, sine_error=float(np.max(abs(a-sine))),
                           recurrence_error=float(np.max(abs(a-recurrence))),
                           scalar_integral_error=float(abs(bessel[7]-bessel_real_quadrature(7, 2*time)))))
    norm_cases = []
    for time in (0., .7, 4., 13.):
        a = halfline_amplitudes(160, time)
        dt = 1e-5
        derivative = (halfline_amplitudes(160, time+dt)-halfline_amplitudes(160, time-dt))/(2*dt)
        adjacency = np.r_[a[1:], 0]+np.r_[0, a[:-1]]
        norm_cases.append(dict(time=time, norm=float(np.vdot(a, a).real),
                               schrodinger_error=float(np.linalg.norm(1j*derivative-adjacency))))
    prefix_cases = []
    for count in (1, 4, 8):
        for time in (5., 20., 100.):
            p = float(np.vdot(halfline_amplitudes(count, time), halfline_amplitudes(count, time)).real)
            prefix_cases.append(dict(prefix=count, time=time, unfinished_probability=p,
                                     all_later_bound=prefix_bound(count, time)))
    thresholds = [dict(prefix=n, epsilon=.01, time=completion_time(n, .01),
                       bound_at_threshold=prefix_bound(n, completion_time(n, .01)))
                  for n in (1, 4, 8, 16)]
    truncations = []
    for length, time in ((32, 2.), (48, 6.), (96, 15.)):
        exact = halfline_amplitudes(2*length, time)
        finite = np.r_[path_amplitudes(length, time), np.zeros(length)]
        truncations.append(dict(length=length, time=time,
                                sampled_norm_difference=float(np.linalg.norm(exact-finite)),
                                full_norm_bound=truncation_bound(length, time)))
    # Noncommuting complex gates, with B0 and B1 permanently idle after steps 1 and 3.
    gates = [swap_gate(0, 2, 4), on_site(unitary(.43*Y+.21*Z, 1.), 1, 4),
             swap_gate(1, 3, 4)]
    for j in range(6):
        gates.append(on_site(unitary(.3*X+.2*Y, .8+j*.1), j % 2, 4))
    h, w, prefixes = history_matrices(gates)
    length = len(prefixes)
    adjacency = np.diag(np.ones(length-1), 1)+np.diag(np.ones(length-1), -1)
    bare = np.kron(2*np.eye(length)+adjacency, np.eye(16))
    eigenvalues = np.linalg.eigvalsh(h)
    gauge = dict(conjugacy_error=float(np.linalg.norm(h-w@bare@w.conj().T, 2)),
                 unitary_error=float(np.linalg.norm(w.conj().T@w-np.eye(len(w)), 2)),
                 min_eigenvalue=float(eigenvalues[0]), max_eigenvalue=float(eigenvalues[-1]))
    rng = np.random.default_rng(401)
    direct_cases = []
    for reference_size in (3, 4):
        initial = initial_payload(2, reference_size, rng if reference_size == 3 else None)
        histories = [v@initial for v in prefixes]
        ideal = output_density(histories[3], 4, 4, reference_size)
        persistence = max(float(np.linalg.norm(output_density(s, 4, 4, reference_size)-ideal, 2))
                          for s in histories[3:])
        for time in (.6, 2.3):
            initial_clock = np.zeros((length*16, reference_size), complex)
            initial_clock[:16] = initial
            direct = (unitary(h, time) @ initial_clock).reshape(length, 16, reference_size)
            amplitudes = path_amplitudes(length, time)*np.exp(-2j*time)
            ansatz = np.array([a*s for a, s in zip(amplitudes, histories)])
            reduced = sum(output_density(s, 4, 4, reference_size) for s in direct)
            mixture = sum(abs(a)**2*output_density(s, 4, 4, reference_size)
                          for a, s in zip(amplitudes, histories))
            direct_cases.append(dict(reference_size=reference_size, time=time,
                                     history_error=float(np.linalg.norm(direct-ansatz)),
                                     partial_trace_mixture_error=float(np.linalg.norm(reduced-mixture, 2)),
                                     persistence_error=persistence,
                                     output_trace_distance=trace_distance(reduced, ideal),
                                     unfinished_probability=float(sum(abs(amplitudes[:3])**2))))
    # A fixed rule U_n=SWAP(A_n,B_n), not a new processor selected for each prefix.
    stream = []
    for m in (1, 2, 3):
        size = 2**m
        initial = initial_payload(m, size)
        states = [initial]
        for j in range(m):
            states.append(swap_gate(j, m+j, 2*m) @ states[-1])
        reductions = [output_density(s, size, size, size) for s in states]
        ideal = reductions[-1]
        for time in (3., 10., 30.):
            probabilities = abs(halfline_amplitudes(m, time))**2
            p = float(sum(probabilities))
            actual = (1-p)*ideal+sum(q*r for q, r in zip(probabilities, reductions[:-1]))
            stream.append(dict(prefix=m, time=time, bell_reference_trace_distance=trace_distance(actual, ideal),
                               unfinished_probability=p, all_later_bound=prefix_bound(m, time)))
    # Removing the frozen-output contract: first SWAP, then X on the same archive forever.
    erased = []
    for time in (5., 20., 80.):
        a = halfline_amplitudes(512, time)
        even = float(sum(abs(a[::2])**2))
        analytic = .5+float(bessel_sequence(2, 4*time)[1].real)/(4*time)
        erased.append(dict(time=time, wrong_archive_probability=even,
                           parity_formula=analytic, distance_from_half=abs(even-.5),
                           half_limit_envelope=1/(4*time), unfinished_probability=float(abs(a[0])**2)))
    return dict(round=401,
                scope='Fixed half-line Feynman clock, frozen finite output prefixes, declared infinite internal workspace; not a physical dimension derivation.',
                kernel_comparisons=kernel, norm_and_equation=norm_cases,
                prefix_bounds=prefix_cases, one_percent_thresholds=thresholds,
                finite_cutoff_checks=truncations, gauge_check=gauge,
                direct_unknown_reference_cases=direct_cases,
                fresh_archive_stream=stream, unfrozen_archive_counterexample=erased,
                same_fixed_hamiltonian_for_all_finite_prefixes=True,
                infinite_internal_workspace_assumed=True,
                frozen_output_contract_required=True,
                postselection_used=False, finite_time_exact_completion_claimed=False,
                arbitrary_repeated_readout_guaranteed=False,
                whole_state_convergence_claimed=False,
                full_sustained_cognition_model_constructed=False,
                spatial_dimension_generated=False)


class Audit(unittest.TestCase):
    def test_halfline_kernel_against_independent_integral(self):
        for row in report()['kernel_comparisons']:
            for key in ('sine_error', 'recurrence_error', 'scalar_integral_error'):
                self.assertLess(row[key], 3e-12)

    def test_initial_condition_norm_and_schrodinger_equation(self):
        for row in report()['norm_and_equation']:
            self.assertAlmostEqual(row['norm'], 1., places=11)
            self.assertLess(row['schrodinger_error'], 2e-9)

    def test_prefix_probability_and_all_later_envelope(self):
        for row in report()['prefix_bounds']:
            self.assertLessEqual(row['unfinished_probability'], row['all_later_bound']+1e-12)
        for row in report()['one_percent_thresholds']:
            self.assertAlmostEqual(row['bound_at_threshold'], row['epsilon'])

    def test_finite_cutoff_has_explicit_window_error(self):
        for row in report()['finite_cutoff_checks']:
            self.assertLessEqual(row['sampled_norm_difference'], row['full_norm_bound']+3e-13)
            self.assertLess(row['full_norm_bound'], 1e-4)

    def test_positive_fixed_history_hamiltonian(self):
        row = report()['gauge_check']
        self.assertLess(row['conjugacy_error'], 1e-12)
        self.assertLess(row['unitary_error'], 1e-12)
        self.assertGreaterEqual(row['min_eigenvalue'], -1e-12)
        self.assertLessEqual(row['max_eigenvalue'], 4.+1e-12)

    def test_full_hamiltonian_unknown_reference_unconditional_output(self):
        for row in report()['direct_unknown_reference_cases']:
            for key in ('history_error', 'partial_trace_mixture_error', 'persistence_error'):
                self.assertLess(row[key], 1e-12)
            self.assertLessEqual(row['output_trace_distance'], row['unfinished_probability']+1e-12)

    def test_correlated_finite_prefix_of_fixed_infinite_swap_stream(self):
        for row in report()['fresh_archive_stream']:
            self.assertLessEqual(row['bell_reference_trace_distance'], row['unfinished_probability']+1e-12)
            self.assertLessEqual(row['unfinished_probability'], row['all_later_bound']+1e-12)

    def test_unfinished_probability_does_not_replace_archive_contract(self):
        for row in report()['unfrozen_archive_counterexample']:
            self.assertAlmostEqual(row['wrong_archive_probability'], row['parity_formula'], places=11)
            self.assertLessEqual(row['distance_from_half'], row['half_limit_envelope']+1e-12)
            self.assertGreater(row['wrong_archive_probability'], .45)
            self.assertLess(row['unfinished_probability'], .001)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = dict(report())
    result['checks'] = dict(run=checks.testsRun, failures=len(checks.failures), errors=len(checks.errors))
    result['runtime'] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
