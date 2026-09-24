"""Round 437: exchange correlation capacity is not participation capacity.

All-N statements are proved in the note. Finite exact matrices, rational
certificates and independent spectral calculations audit their interfaces.
No graph samples, optimized threshold, or final Heisenberg compiler is used.
"""
import argparse
from fractions import Fraction
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import exchange_relation_audit as core

TARGET = Path(__file__).with_name('exchange_participation_capacity_audit_results.json')
OBS = {}


def norm(a):
    return float(np.linalg.norm(a, 2))


def full_model(nodes):
    edges = list(itertools.combinations(range(nodes), 2))
    m = len(edges)
    size = 2**(nodes+m)
    h = np.zeros((size, size), dtype=np.int64)
    for col in range(size):
        for e, (i, j) in enumerate(edges):
            mask = 1 << (m-1-e)
            h[col ^ mask, col] += 1
            if col & mask:
                h[col, col] += 1
                mi, mj = 1 << (nodes+m-1-i), 1 << (nodes+m-1-j)
                target = col ^ mi ^ mj if bool(col & mi) != bool(col & mj) else col
                h[target, col] += 1
    return h, edges


def symmetric_basis(nodes):
    v = np.zeros((2**nodes, nodes+1))
    for k in range(nodes+1):
        for label in range(2**nodes):
            if label.bit_count() == k:
                v[label, k] = 1/math.sqrt(math.comb(nodes, k))
    return v


def edge_hamiltonian(edges):
    size = 2**edges
    out = np.diag([2*k.bit_count() for k in range(size)]).astype(np.int64)
    for col in range(size):
        for e in range(edges):
            out[col ^ (1 << e), col] += 1
    return out


def activation_bound(time):
    assert 0 <= time <= math.pi
    amplitude = max(0., math.sin(time)-2*(1-math.cos(time)))
    return amplitude**2


class Audit(unittest.TestCase):
    def close(self, a, b, tolerance=2e-11):
        self.assertLess(float(np.linalg.norm(a-b)), tolerance)

    def test_01_known_star_capacity_and_sharp_saturation(self):
        rows = []
        for k in range(1, 7):
            swaps = [core.swap(k+1, 0, j) for j in range(1, k+1)]
            star = sum(swaps)
            values, vectors = np.linalg.eigh(star)
            ground = vectors[:, np.abs(values+1) < 1e-9]
            rho = ground@ground.conj().T/k
            self.assertEqual(ground.shape[1], k)
            self.assertAlmostEqual(float(values[0]), -1.)
            self.assertGreaterEqual(float(values[0]+1), -1e-11)
            singlets = [float(np.trace(rho@(np.eye(2**(k+1))-s)/2).real) for s in swaps]
            expected = (k+1)/(2*k)
            self.close(np.array(singlets), np.full(k, expected))
            self.assertAlmostEqual(sum(singlets), (k+1)/2)
            rows.append(dict(partners=k, minimum_star_swap_eigenvalue=float(values[0]),
                ground_dimension=k, equal_singlet_probabilities=singlets,
                exact_probability=str(Fraction(k+1, 2*k))))
        OBS['known_star_bound'] = dict(inequality='sum_j SWAP_0j >= -I',
            singlet_sum_bound='(k+1)/2', finite_spectra=rows,
            provenance='standard spin addition / singlet monogamy, not a new theorem',
            analytic_all_k_saturation=True)

    def test_02_strong_relation_threshold_and_weak_partner_count(self):
        rng = np.random.default_rng(437)
        n = 5
        swaps = [core.swap(n, 0, j) for j in range(1, n)]
        all_budgets = []
        for _ in range(16):
            v = rng.normal(size=(2**n, 3))+1j*rng.normal(size=(2**n, 3))
            rho = v@v.conj().T
            rho /= np.trace(rho)
            p = [float(np.trace(rho@(np.eye(2**n)-s)/2).real) for s in swaps]
            budget = sum(max(0., 2*x-1) for x in p)
            self.assertLessEqual(budget, 1+1e-11)
            all_budgets.append(budget)
        # Exact sharp family supplies arbitrarily many weakly entangled partners.
        examples = []
        for k in (1, 2, 3, 4, 6, 100):
            p = Fraction(k+1, 2*k)
            self.assertEqual(k*(2*p-1), 1)
            self.assertGreater(p, Fraction(1, 2))
            examples.append(dict(partners=k, singlet_probability=str(p),
                excess=str(p-Fraction(1, 2)), positive_excess_budget=1))
        delta = Fraction(1, 6)
        self.assertEqual(1/(2*delta), 3)
        OBS['correlation_capacity_scope'] = dict(
            positive_part_budget='sum_j max(0, 2 p_singlet(0,j)-1) <= 1',
            fixed_threshold_partner_bound='floor(1/(2 delta))',
            threshold_needed_for_three=str(delta), threshold_not_selected_by_exchange=True,
            arbitrary_weak_partner_family=examples, random_mixed_state_budgets=all_budgets,
            Hamiltonian_support_or_edge_register_occupancy_not_bounded_by_this_theorem=True)

    def test_03_environment_uniform_activation_lower_bound(self):
        h, edges = full_model(3)
        values, vectors = np.linalg.eigh(h)
        v0 = np.eye(64)[:, ::8]
        rows = []
        for t in (.1, .5, .8):
            evolved = (vectors*np.exp(-1j*t*values))@vectors.conj().T@v0
            minima = []
            for e in range(len(edges)):
                selected = [j for j in range(64) if j & (1 << (2-e))]
                block = evolved[selected, :]
                minimum = float(np.linalg.eigvalsh(block.conj().T@block)[0])
                self.assertGreaterEqual(minimum, activation_bound(t)-1e-12)
                minima.append(minimum)
            rows.append(dict(time=t, analytic_uniform_lower=activation_bound(t),
                full_unknown_data_operator_minima=minima))
        # Check the lemma independently on environments with very large drift.
        rng = np.random.default_rng(43703)
        environments = []
        for d, scale in ((3, 1.), (5, 100.)):
            w = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
            k = scale*(w+w.conj().T)/2
            q, _ = np.linalg.qr(w)
            a = (q*np.linspace(0., 2., d))@q.conj().T
            hs = np.block([[k, np.eye(d)], [np.eye(d), k+a]])
            u = core.evolve(hs, .5)
            lower = float(np.linalg.eigvalsh(u[d:, :d].conj().T@u[d:, :d])[0])
            self.assertGreater(lower, activation_bound(.5)-1e-12)
            environments.append(dict(dimension=d, drift_norm=norm(k), occupation_minimum=lower))
        rational_amplitude = Fraction(1, 2)-Fraction(1, 2)**3/6-Fraction(1, 2)**2
        self.assertEqual(rational_amplitude, Fraction(11, 48))
        rational_probability = rational_amplitude**2
        self.assertGreater(rational_probability, Fraction(1, 20))
        self.assertGreater(activation_bound(.5), float(rational_probability))
        OBS['all_N_all_data_activation'] = dict(
            exact_certificate_at_half=str(rational_probability),
            formula='[max(0,sin(t)-2(1-cos(t)))]^2, 0 <= t <= pi',
            arbitrary_reference_included=True, independent_of_environment_norm_and_N=True,
            expected_degree_lower_at_half='(N-1)*121/2304',
            explicit_operator_tests=rows, unrelated_large_environment_tests=environments)

    def test_04_complete_symmetric_sector_factorization(self):
        records = []
        for n in (2, 3, 4):
            h, edges = full_model(n)
            m = len(edges)
            s = symmetric_basis(n)
            embed = np.kron(s, np.eye(2**m))
            effective = np.kron(np.eye(n+1), edge_hamiltonian(m))
            error = float(np.linalg.norm(h@embed-embed@effective))
            self.close(h@embed, embed@effective)
            records.append(dict(nodes=n, data_dimension=2**n, relation_qubits=m,
                whole_dimension=len(h), checked_isometry_columns=embed.shape[1], residual=error))
        h, _ = full_model(3)
        s = symmetric_basis(3)
        z = np.zeros((8, 1)); z[0, 0] = 1
        initial_map = np.kron(s, z)
        time = .5
        one_edge = core.evolve(np.array([[0., 1.], [1., 2.]]), time)[:, :1]
        edge_state = np.kron(np.kron(one_edge, one_edge), one_edge)
        actual = core.evolve(h, time)@initial_map
        ideal = np.kron(s, edge_state)
        self.close(actual, ideal)
        coefficients = np.zeros((4, 2), complex)
        coefficients[0, 0] = 1/math.sqrt(2)
        coefficients[3, 1] = 1j/math.sqrt(2)
        self.close(actual@coefficients, ideal@coefficients)
        self.close((actual@coefficients).conj().T@(actual@coefficients), coefficients.conj().T@coefficients)
        OBS['symmetric_sector'] = dict(generator_checks=records,
            full_complex_unknown_symmetric_data_and_reference_retained=True,
            exact_intertwining_not_only_chosen_known_input=True,
            single_edge_H=[[0, 1], [1, 2]],
            edge_probability='sin(sqrt(2)*t)^2/2',
            edge_probability_at_half=math.sin(math.sqrt(2)/2)**2/2,
            complete_generator_identity_verified=True)

    def test_05_joint_graph_statistics_and_exact_sparse_tail(self):
        # N=3 direct full Hamiltonian measurement versus exact product law.
        h, _ = full_model(3)
        initial = np.zeros(64); initial[0] = 1
        t = math.pi/(2*math.sqrt(2))
        final = core.evolve(h, t)@initial
        edge_probabilities = np.sum(abs(final.reshape(8, 8))**2, axis=0)
        self.close(edge_probabilities, np.full(8, 1/8))
        # For any graph max degree <= 3 implies at most floor(3N/2) edges.
        # At t=pi/(2 sqrt(2)), every graph has probability 2^-M exactly.
        n = 32
        m = math.comb(n, 2)
        edge_cap = 3*n//2
        tail = Fraction(sum(math.comb(m, k) for k in range(edge_cap+1)), 2**m)
        self.assertLess(tail, Fraction(1, 10**70))
        self.assertEqual(m, 496)
        OBS['dense_graph_witness'] = dict(
            measurement_time='pi/(2 sqrt(2))', probability_per_possible_edge='1/2',
            whole_graph_law='G(N,1/2), a measurement distribution, not an assumed classical graph',
            direct_three_node_joint_probabilities=edge_probabilities.tolist(),
            witness_nodes=n, possible_edges=m, sparse_edge_cap=edge_cap,
            probability_maximum_degree_at_most_three_upper_fraction=str(tail),
            probability_maximum_degree_at_most_three_upper=float(tail),
            exact_expected_degree=str(Fraction(n-1, 2)),
            fixed_degree_sparse_probability_tends_to_zero=True,
            no_asymptotic_monte_carlo_or_graph_sampling_used=True)

    def test_06_internal_resource_account_and_simulator_robustness(self):
        rows = []
        for n in (2, 3, 4):
            h, edges = full_model(n)
            m = len(edges)
            cols = np.arange(0, len(h), 2**m)
            self.assertTrue(np.array_equal(h[np.ix_(cols, cols)], np.zeros((2**n, 2**n))))
            self.assertTrue(np.array_equal(h[:, cols].T@h[:, cols], m*np.eye(2**n)))
            rows.append(dict(nodes=n, relation_qubits=m, blank_mean_energy=0,
                blank_energy_variance=m, total_qubits=n+m))
        h1 = np.array([[0., 1.], [1., 2.]])
        t = math.pi/(2*math.sqrt(2))
        psi = core.evolve(h1, t)[:, 0]
        x = np.array([[0., 1.], [1., 0.]])
        occupation = float(abs(psi[1])**2)
        flip_energy = float(np.vdot(psi, x@psi).real)
        self.assertAlmostEqual(occupation, .5)
        self.assertAlmostEqual(flip_energy, -1.)
        self.assertAlmostEqual(flip_energy+2*occupation, 0.)
        n = 32; m = math.comb(n, 2)
        tail = Fraction(sum(math.comb(m, k) for k in range(3*n//2+1)), 2**m)
        # Reuse 436 trace-distance contract; this is not a new hardware run.
        eta = epsilon = Fraction(1, 10000)
        self.assertLess(t, 2.)
        trace_distance_upper = 2*eta+2*epsilon
        dense_lower = 1-tail-trace_distance_upper
        self.assertGreater(dense_lower, Fraction(9995, 10000))
        OBS['resource_and_simulation_scope'] = dict(
            blank_state_energy_identities=rows, all_possible_edge_carriers_already_allocated=True,
            no_constant_total_hardware_per_actor_claimed=True,
            occupation_per_edge_at_witness_time=occupation,
            off_diagonal_energy_per_edge_at_witness_time=flip_energy,
            conserved_total_mean_energy_per_edge=flip_energy+2*occupation,
            initial_energy_above_sector_ground='(sqrt(2)-1)*M',
            initial_energy_above_full_ground_at_least_sector_value=True,
            final_static_simulator_trace_distance_upper=str(trace_distance_upper),
            simulated_nonsparse_event_probability_lower=float(dense_lower),
            cited_436_contract_applied_not_final_Heisenberg_hardware_run=True,
            only_finite_N_finite_time_positive_error_simulation_claimed=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=437, baseline_round=436, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(known_singlet_monogamy_reused_not_claimed_as_new=True,
            exchange_correlation_budget_distinguished_from_participation_capacity=True,
            all_N_all_unknown_data_activation_lower_bound_proved=True,
            exact_joint_dense_graph_counterexample_with_reference_proved=True,
            additional_relation_hardware_and_initial_energy_accounted=True,
            naive_433_automatic_sparse_limit_refuted=True,
            all_exchange_based_sparse_mechanisms_refuted=False,
            degree_three_or_physical_space_generated=False,
            final_Heisenberg_coupling_list_or_hardware_run_completed=False,
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
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
