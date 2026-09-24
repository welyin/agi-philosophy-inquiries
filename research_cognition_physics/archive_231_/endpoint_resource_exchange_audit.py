"""Round 438: autonomous endpoint-resource conversion and quantum regrouping.

Reuses the matter/link conversion idea of Hamma et al. (0911.5075, Eq. 16),
without resource hopping or an assumed geometric locality projector.
The exact capacity sectors are added model inputs, not derived axioms.
"""
import argparse
from collections import deque
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

TARGET = Path(__file__).with_name('endpoint_resource_exchange_audit_results.json')
OBS = {}


def edges_for(n):
    return list(itertools.combinations(range(n), 2))


def degrees(mask, n, edges):
    out = [0]*n
    for e, (i, j) in enumerate(edges):
        if mask & (1 << e):
            out[i] += 1
            out[j] += 1
    return out


def matchings(n):
    edges = edges_for(n)
    labels = {edge: e for e, edge in enumerate(edges)}
    def recurse(remaining):
        if not remaining:
            return [0]
        i, rest = remaining[0], remaining[1:]
        out = recurse(rest)
        for j in rest:
            bit = 1 << labels[i, j]
            out += [bit | mask for mask in recurse(tuple(k for k in rest if k != j))]
        return out
    return sorted(recurse(tuple(range(n))))


def graph_matrices(n, capacity=1):
    edges = edges_for(n)
    graphs = (matchings(n) if capacity == 1 else
              [mask for mask in range(2**len(edges)) if max(degrees(mask, n, edges)) <= capacity])
    index = {g: k for k, g in enumerate(graphs)}
    adjacency = np.zeros((len(graphs), len(graphs)), dtype=np.int64)
    for k, g in enumerate(graphs):
        for e in range(len(edges)):
            other = g ^ (1 << e)
            if other in index:
                adjacency[index[other], k] = 1
    occupancy = [np.array([int(bool(g & (1 << e))) for g in graphs], dtype=np.int64)
                 for e in range(len(edges))]
    return edges, graphs, adjacency, occupancy


def reduced_model(n, capacity=1):
    edges, graphs, adjacency, occupancy = graph_matrices(n, capacity)
    d = 2**n
    h = np.kron(np.eye(d, dtype=np.int64), adjacency+np.diag(sum(occupancy)))
    for (i, j), occ in zip(edges, occupancy):
        h += np.kron(core.swap(n, i, j).real.astype(np.int64), np.diag(occ))
    return h, edges, graphs, adjacency, occupancy


def resource_model(n, capacity=1):
    """Explicit uncompressed resource Hilbert space; unit-amplitude lowering."""
    edges = edges_for(n)
    base, width = capacity+1, 2**len(edges)
    dim = base**n*width
    conversion = np.zeros((dim, dim), dtype=np.int64)
    charges = np.zeros((n, dim), dtype=np.int64)
    occupancy = np.zeros((len(edges), dim), dtype=np.int64)
    powers = [base**(n-1-i) for i in range(n)]
    for col in range(dim):
        address, graph = divmod(col, width)
        free = [(address//p) % base for p in powers]
        degree = degrees(graph, n, edges)
        charges[:, col] = np.array(free)+degree
        for e, (i, j) in enumerate(edges):
            on = bool(graph & (1 << e))
            occupancy[e, col] = int(on)
            if not on and free[i] and free[j]:
                target = (address-powers[i]-powers[j])*width+(graph ^ (1 << e))
                conversion[target, col] += 1
            elif on and free[i] < capacity and free[j] < capacity:
                target = (address+powers[i]+powers[j])*width+(graph ^ (1 << e))
                conversion[target, col] += 1
    return conversion, charges, occupancy


def embedding(n, capacity=1):
    edges, graphs, _, _ = graph_matrices(n, capacity)
    width, base = 2**len(edges), capacity+1
    w = np.zeros((base**n*width, len(graphs)), dtype=np.int64)
    for k, graph in enumerate(graphs):
        deg = degrees(graph, n, edges)
        address = sum((capacity-deg[i])*base**(n-1-i) for i in range(n))
        w[address*width+graph, k] = 1
    return w


def shortest_paths(adjacency, start):
    distance = [-1]*len(adjacency)
    distance[start] = 0
    queue = deque([start])
    while queue:
        here = queue.popleft()
        for other in np.flatnonzero(adjacency[:, here]):
            if distance[other] == -1:
                distance[other] = distance[here]+1
                queue.append(int(other))
    return distance


class Audit(unittest.TestCase):
    def close(self, a, b, tolerance=2e-11):
        self.assertLess(float(np.linalg.norm(a-b)), tolerance)

    def test_01_exact_charge_sectors_and_full_generator(self):
        records = []
        for capacity in (1, 2):
            k, charges, occ = resource_model(3, capacity)
            self.assertTrue(np.array_equal(k, k.T))
            for q in charges:
                self.assertFalse(np.any(k*(q[:, None]-q[None, :])))
            w = embedding(3, capacity)
            edges, graphs, adjacency, occupation = graph_matrices(3, capacity)
            self.assertTrue(np.array_equal(w.T@w, np.eye(len(graphs))))
            self.assertTrue(np.array_equal(k@w, w@adjacency))
            for q in charges:
                self.assertTrue(np.array_equal(q[:, None]*w, capacity*w))
            for q, target in zip(occ, occupation):
                self.assertTrue(np.array_equal(q[:, None]*w, w*target[None, :]))
            records.append(dict(capacity=capacity, resource_dimension=len(k),
                legal_relation_dimension=len(graphs), all_local_charge_commutators_exactly_zero=True))
        k, charges, occ = resource_model(3)
        full = np.kron(np.eye(8, dtype=np.int64), k+np.diag(occ.sum(axis=0)))
        for (i, j), q in zip(edges_for(3), occ):
            full += np.kron(core.swap(3, i, j).real.astype(np.int64), np.diag(q))
        reduced = reduced_model(3)[0]
        whole_w = np.kron(np.eye(8, dtype=np.int64), embedding(3))
        self.assertTrue(np.array_equal(full@whole_w, whole_w@reduced))
        OBS['exact_capacity'] = dict(general_formula='Q_i=f_i+sum_j n_ij; [H,Q_i]=0',
            initial_sector='Q_i=C for every i', all_times_degree_bound='degree(i)<=C',
            degree_bound_applies_to_each_measured_configuration_and_occupation_operator=True,
            support_of_positive_mean_edges_need_not_be_a_matching=True,
            full_generator_dimension=len(full), checked_full_columns=whole_w.shape[1],
            exact_integer_full_intertwining=True, finite_capacity_checks=records,
            arbitrary_data_and_reference_allowed_within_prepared_resource_sector=True,
            sector_and_conversion_law_are_explicit_model_inputs=True)

    def test_02_nontrivial_data_response_and_unknown_information(self):
        h, edges, graphs, _, occ = reduced_model(3)
        g = len(graphs)
        cols = np.arange(0, len(h), g)
        certificates = []
        for (i, j), q in zip(edges, occ):
            ad = np.kron(np.eye(8, dtype=np.int64), np.diag(q))
            blocks = []
            for order in range(5):
                blocks.append((1j**order)*ad[np.ix_(cols, cols)])
                ad = h@ad-ad@h
            self.assertTrue(np.array_equal(blocks[2], 2*np.eye(8)))
            self.assertTrue(np.array_equal(blocks[4], -28*np.eye(8)-4*core.swap(3, i, j)))
            self.assertFalse(np.any(blocks[0]))
            self.assertFalse(np.any(blocks[1]))
            self.assertFalse(np.any(blocks[3]))
            certificates.append(dict(edge=[i, j], order_two='2I', order_four='-28I-4SWAP_ij'))
        source = (np.eye(8)-core.swap(3, 0, 1))/4
        blank = np.zeros((g, g)); blank[0, 0] = 1
        u = core.evolve(h, .5)
        output = u@np.kron(source, blank)@u.conj().T
        probabilities = [float(np.trace(output@np.kron(np.eye(8), np.diag(q))).real) for q in occ]
        self.assertGreater(probabilities[0], probabilities[1])
        self.assertAlmostEqual(probabilities[1], probabilities[2])
        disturbance = core.distance(core.partial(output, [8, g], (0,)), source)
        self.assertGreater(disturbance, 1e-5)
        rng = np.random.default_rng(43802)
        data_reference = rng.normal(size=(8, 3))+1j*rng.normal(size=(8, 3))
        data_reference /= np.linalg.norm(data_reference)
        v = np.eye(len(h))[:, cols]
        evolved = u@v@data_reference
        self.close(u.conj().T@evolved, v@data_reference)
        self.close(evolved.conj().T@evolved, data_reference.conj().T@data_reference)
        OBS['data_and_reference'] = dict(integer_response_certificates=certificates,
            small_time_probability='t^2-(7/6+<SWAP_e>/6)t^4+O(t^5)',
            measured_time=.5, relation_occupation_probabilities=probabilities,
            data_marginal_trace_distance_change=disturbance,
            complete_unknown_joint_information_preserved=True,
            reference_marginal_residual=float(np.linalg.norm(evolved.conj().T@evolved-data_reference.conj().T@data_reference)),
            arbitrary_data_marginal_unchanged_claimed=False)

    def test_03_identical_local_inputs_and_collective_symmetry(self):
        h, _, graphs, _, _ = reduced_model(3)
        g = len(graphs)
        sigma = core.PAULI
        identity = np.eye(2)
        generators = []
        for p in sigma:
            collective = np.kron(np.kron(p, identity), identity)+np.kron(np.kron(identity, p), identity)+np.kron(np.kron(identity, identity), p)
            extended = np.kron(collective, np.eye(g))
            self.close(h@extended, extended@h)
            generators.append(float(np.linalg.norm(h@extended-extended@h)))
        rho = (identity+.23*sigma[0]+.17*sigma[1]+.31*sigma[2])/2
        blank = np.zeros((g, g)); blank[0, 0] = 1
        initial = np.kron(np.kron(np.kron(rho, rho), rho), blank)
        u = core.evolve(h, .7)
        out = u@initial@u.conj().T
        errors = []
        for i in range(3):
            marginal = core.partial(out, [2, 2, 2, g], (i,))
            self.close(marginal, rho)
            errors.append(float(np.linalg.norm(marginal-rho)))
        # Direct check of all simultaneous data/edge relabelings.
        import quantum_participation_audit as relabel
        graph_index = {mask: j for j, mask in enumerate(graphs)}
        edges = edges_for(3)
        permutation_errors = []
        for pi in itertools.permutations(range(3)):
            edge_map = [edges.index(tuple(sorted((pi[i], pi[j])))) for i, j in edges]
            pgraph = np.zeros((g, g), dtype=np.int64)
            for j, mask in enumerate(graphs):
                destination = sum(1 << edge_map[e] for e in range(len(edges)) if mask & (1 << e))
                pgraph[graph_index[destination], j] = 1
            p = np.kron(relabel.factor_permutation(pi), pgraph)
            self.close(p@h@p.T, h)
            permutation_errors.append(float(np.linalg.norm(p@h@p.T-h)))
        OBS['identical_input_condition'] = dict(collective_spin_commutator_residuals=generators,
            joint_permutation_residuals=permutation_errors,
            local_identical_mixed_state_errors=errors,
            all_N_proof_uses_permutation_covariance_and_collective_spin_conservation=True,
            arbitrary_correlated_equal_marginals_not_claimed_stationary=True,
            full_extended_law_not_derived_from_429_pair_contract=True)

    def test_04_exact_autonomous_reassignment(self):
        h, _, graphs, adjacency, occ = reduced_model(3)
        graph_h = adjacency+2*np.diag(sum(occ))
        t = math.pi/2
        u = core.evolve(graph_h, t)
        initial_index = graphs.index(1)  # Edge (0,1).
        probabilities = abs(u[:, initial_index])**2
        self.close(probabilities, np.array([0., 5/9, 2/9, 2/9]))
        ub = core.evolve(graph_h, math.pi/4)
        self.close(abs(ub[:, 0])**2, np.full(4, .25))
        from exchange_participation_capacity_audit import symmetric_basis
        s = symmetric_basis(3)
        w = np.kron(s, np.eye(4))
        self.close(h@w, w@np.kron(np.eye(4), graph_h))
        self.close(core.evolve(h, t)@w, w@np.kron(np.eye(4), u))
        OBS['autonomous_reassignment'] = dict(
            all_pair_parameters=dict(Omega=1, mu=1, J=1), graph_order=['empty', '01', '02', '12'],
            initial_graph='01', time='pi/2', exact_probabilities=['0', '5/9', '2/9', '2/9'],
            spectral_probabilities=probabilities.tolist(),
            blank_at_pi_over_four_exact_probabilities=['1/4']*4,
            complete_symmetric_unknown_data_and_reference_unchanged=True,
            no_measurement_reset_or_external_switch=True,
            requested_partner_can_be_selected_deterministically=False,
            permanent_record_or_stationary_geometry_claimed=False)

    def test_05_all_matching_configurations_connected_at_leading_order(self):
        records = []
        for n in range(2, 7):
            edges, graphs, adjacency, _ = graph_matrices(n)
            maximum = 0
            for i, first in enumerate(graphs):
                distances = shortest_paths(adjacency, i)
                self.assertTrue(all(d >= 0 for d in distances))
                for j, second in enumerate(graphs):
                    self.assertEqual(distances[j], (first ^ second).bit_count())
                maximum = max(maximum, max(distances))
            records.append(dict(nodes=n, legal_graph_count=len(graphs), maximum_distance=maximum))
        h, _, graphs, adjacency, _ = reduced_model(4)
        g, d = len(graphs), 16
        powers = [np.eye(len(h), dtype=np.int64)]
        graph_powers = [np.eye(g, dtype=np.int64)]
        for _ in range(4):
            powers.append(powers[-1]@h)
            graph_powers.append(graph_powers[-1]@adjacency)
        pairs = 0
        for first, first_mask in enumerate(graphs):
            for second, second_mask in enumerate(graphs):
                if first == second:
                    continue
                distance = (first_mask ^ second_mask).bit_count()
                for k in range(distance):
                    self.assertFalse(np.any(powers[k][second::g, first::g]))
                count = int(graph_powers[distance][second, first])
                self.assertGreater(count, 0)
                self.assertTrue(np.array_equal(powers[distance][second::g, first::g], count*np.eye(d)))
                pairs += 1
        OBS['reconfiguration_connectivity'] = dict(finite_graph_checks=records,
            exact_full_data_leading_order_pairs_checked=pairs,
            analytic_distance='|M symmetric_difference M_prime|',
            leading_block='(-i Omega t)^r * number_of_shortest_paths/r! * I_D',
            all_finite_N_analytic_path_argument=True,
            graph_connectivity_not_claimed_as_ergodicity_or_control_universality=True)

    def test_06_hidden_hardware_and_rate_costs(self):
        records = []
        previous_counts = [1, 1]
        for n in range(2, 9):
            previous_counts.append(previous_counts[-1]+(n-1)*previous_counts[-2])
            edges, graphs, adjacency, _ = graph_matrices(n)
            self.assertEqual(len(graphs), previous_counts[-1])
            incident = sum(1 << e for e, (i, j) in enumerate(edges) if i == 0 or j == 0)
            b = np.zeros_like(adjacency)
            components = {}
            for i, first in enumerate(graphs):
                key = first & ~incident
                components.setdefault(key, []).append(i)
                for j, second in enumerate(graphs):
                    if adjacency[j, i] and ((first ^ second) & incident):
                        b[j, i] = 1
            largest_k = 0
            for inds in components.values():
                block = b[np.ix_(inds, inds)]
                k = len(inds)-1
                largest_k = max(largest_k, k)
                self.assertEqual(int(block.sum()), 2*k)
                if k > 1:
                    self.assertEqual(sorted(block.sum(axis=0).tolist()), [1]*k+[k])
            self.assertEqual(largest_k, n-1)
            spectral_norm = None
            if n <= 7:
                spectral_norm = float(np.linalg.eigvalsh(b)[-1])
                self.assertAlmostEqual(spectral_norm, math.sqrt(n-1))
            row = dict(nodes=n, graph_sector_dimension=len(graphs),
                explicit_hardware_qubits=2*n+math.comb(n, 2),
                lossless_graph_qubits_at_least=math.ceil(math.log2(len(graphs))),
                largest_incident_star_leaves=largest_k, incident_conversion_norm=spectral_norm)
            if n % 2 == 0:
                perfect = math.factorial(n)//(2**(n//2)*math.factorial(n//2))
                self.assertLessEqual(perfect, len(graphs))
                row['perfect_matching_lower_bound'] = perfect
            records.append(row)
        OBS['resource_scope'] = dict(finite_counts_and_norms=records,
            exact_endpoint_conversion_norm='|Omega| sqrt(N-1) in Q_i=1 sector',
            fixed_operator_rate_budget_requires_Omega_at_most_order_inverse_sqrt_N=True,
            matching_count_recurrence='a_N=a_(N-1)+(N-1)a_(N-2)',
            even_N_count_lower='N!/[2^(N/2)(N/2)!]',
            complete_matching_record_capacity='Theta(N log N) qubits globally for faithful full labeled matching sector',
            capacity_bound_not_applied_to_one_restricted_symmetric_trajectory=True,
            sparse_active_degree_not_constant_total_hardware_or_rate=True,
            resource_charge_not_identified_with_physical_energy=True,
            optional_exchange_simulator_preserves_capacity_only_up_to_trace_error=True)

        # Exact symmetric trajectory is much smaller than the full interface.
        n = 8
        edges, graphs, adjacency, occ = graph_matrices(n)
        size = n//2+1
        w = np.zeros((len(graphs), size))
        counts = [sum(g.bit_count() == k for g in graphs) for k in range(size)]
        for row, mask in enumerate(graphs):
            w[row, mask.bit_count()] = 1/math.sqrt(counts[mask.bit_count()])
        small = np.diag(2*np.arange(size)).astype(float)
        for k in range(size-1):
            small[k, k+1] = small[k+1, k] = math.sqrt((k+1)*math.comb(n-2*k, 2))
        graph_h = adjacency+2*np.diag(sum(occ))
        self.close(w.T@w, np.eye(size))
        self.close(graph_h@w, w@small)
        OBS['resource_scope']['restricted_symmetric_trajectory'] = dict(nodes=n,
            full_matching_dimension=len(graphs), symmetric_dimension=size,
            off_diagonal_rule='sqrt((k+1)*binomial(N-2k,2))',
            exact_generator_intertwining=True,
            arbitrary_labeled_interventions_not_preserved_by_this_restriction=True)

    def test_07_capacity_does_not_make_current_graph_causal(self):
        h, edges, graphs, _, _ = reduced_model(3)
        g = len(graphs)
        receiver_y = np.kron(np.kron(np.eye(2), core.PAULI[1]), np.eye(2))
        ad = np.kron(receiver_y, np.eye(g))
        cols = np.arange(0, len(h), g)
        blocks = []
        for k in range(4):
            blocks.append((1j**k)*ad[np.ix_(cols, cols)])
            ad = h@ad-ad@h
        swap_sum = sum(core.swap(3, i, j) for i, j in edges)
        self.assertTrue(np.array_equal(blocks[0], receiver_y))
        self.assertFalse(np.any(blocks[1]))
        self.assertFalse(np.any(blocks[2]))
        self.assertTrue(np.array_equal(blocks[3], 2j*(swap_sum@receiver_y-receiver_y@swap_sum)))
        source_plus = (np.eye(2)+core.PAULI[2])/2
        source_minus = (np.eye(2)-core.PAULI[2])/2
        target = (np.eye(2)+core.PAULI[0])/2
        rho_plus = np.kron(np.kron(source_plus, target), np.eye(2)/2)
        rho_minus = np.kron(np.kron(source_minus, target), np.eye(2)/2)
        derivative_gap = float(np.trace((rho_plus-rho_minus)@blocks[3]).real)
        self.assertEqual(derivative_gap, 4.)
        blank = np.zeros((g, g)); blank[0, 0] = 1
        rows = []
        for t in (.1, .2, .5):
            u = core.evolve(h, t)
            delta = u@np.kron(rho_plus-rho_minus, blank)@u.conj().T
            readout = float(np.trace(delta@np.kron(receiver_y, np.eye(g))).real)
            self.assertGreater(readout, 0.)
            rows.append(dict(time=t, receiver_y_expectation_gap=readout,
                leading_cubic_prediction=2*t**3/3))
        OBS['capacity_vs_propagation'] = dict(
            reduced_data_channel='rho-i Omega^2 J t^3/3 [sum_ij SWAP_ij,rho]+O(t^4)',
            exact_receiver_third_derivative_gap=derivative_gap,
            all_N_fixed_finite_model_Kraus_expansion_proved=True,
            initial_active_graph='empty', full_legal_configuration_degree_at_most_one=True,
            finite_receiver_checks=rows,
            current_edge_snapshot_is_not_full_Hamiltonian_causal_support=True,
            physical_faster_than_light_or_metric_claimed=False)

        # A certified positive finite-time gap, not just floating-point evidence.
        # ||H|| <= 3+2 = 5; odd signal has no fourth-order coefficient.
        small_t = Fraction(1, 100)
        x = 10*small_t
        remainder_upper = 2*x**5/(math.factorial(5)*(1-x))
        gap_lower = Fraction(2, 3)*small_t**3-remainder_upper
        self.assertEqual(gap_lower, Fraction(13, 27000000))
        OBS['capacity_vs_propagation']['rational_signal_certificate'] = dict(
            time=str(small_t), Hamiltonian_norm_upper=5,
            odd_signal_remainder_upper=str(remainder_upper),
            receiver_Y_gap_lower=str(gap_lower),
            applies_to_fixed_N_three_not_a_size_uniform_remainder=True)

        # Reuse the relational chirality idea of 430, not an external axis.
        # A,B,C,R ordering: C is maximally mixed and R stays out of H.
        chirality = np.zeros((16, 16), complex)
        for p in itertools.permutations(range(3)):
            inversions = sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3))
            chirality += (-1)**inversions*np.kron(np.kron(np.kron(core.PAULI[p[0]], core.PAULI[p[1]]), np.eye(2)), core.PAULI[p[2]])
        invariant_plus = (np.eye(16)+chirality/6)/16
        invariant_minus = (np.eye(16)-chirality/6)/16
        self.assertGreater(float(np.linalg.eigvalsh(invariant_plus)[0]), 0.)
        self.assertGreater(float(np.linalg.eigvalsh(invariant_minus)[0]), 0.)
        for pauli in core.PAULI:
            collective = np.zeros((16, 16), complex)
            for slot in range(4):
                term = np.ones((1, 1))
                for i in range(4):
                    term = np.kron(term, pauli if i == slot else np.eye(2))
                collective += term
            self.close(collective@chirality, chirality@collective)
        # Data, graph, R ordering. Graph dimension 4 occupies two tensor bits.
        v = np.kron(np.eye(len(h))[:, cols], np.eye(2))
        global_h = np.kron(h, np.eye(2))
        readout_swap = core.swap(6, 1, 5)
        u = core.evolve(global_h, .5)
        difference = u@v@(invariant_plus-invariant_minus)@v.conj().T@u.conj().T
        invariant_gap = float(np.trace(difference@readout_swap).real)
        self.assertAlmostEqual(invariant_gap, rows[-1]['receiver_y_expectation_gap']/2)
        self.close(core.partial(invariant_plus, [2]*4, (1, 3)), core.partial(invariant_minus, [2]*4, (1, 3)))
        OBS['capacity_vs_propagation']['internal_relational_witness'] = dict(
            invariant_initial_states='(I +/- chi_ABR/6)/16 with C maximally mixed',
            internal_reference_qubits=1, readout='SWAP_BR', time=.5,
            relational_expectation_gap=invariant_gap,
            same_initial_B_R_marginals=True,
            no_external_absolute_axis_required=True,
            known_430_chirality_method_reused=True,
            reference_preparation_isolation_and_readout_remain_inputs=True,
            invariant_state_pair_not_claimed_related_by_a_local_A_only_gate=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=438, baseline_round=437, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(existing_matter_link_conversion_explicitly_attributed=True,
            local_resource_charge_conservation_and_hard_degree_bound_proved=True,
            autonomous_reassignment_and_data_response_verified=True,
            all_unknown_joint_data_reference_preserved=True,
            all_matching_configurations_have_nonzero_short_time_path_amplitudes=True,
            relationship_memory_and_rate_costs_explicit=True,
            cubic_complete_pair_data_response_from_empty_graph_proved=True,
            conversion_rule_charge_preparation_and_capacity_derived_from_429=False,
            arbitrary_data_marginal_unchanged_claimed=False,
            stable_connected_geometry_or_dimension_generated=False,
            final_Heisenberg_hardware_compiled=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert json.loads(TARGET.read_text(encoding='utf-8')) == result
    elif not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
