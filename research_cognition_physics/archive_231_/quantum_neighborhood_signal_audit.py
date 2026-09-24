"""Round 443: physical subject readout in the unchanged six-tree model.

The obsolete graph trace counterexample reproduces a problem acknowledged in
Arrighi--Durbec--Wilson (2024).  The new finite certificate concerns actual
data transmission, including every coherent far-graph input and its reference.
"""
import argparse
from fractions import Fraction as Q
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import branching_tree_distance_audit as old

TARGET = Path(__file__).with_name('quantum_neighborhood_signal_audit_results.json')
OBS = {}
W2 = np.array([[1, 1, 1, 1], [1, -1, 1, -1],
               [1, 1, -1, -1], [1, -1, -1, 1]], dtype=np.int64).T


def bits(z, n=6):
    return tuple((z >> (n-1-i)) & 1 for i in range(n))


def ball(tree, center, radius, n=6):
    adj = old.adjacency(tree, n)
    region = {center}
    for _ in range(radius):
        region |= set().union(*(adj[v] for v in region))
    return frozenset(region)


def port_edges(tree):
    """A fixed valid port assignment for the two legacy counterexample trees."""
    def port(v, w):
        if v not in (1, 2):
            return 0
        if w in (1, 2):
            return 1
        if w in (3, 4):
            return 2
        return 0
    return tuple(sorted((a, port(a, b), b, port(b, a)) for a, b in tree))


def legacy_disk(tree, region):
    adj = old.adjacency(tree, 6)
    border = set().union(*(adj[v] for v in region))-set(region)
    vertices = tuple(sorted(set(region)|border))
    edges = tuple(e for e in port_edges(tree) if e[0] in vertices and e[2] in vertices)
    labels = tuple((v, 0 if v in region else None) for v in vertices)
    return vertices, edges, labels


def legacy_complement(tree, region):
    vertices = tuple(v for v in range(6) if v not in region)
    edges = tuple(e for e in port_edges(tree) if e[0] in vertices and e[2] in vertices)
    return vertices, edges, tuple((v, 0) for v in vertices)


def split_record(tree, z, center=5, radius=1):
    """Each physical edge qubit goes to exactly one side, including zero bits."""
    region = ball(tree, center, radius)
    data = [('D', v, bits(z)[v]) for v in range(6)]
    edges = [('E', a, b, int((a, b) in tree))
             for a, b in itertools.combinations(range(6), 2)]
    local = tuple(x for x in data if x[1] in region) + tuple(
        x for x in edges if x[1] in region or x[2] in region)
    external = tuple(x for x in data if x[1] not in region) + tuple(
        x for x in edges if x[1] not in region and x[2] not in region)
    return local, external


def sector():
    trees, graphs, full, flip, distance = old.six_vertex_sector()
    ids = np.array([(1 << (5-v))*6+g for v in range(6) for g in range(6)])
    h = full[np.ix_(ids, ids)]
    far = np.flatnonzero(distance == 3)
    output = np.arange(30, 36)
    return trees, full, flip, ids, h, far, output


def block_powers(h, far, output, maximum):
    power = np.eye(36, dtype=object)[:, far]
    result = []
    for k in range(maximum+1):
        result.append(power[output].copy())
        power = h.astype(object) @ power
    return result


def exact_taylor_certificate(h, far, output, degree=40):
    powers = block_powers(h, far, output, degree)
    real = np.full((6, 4), Q(0), dtype=object)
    imag = real.copy()
    for k, block in enumerate(powers):
        term = (block @ W2) * Q(1, 2*math.factorial(k))
        if k % 4 == 0:
            real += term
        elif k % 4 == 1:
            imag -= term
        elif k % 4 == 2:
            real -= term
        else:
            imag += term
    gram_real = real.T @ real + imag.T @ imag
    gram_imag = real.T @ imag - imag.T @ real
    tail = Q(9**(degree+1), math.factorial(degree+1)) / (1-Q(9, degree+2))
    return gram_real, gram_imag, tail


class Audit(unittest.TestCase):
    def close(self, a, b, tol=3e-12):
        self.assertLess(float(np.max(np.abs(np.asarray(a)-np.asarray(b)))), tol)

    def test_01_known_legacy_trace_failure_in_current_sector(self):
        g = frozenset({(0, 1), (1, 2), (2, 5), (1, 3), (2, 4)})
        h = frozenset({(0, 2), (1, 2), (1, 5), (1, 3), (2, 4)})
        trees, *_ = sector()
        self.assertIn(g, trees)
        self.assertIn(h, trees)
        for tree in (g, h):
            ports = [(a, p) for a, p, b, q in port_edges(tree)]
            ports += [(b, q) for a, p, b, q in port_edges(tree)]
            self.assertEqual(len(ports), len(set(ports)))
        terms = []
        for left in (g, h):
            for right in (g, h):
                region = ball(left, 0, 1)|ball(right, 0, 1)
                if legacy_complement(left, region) == legacy_complement(right, region):
                    terms.append((legacy_disk(left, region), legacy_disk(right, region)))
        basis = list(dict.fromkeys(x for pair in terms for x in pair))
        self.assertEqual(len(basis), 4)
        twice = np.zeros((4, 4), dtype=np.int64)
        for x, y in terms:
            twice[basis.index(x), basis.index(y)] += 1
        self.assertEqual(sorted(np.linalg.eigvalsh(twice).tolist()), [-1, 1, 1, 1])
        self.assertEqual(int(np.trace(twice)), 2)
        OBS['known_2017_trace_failure'] = dict(
            attributed_to_2024_correction=True, radius=1,
            fixed_port_assignment_verified=True,
            output_eigenvalues=['-1/2', '1/2', '1/2', '1/2'],
            new_discovery_claimed=False)

    def test_02_branchwise_isometry_and_subject_readout_all_matrix_units(self):
        trees, *_ = sector()
        records = [split_record(tree, z) for z in range(64) for tree in trees]
        self.assertEqual(len(set(records)), 384)  # V*V=I exactly.
        for z in range(64):
            for g, tree in enumerate(trees):
                local, outside = records[6*z+g]
                self.assertEqual(len(local)+len(outside), 21)
                self.assertFalse(set(local)&set(outside))
                self.assertEqual(sorted(x for x in local+outside if x[0] == 'D'),
                                 [('D', v, bits(z)[v]) for v in range(6)])
                self.assertEqual(len([x for x in local+outside if x[0] == 'E']), 15)
        stripped = [tuple(x for x in local if x[:2] != ('D', 5))
                    for local, outside in records]
        retained = 0
        for x, (lx, ex) in enumerate(records):
            zx, gx = divmod(x, 6)
            for y, (ly, ey) in enumerate(records):
                zy, gy = divmod(y, 6)
                via_split = ex == ey and stripped[x] == stripped[y]
                direct = gx == gy and zx//2 == zy//2
                self.assertEqual(via_split, direct)
                retained += int(via_split)
        OBS['finite_neighborhood_channel'] = dict(
            input_dimension=384, injective_basis_pairs=384,
            local_basis_dimension=len(set(l for l, e in records)),
            environment_basis_dimension=len(set(e for l, e in records)),
            exact_readout_matrix_units_checked=384**2,
            nonzero_subject_matrix_units=retained,
            data_qubits=6, potential_edge_qubits=15,
            implicit_global_routing_cost_claimed_zero=False,
            arbitrary_reference_covered_by_isometry_proof=True)

    def test_03_one_excitation_block_of_unchanged_full_H(self):
        trees, full, flip, ids, h, far, output = sector()
        embedding = np.eye(384, dtype=np.int64)[:, ids]
        self.assertTrue(np.array_equal(full @ embedding, embedding @ h))
        expected = np.kron(np.eye(6, dtype=np.int64), flip)
        for g, tree in enumerate(trees):
            a = np.zeros((6, 6), dtype=np.int64)
            for v, w in tree:
                a[v, w] = a[w, v] = 1
            laplacian = np.diag(a.sum(axis=1))-a
            index = np.arange(g, 36, 6)
            expected[np.ix_(index, index)] += 5*np.eye(6, dtype=np.int64)-laplacian
        self.assertTrue(np.array_equal(h, expected))
        ground = np.eye(384, dtype=np.int64)[:, :6]
        self.assertTrue(np.array_equal(full @ ground,
            ground @ (flip+5*np.eye(6, dtype=np.int64))))
        self.assertTrue(np.array_equal(far, [0, 1, 4, 5]))
        self.assertTrue(np.array_equal(flip.sum(axis=1), np.full(6, 4)))
        self.assertLessEqual(float(np.linalg.norm(h, 2)), 9+2e-12)
        OBS['unchanged_autonomous_model'] = dict(
            full_dimension=384, exact_invariant_one_excitation_dimension=36,
            graph_dimension=6, far_graph_indices=far.tolist(),
            J=1, kappa=1, omitted_scalar_mu_edge_phase=5,
            norm_bound=9, receiver_label=5, sender_label=0,
            preparation='all six data qubits zero; sender optionally applies X')

    def test_04_exact_coherent_signal_orders(self):
        trees, full, flip, ids, h, far, output = sector()
        blocks = block_powers(h, far, output, 5)
        for k in range(3):
            self.assertFalse(np.any(blocks[k]))
        m3 = np.array([[1, 0, 1, 0], [0, 1, 0, 1], [1, 1, 2, 2],
                       [2, 2, 1, 1], [1, 0, 1, 0], [0, 1, 0, 1]])
        self.assertTrue(np.array_equal(blocks[3], m3))
        self.assertTrue(np.array_equal(W2.T @ m3.T @ m3 @ W2,
                                        4*np.diag([22, 4, 2, 0])))
        self.assertFalse(np.any(blocks[4] @ W2[:, 3]))
        dark_m5 = blocks[5] @ W2[:, 3]
        self.assertTrue(np.array_equal(dark_m5, [-2, 2, 0, 0, 2, -2]))
        # Verify coefficients in noncommutative cubic expansion, not sampling J,kappa.
        a = h-np.kron(np.eye(6, dtype=np.int64), flip)
        f = h-a
        cubic = [np.zeros((6, 4), dtype=np.int64) for _ in range(4)]
        for word in itertools.product((0, 1), repeat=3):
            product = np.eye(36, dtype=np.int64)
            for letter in word:
                product = product @ (a if letter == 0 else f)
            cubic[sum(word)] += product[np.ix_(output, far)]
        self.assertFalse(np.any(cubic[2]))
        self.assertFalse(np.any(cubic[3]))
        self.assertTrue(np.array_equal(cubic[0], np.eye(6, dtype=np.int64)[:, far]))
        self.assertTrue(np.array_equal(cubic[0]+cubic[1], m3))
        e, b = cubic[:2]
        self.assertTrue(np.array_equal(W2.T @ e.T @ e @ W2, 4*np.eye(4, dtype=int)))
        self.assertTrue(np.array_equal(W2.T @ (e.T @ b+b.T @ e) @ W2,
                                        4*np.diag([2, 2, -2, -2])))
        self.assertTrue(np.array_equal(W2.T @ b.T @ b @ W2,
                                        4*np.diag([19, 1, 3, 1])))
        OBS['coherence_changes_actual_signal'] = dict(
            first_three_blocks_zero=True, cubic_block=m3.tolist(),
            cubic_gram_character_eigenvalues=[22, 4, 2, 0],
            mixed_far_leading_probability='7*t^6/36',
            uniform_coherent_far_leading_probability='11*t^6/18',
            dark_coherent_far_leading_probability='t^10/3600',
            first_nonzero_dark_amplitude_order=5,
            general_cubic_form='J^3*E + J^2*kappa*B; E=far inclusion, B=M3-E',
            general_cubic_gram_eigenvalues=[
                'J^4*(J^2+2*J*kappa+19*kappa^2)', 'J^4*(J+kappa)^2',
                'J^4*(J^2-2*J*kappa+3*kappa^2)', 'J^4*(J-kappa)^2'],
            dark_cubic_cancellation_requires_J_equal_kappa_except_trivial_J_zero=True)

    def test_05_exact_rational_all_graph_and_reference_signal_certificate(self):
        _, _, _, _, h, far, output = sector()
        real, imag, tail = exact_taylor_certificate(h, far, output)
        for i in range(4):
            for j in range(4):
                self.assertEqual(imag[i, j], 0)
                if i != j:
                    self.assertEqual(real[i, j], 0)
        minimum = min(real[i, i] for i in range(4))
        self.assertGreater(minimum, Q(101, 10000)**2)
        self.assertLess(tail, Q(1, 10**10))
        self.assertGreater(Q(101, 10000)-tail, Q(1, 100))
        OBS['all_far_graph_signal_certificate'] = dict(
            time='1', taylor_degree=40, norm_bound=9,
            exact_polynomial_gram_diagonal=True,
            polynomial_gram_diagonal=[float(real[i, i]) for i in range(4)],
            exact_minimum_polynomial_gram=str(minimum),
            exact_operator_remainder_bound=str(tail),
            remainder_upper_bound='1/10000000000',
            certified_receiver_excitation_probability_lower_bound='1/10000',
            certified_receiver_trace_distance_lower_bound='1/10000',
            unknown_graph_reference_dimension_independent=True,
            unknown_other_subject_data_covered=False,
            arbitrary_large_graphs_covered=False)

    def test_06_complete_subject_readout_with_entangled_internal_reference(self):
        _, full, _, ids, h, far, output = sector()
        eigen, vectors = np.linalg.eigh(h)
        table = []
        for time in (.1, 1.):
            u = (vectors*np.exp(-1j*eigen*time)) @ vectors.T
            k = u[np.ix_(output, far)]
            gram = k.conj().T @ k
            char = W2.T @ gram @ W2 / 4
            self.close(char, np.diag(np.diag(char)))
            table.append(dict(time=time,
                probability_character_states=np.diag(char).real.tolist(),
                probability_uniform_mixture=float(np.trace(gram).real/4),
                minimum_over_all_graph_inputs=float(np.linalg.eigvalsh(gram)[0])))
        eigen_full, vectors_full = np.linalg.eigh(full)
        u_full = (vectors_full*np.exp(-1j*eigen_full)) @ vectors_full.T
        states = []
        for message in (0, 1):
            psi = np.zeros((384, 4), dtype=complex)
            for r, g in enumerate(far):
                psi[message*32*6+g, r] = .5
            evolved = (u_full @ psi).reshape(32, 2, 6, 4)
            reduced = np.einsum('ebgr,ecgs->brcs', evolved, evolved.conj()).reshape(8, 8)
            self.close(np.trace(reduced), 1)
            self.assertGreater(float(np.linalg.eigvalsh(reduced)[0]), -3e-12)
            states.append(reduced)
        self.close(states[0][:4, :4], np.eye(4)/4)
        self.close(states[0][4:, 4:], 0)
        p = float(np.trace(states[1][4:, 4:]).real)
        delta = states[1]-states[0]
        trace_distance = float(np.sum(np.abs(np.linalg.eigvalsh(delta)))/2)
        self.close(p, trace_distance)
        self.close(p, table[-1]['probability_uniform_mixture'])
        self.assertGreater(table[-1]['minimum_over_all_graph_inputs'], 1e-4)
        OBS['numerical_readout'] = dict(
            table=table, internal_reference_dimension=4,
            full_data_graph_H_used=True,
            initial_graph_reference_maximally_entangled=True,
            receiver_reference_trace_distance=trace_distance,
            receiver_excitation_probability=p,
            no_graph_postselection=True)


def run():
    OBS.clear()
    buffer = io.StringIO()
    result = unittest.TextTestRunner(stream=buffer, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(buffer.getvalue())
    return dict(round=443, baseline_round=442, date='2026-09-24',
        runtime=dict(python=platform.python_version(), numpy=np.__version__),
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        observations=OBS.copy(),
        scope=dict(
            unchanged_440_442_autonomous_H=True,
            physical_receiver_data_signal_certified=True,
            all_coherent_far_graph_inputs_and_reference_covered=True,
            finite_channel_and_actual_subject_readout_compatible=True,
            obsolete_trace_failure_attributed_to_2024=True,
            exact_rational_operator_certificate=True,
            universal_branching_graph_propagation_bound=False,
            unknown_other_subject_data_preserved=False,
            graph_neighborhood_routing_autonomously_implemented=False,
            nested_dynamic_ball_channels_assumed=False,
            new_internal_measurement_or_encoding_device_constructed=False,
            initial_geometry_selected=False,
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
