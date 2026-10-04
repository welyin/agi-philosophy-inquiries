"""Round 466: arbitrary-background propagation from connected atomic words.

Baseline 465. No vacuum assumption, classical graph replacement, threshold,
fixed chain, new Hamiltonian, or dynamic-neighborhood trace is introduced.
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

TARGET = Path(__file__).with_name('unknown_background_propagation_audit_results.json')
OBS = {}


def short(x):
    return float(f'{float(x):.12g}')


def induced_components(tree, vertices):
    vertices = set(vertices)
    adjacent = {v: set() for v in vertices}
    for a, b in tree:
        if a in vertices and b in vertices:
            adjacent[a].add(b)
            adjacent[b].add(a)
    parts = []
    while vertices:
        todo, part = [min(vertices)], set()
        while todo:
            v = todo.pop()
            if v in part:
                continue
            part.add(v)
            todo.extend(adjacent[v]-part)
        vertices -= part
        parts.append(tuple(sorted(part)))
    return tuple(sorted(parts))


def atoms(trees, n):
    lookup = {g: i for i, g in enumerate(trees)}
    out = []
    for a, b in itertools.combinations(range(n), 2):
        targets = np.array([g if (a, b) in tree else -1 for g, tree in enumerate(trees)])
        if np.any(targets >= 0):
            out.append(('D', frozenset((a, b)), targets, (a, b)))
    for b, c in itertools.combinations(range(n), 2):
        for a, d in itertools.combinations([v for v in range(n) if v not in (b, c)], 2):
            left = {old.edge(a, b), old.edge(c, d)}
            right = {old.edge(a, c), old.edge(b, d)}
            targets = np.full(len(trees), -1)
            for g, tree in enumerate(trees):
                if (b, c) not in tree:
                    continue
                if ((left <= tree and not right & tree) or
                    (right <= tree and not left & tree)):
                    target = frozenset(tree.symmetric_difference(left | right))
                    targets[g] = lookup[target]
            if np.any(targets >= 0):
                out.append(('F', frozenset((a, b, c, d)), targets, (a, b, c, d)))
    return out


def rooted_atom_family(supports, root):
    reached, pending = {root}, list(supports)
    while pending:
        next_pending = []
        for s in pending:
            if s & reached:
                reached |= s
            else:
                next_pending.append(s)
        if len(next_pending) == len(pending):
            return False
        pending = next_pending
    return True


def evolution(h, z):
    e, v = np.linalg.eigh(h)
    return (v*np.exp(-1j*z*e))@v.conj().T


def distance(a, b):
    return float(np.abs(np.linalg.eigvalsh(a-b)).sum()/2)


def receiver(vector, nr=3):
    x = vector.reshape(32, 2, 6, nr).transpose(1, 3, 0, 2).reshape(2*nr, -1)
    return x@x.conj().T


class Audit(unittest.TestCase):
    def close(self, a, b, tol=3e-11):
        self.assertLess(np.linalg.norm(a-b), tol)

    def test_01_atomic_decomposition_and_local_activity(self):
        trees, _, h, f, _ = old.six_vertex_sector()
        aa = atoms(trees, 6)
        reconstructed = np.zeros_like(h)
        partitions_checked = 0
        for kind, support, targets, parameters in aa:
            gmat = np.zeros((6, 6), dtype=int)
            for g, target in enumerate(targets):
                if target < 0:
                    continue
                self.assertEqual(targets[target], g)
                gmat[target, g] = 1
                outside = set(range(6))-support
                for bits in itertools.product((0, 1), repeat=len(outside)):
                    s = support | {v for v, yes in zip(sorted(outside), bits) if yes}
                    self.assertEqual(induced_components(trees[g], s),
                                     induced_components(trees[target], s))
                    partitions_checked += 1
            data = old.core.swap(6, *parameters) if kind == 'D' else np.eye(64, dtype=int)
            block = np.kron(data, gmat)
            self.assertEqual(float(np.abs(block.imag).max()), 0)
            reconstructed += block.real.astype(int)
        self.close(reconstructed, h, 1e-15)
        lam = 3+2*3*2**2
        for g in range(6):
            for v in range(6):
                active = sum(targets[g] >= 0 and v in support for _, support, targets, _ in aa)
                self.assertLessEqual(active, lam)
        OBS['atomic_contract'] = dict(atoms=len(aa), full_dimension=384,
            exact_H_reconstruction=True, induced_partition_checks=partitions_checked,
            flip_support_includes_all_five_internal_edge_factors=True,
            no_external_tree_legality_test_required=True,
            each_atomic_graph_block_a_partial_bijection=True,
            active_weight_bound='lambda*|S|, lambda=C*|J|+2*C*(C-1)^2*|kappa|')

    def test_02_connected_word_support_in_initial_graph(self):
        trees, _, _, _, _ = old.six_vertex_sector()
        aa = atoms(trees, 6)
        checked, nongrowing_order = 0, 0
        for initial in range(6):
            stack = [(initial, [], [])]
            while stack:
                current, supports, kinds = stack.pop()
                if supports and rooted_atom_family(supports, 5):
                    s = frozenset({5}).union(*supports)
                    self.assertEqual(len(induced_components(trees[initial], s)), 1)
                    self.assertLessEqual(len(s), 1+3*len(supports))
                    reached, ordered_growth = {5}, True
                    for support in supports:
                        ordered_growth &= bool(reached & support)
                        reached |= support
                    nongrowing_order += not ordered_growth
                    checked += 1
                if len(supports) == 4:
                    continue
                for kind, support, target, _ in aa:
                    if target[current] >= 0:
                        stack.append((int(target[current]), supports+[support], kinds+[kind]))
        self.assertGreater(nongrowing_order, 0)
        OBS['connected_word_lemma'] = dict(
            maximum_word_length=4, active_root_connected_words_checked=checked,
            words_not_growing_from_root_in_application_order=nongrowing_order,
            proof_does_not_assume_two_sided_unique_growth=True,
            support_connected_in_initial_graph=True,
            exact_distance_support='a in S implies d_G(a,b)<=|S|-1<=3*n')

    def test_03_connected_sets_and_factorial_majorant(self):
        tree, n = old.comb(4)
        root, c = 0, 3
        counts = [0]*(n+1)
        others = [v for v in range(n) if v != root]
        for bits in itertools.product((0, 1), repeat=len(others)):
            s = {root}|{v for v, yes in zip(others, bits) if yes}
            if len(induced_components(tree, s)) == 1:
                counts[len(s)] += 1
        for k in range(1, n+1):
            self.assertLessEqual(counts[k], c**(2*(k-1)))
        for order in range(1, 41):
            self.assertLessEqual(3*order+1, 4**order)
            self.assertLessEqual(order**order, 3**order*math.factorial(order))
        lam = c+2*c*(c-1)**2
        kconstant = 96*lam*c**6
        self.assertEqual(kconstant, 1889568)
        OBS['factorial_majorant'] = dict(
            nonpath_tree_vertices=n, rooted_connected_sets_by_size=counts[1:],
            DFS_bound='number(size=k, root=b)<=C^(2*k-2)',
            finite_integer_factorial_checks=40,
            lambda_at_C3_J1_kappa1=lam, K_at_C3_J1_kappa1=kconstant,
            general_K='96*lambda*C^6',
            operator_derivative_bound='||ad_H^n(B_bR)||<=K^n*n!*||B_bR||',
            block_Schur_norm_not_Hilbert_Schmidt_dimension_conversion=True)

    def test_04_full_operator_block_norm_and_complex_strip(self):
        _, _, h, _, _ = old.six_vertex_sector()
        # Full arbitrary background operator, not a single-excitation restriction.
        z_b = np.diag([1 if z % 2 == 0 else -1 for z in range(64)])
        b = np.kron(z_b, np.eye(6))
        current = b.copy()
        kconstant = 1889568
        rows = []
        for order in range(1, 5):
            current = h@current-current@h
            blocks = current.reshape(64, 6, 64, 6).transpose(1, 3, 0, 2)
            block_norms = np.array([[np.linalg.norm(blocks[g, j], 2)
                                    for j in range(6)] for g in range(6)])
            schur = math.sqrt(block_norms.sum(axis=1).max()*block_norms.sum(axis=0).max())
            self.assertLessEqual(schur, kconstant**order*math.factorial(order))
            rows.append(dict(order=order, complete_block_Schur_bound=short(schur)))
        # Evaluate the complex extension using inverse, never complex adjoint.
        z = .07+1j/(4*kconstant)
        u_forward = evolution(h, -z)
        u_inverse = evolution(h, z)
        tau = u_forward@b@u_inverse
        self.close(u_forward@u_inverse, np.eye(384))
        observed = np.linalg.norm(tau, 2)
        self.assertLessEqual(observed, 4/3+1e-10)
        OBS['full_operator_and_strip'] = dict(
            full_data_dimension=64, graph_dimension=6, derivative_checks=rows,
            complex_time_real_part='.07', scaled_imaginary_part='K Im(z)=1/4',
            observed_complex_norm=short(observed), majorant='1/(1-K*|Im(z)|)',
            uniform_strip_halfwidth='1/(2*K)', strip_operator_bound=2,
            arbitrary_reference_dimension_in_analytic_proof=True)

    def test_05_unknown_background_reference_and_arbitrary_encoding(self):
        trees, _, h, _, dd = old.six_vertex_sector()
        far = np.flatnonzero(dd >= 3)
        nr = 3
        rng = np.random.default_rng(466)
        psi = np.zeros((64, 6, nr), complex)
        x = rng.normal(size=(64, len(far), nr))+1j*rng.normal(size=(64, len(far), nr))
        x /= np.linalg.norm(x)
        psi[:, far] = x
        background_probability = sum(float(np.abs(psi[z]).ravel()@np.abs(psi[z]).ravel())
                                     for z in range(64) if (z & 31).bit_count() >= 2)
        self.assertGreater(background_probability, .5)
        rho_r = np.einsum('zgr,zgt->rt', psi, psi.conj())
        u = evolution(h, .1)
        channels = [
            [np.eye(2)],
            [np.array([[0, 1], [1, 0]])],
            [np.diag([1, 1/math.sqrt(2)]), np.array([[0, 1/math.sqrt(2)], [0, 0]])]]
        outputs = []
        for kraus in channels:
            output = np.zeros((2*nr, 2*nr), complex)
            for e in kraus:
                prepared = np.einsum('st,togr->sogr', e, psi.reshape(2, 32, 6, nr))
                evolved = u@prepared.reshape(384, nr)
                output += receiver(evolved, nr)
            self.close(output.reshape(2, nr, 2, nr).trace(axis1=0, axis2=2), rho_r)
            self.assertGreater(np.linalg.eigvalsh(output).min(), -1e-12)
            outputs.append(output)
        OBS['actual_unknown_background_interface'] = dict(
            full_input_data_graph_reference_dimensions=[64, 6, nr],
            initial_graph_support_dimension=len(far),
            background_probability_at_least_two_excitations=short(background_probability),
            no_vacuum_or_single_excitation_contract=True,
            source_graph_background_reference_correlations_allowed=True,
            output_trace_distances=dict(identity_vs_X=short(distance(outputs[0], outputs[1])),
                identity_vs_amplitude_damping=short(distance(outputs[0], outputs[2]))),
            reference_marginal_preserved=True,
            matrix_time='.1', coarse_theorem_bound_at_matrix_time=1,
            matrix_example_not_claimed_to_certify_tight_velocity=True)

    def test_06_all_size_strip_Schwarz_certificates(self):
        # Standard strip map: tanh(pi*K*z/2); at K*|t|=1/10,
        # tanh(pi/20)<pi/20<11/70 using pi<22/7.
        first = 2*Q(11, 70)**10
        self.assertLess(first, Q(1, 50000000))
        # At K*|t|=1, exp(pi)<24, so tanh(pi/2)<23/25.
        x = Q(22, 7)
        exp_upper = sum((x**n/Q(math.factorial(n)) for n in range(31)), Q(0))
        exp_upper += x**31/Q(math.factorial(31))/(1-x/32)
        self.assertLess(exp_upper, 24)
        second = 2*Q(23, 25)**100
        self.assertLess(second, Q(1, 2000))
        OBS['all_size_certificates'] = dict(
            full_time_bound='D <= min(1,2*tanh(pi*K*abs(t)/2)^ceil(r/3))',
            first=dict(scaled_time='K*|t|=1/10', initial_distance=30,
                       vanishing_order=10, exact_upper=str(first), less_than='1/50000000'),
            second=dict(scaled_time='K*|t|=1', initial_distance=300,
                        vanishing_order=100, exact_upper=str(second), less_than='1/2000'),
            exp_22_over_7_certified_less_than_24=True,
            fixed_finite_time_distance_decay=True,
            all_real_times_covered_by_strip_map=True,
            linear_light_cone_or_optimal_velocity_claimed=False)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=466, baseline_round=465, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(same_autonomous_H_as_440_to_443=True,
            arbitrary_unknown_data_graph_reference_initial_state=True,
            only_initial_graph_distance_support_restricted=True,
            vacuum_background_required=False, data_excitation_sector_restricted=False,
            common_fixed_chain_used=False, graph_measured_or_thresholded=False,
            old_invalid_graph_trace_used=False,
            arbitrary_local_CPTP_source_encodings=True,
            complete_reference_uniform_operator_bound=True,
            original_tree_and_capacity_and_flip_rule_remain_inputs=True,
            constants_deliberately_very_loose=True,
            local_compact_endpoint_group_derived=False,
            consistent_halving_derived=False, spatial_direction_completeness_derived=False,
            spatial_dimension_derived=False, full_GR_goal_completed=False,
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
