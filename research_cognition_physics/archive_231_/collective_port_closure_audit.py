"""Round 505: common coherent port spaces of the unaveraged tree family.

The exact contract keeps every graph and passive reference. The approximate
obstruction concerns the declared uniform leaf mode and two-role compression.
No spatial dimension is inferred from a Hilbert-space dimension.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import coherent_graph_mean_obstruction as old

HERE = Path(__file__).resolve().parent
TARGET = HERE/'collective_port_closure_audit_results.json'
OBS = {}


def laplacian(tree, n):
    a = np.zeros((n, n), dtype=np.int64)
    for u, v in tree:
        a[u, v] = a[v, u] = 1
    return np.diag(a.sum(axis=1))-a


def path_tree(i, p, q, a, b):
    """A legal tree with p--q, p--a and q--b, for any distinct role labels."""
    n = 2*i+2
    order = [p, q]+[v for v in range(i) if v not in (p, q)]
    tree = {old.tree_tools.edge(u, v) for u, v in zip(order, order[1:])}
    tree.update({old.tree_tools.edge(p, a), old.tree_tools.edge(q, b)})
    degree = [0]*n
    for u, v in tree:
        degree[u] += 1
        degree[v] += 1
    leaves = iter(v for v in range(i, n) if v not in (a, b))
    for v in range(i):
        for _ in range(3-degree[v]):
            tree.add(old.tree_tools.edge(v, next(leaves)))
    return frozenset(tree)


def modular_constraint_rank(matrices, prime=1000003):
    """Exact rank over F_p; combined with explicit rational kernel witnesses."""
    n = matrices[0].shape[0]
    eye = np.eye(n, dtype=np.int64)
    pivots = {}
    for a in matrices:
        block = np.kron(eye, a)-np.kron(a.T, eye)
        for row in block:
            row = row.copy() % prime
            for col, pivot in sorted(pivots.items()):
                if row[col]:
                    row = (row-row[col]*pivot) % prime
            nonzero = np.flatnonzero(row)
            if len(nonzero):
                col = int(nonzero[0])
                pivots[col] = row*pow(int(row[col]), -1, prime) % prime
    return len(pivots)


def moments(i):
    b = Q(3*i, 2*i-1)
    c = Q(i+2, i)
    return b, c, b-c


def certificate():
    t = Q(1, 4096)
    factor = Q(2, 15)-192*t
    assert factor == Q(83, 960) and factor > Q(1, 16)
    return t, factor, t*t/16


class Audit(unittest.TestCase):
    def test_01_nni_cross_operators_and_legal_witnesses(self):
        count = 0
        for i in (2, 3, 4, 7):
            n = 2*i+2
            known = set(old.previous.ensemble(i)[0]) if i <= 4 else None
            for p, q in ((0, 1), (i-1, 0)):
                for a, b in ((i, i+1), (n-1, n-2)):
                    tree = path_tree(i, p, q, a, b)
                    other = set(tree)-{old.tree_tools.edge(p, a), old.tree_tools.edge(q, b)}
                    other.update({old.tree_tools.edge(p, b), old.tree_tools.edge(q, a)})
                    other = frozenset(other)
                    lg, lh = laplacian(tree, n), laplacian(other, n)
                    self.assertTrue(np.array_equal(lg.diagonal(), [3]*i+[1]*(i+2)))
                    self.assertIn(other, old.tree_tools.tree_flips(tree, n))
                    if known is not None:
                        self.assertIn(tree, known)
                        self.assertIn(other, known)
                    x, y = np.eye(n, dtype=np.int64)[:, p]-np.eye(n, dtype=np.int64)[:, q], np.eye(n, dtype=np.int64)[:, a]-np.eye(n, dtype=np.int64)[:, b]
                    self.assertTrue(np.array_equal(lh-lg, np.outer(x, y)+np.outer(y, x)))
                    count += 1
        OBS['nni_witnesses'] = dict(exact_integer_identities=count, I_values=[2, 3, 4, 7])

    def test_02_common_commutants_exact_modular_certificates(self):
        rows = []
        for i, kernel_dimension in ((2, 3), (3, 2)):
            n = 2*i+2
            trees = old.previous.ensemble(i)[0]
            matrices = [laplacian(tree, n) for tree in trees]
            witnesses = [np.eye(n, dtype=np.int64), np.ones((n, n), dtype=np.int64)]
            if i == 2:
                r = np.array([2]*i+[-1]*(i+2), dtype=np.int64)
                witnesses.append(np.outer(r, r))
            for lg in matrices:
                for witness in witnesses:
                    self.assertTrue(np.array_equal(lg @ witness, witness @ lg))
            rank = modular_constraint_rank(matrices)
            self.assertEqual(rank, n*n-kernel_dimension)
            # The explicit independent rational witnesses give the reverse bound.
            entries = [(0, 0), (0, 1)]+([(0, i)] if i == 2 else [])
            minor = [[int(w[a, b]) for a, b in entries] for w in witnesses]
            if i == 2:
                determinant = int(np.dot(minor[0], np.cross(minor[1], minor[2])))
            else:
                determinant = minor[0][0]*minor[1][1]-minor[0][1]*minor[1][0]
            self.assertNotEqual(determinant, 0)
            rows.append(dict(I=i, trees=len(trees), prime=1000003,
                modular_constraint_rank=rank, exact_characteristic_zero_commutant_dimension=kernel_dimension,
                rational_kernel_witnesses=['identity', 'all_ones']+(['role_outer_product'] if i == 2 else [])))
        OBS['finite_common_commutants'] = rows

    def test_03_matrix_units_and_role_coupling(self):
        records = []
        for i in (2, 3, 4):
            du, dv = i-1, i+1
            d = du+dv
            basis = np.eye(d, dtype=np.int64)
            def cross(j, a):
                return np.outer(basis[j], basis[du+a])+np.outer(basis[du+a], basis[j])
            for a in range(dv):
                b = (a+1) % dv
                eab = cross(0, a) @ cross(0, b)
                paa = eab @ eab.T
                self.assertTrue(np.array_equal(paa, np.outer(basis[du+a], basis[du+a])))
                for j in range(du):
                    eja = cross(j, a) @ paa
                    self.assertTrue(np.array_equal(eja, np.outer(basis[j], basis[du+a])))
            n = 2*i+2
            lg = laplacian(path_tree(i, 0, 1, i, i+1), n)
            # Rational role difference: internal entries L, leaf entries -I.
            role = np.array([i+2]*i+[-i]*(i+2), dtype=np.int64)
            image = lg @ role
            internal_deviation = i*image[:i]-image[:i].sum()
            if i == 2:
                self.assertTrue(np.array_equal(image, 3*role))
            else:
                self.assertGreater(int(internal_deviation @ internal_deviation), 0)
            records.append(dict(I=i, cross_matrix_units_exact=True,
                role_couples_to_zero_sum_space=bool(i > 2)))
        OBS['algebra_generators'] = records

    def test_04_full_ensemble_collective_moments(self):
        records = []
        for i in (2, 3, 4):
            n, leaves = 2*i+2, i+2
            trees = old.previous.ensemble(i)[0]
            g = len(trees)
            moment, double_parent = 0, 0
            z = np.zeros((n, 2), dtype=np.int64)
            z[:i, 0] = 1
            z[i:, 1] = 1
            target = leaves*np.array([[1, -1], [-1, 1]], dtype=np.int64)
            for tree in trees:
                adj = old.tree_tools.adjacency(tree, n)
                m = [sum(v >= i for v in adj[u]) for u in range(i)]
                self.assertEqual(sum(m), leaves)
                moment += sum(x*x for x in m)
                double_parent += int(i in adj[0] and i+1 in adj[0])
                self.assertTrue(np.array_equal(z.T @ laplacian(tree, n) @ z, target))
            b, c, a = moments(i)
            self.assertEqual(Q(moment, g*leaves), b)
            self.assertEqual(Q(double_parent, g), Q(1, i*(2*i-1)))
            self.assertEqual(a, Q((i-1)*(i-2), i*(2*i-1)))
            records.append(dict(I=i, trees=g, first_order_internal_norm_squared=str(b),
                compressed_coefficient=str(c), leakage_coefficient=str(a),
                specified_double_parent_probability=str(Q(double_parent, g)), all_graph_compressions_equal=True))
        OBS['ensemble_moments'] = records

    def test_05_uniform_time_certificate_and_actual_port_gap(self):
        time_max, factor, floor = certificate()
        for i in (3, 4, 10, 1000000):
            b, c, a = moments(i)
            self.assertEqual(a-Q(2, 15), Q((i-3)*(11*i-10), 15*i*(2*i-1)))
            self.assertGreaterEqual(a, Q(2, 15))
            self.assertLess(b, Q(9, 4))
            self.assertLess(c, Q(16, 9))
        i = 3
        trees, _, _, _, h, _, _, _ = old.model(i)
        n, g, leaves = 8, len(trees), 5
        source = np.zeros(n*g)
        source[i*g:] = 1/np.sqrt(leaves*g)
        rows = []
        b, c, _ = moments(i)
        for t in (time_max, time_max/2):
            final = old.evolve(h, source, t)
            probability = float(np.linalg.norm(final[:i*g])**2)
            compressed = 4*float(c)/(1+float(c))**2*math.sin((1+float(c))*float(t)/2)**2
            x = (4*(i-1)+6)*t
            tail = x**19/(math.factorial(19)*(1-x/20))
            self.assertGreater(probability-compressed-5*float(tail), float(t*t/16))
            rows.append(dict(I=i, time=str(t), actual_internal_probability=old.short(probability),
                two_role_probability=old.short(compressed), probability_gap=old.short(probability-compressed),
                rigorous_gap_floor=str(t*t/16), numerical_tail_bound=str(tail)))
        OBS['probability_certificate'] = dict(time_max=str(time_max), rigorous_factor=str(factor),
            conservative_gap_at_time_max=str(floor), all_I_at_least_3=True, diagnostics=rows)

    def test_06_exact_small_case_and_graph_reference_interface(self):
        trees, f, _, _, h, _, _, _ = old.model(2)
        n, g = 6, len(trees)
        z = np.zeros((n, 2), dtype=np.int64)
        z[:2, 0], z[2:, 1] = 1, 1
        embedding = np.kron(z, np.eye(g, dtype=np.int64))
        m = np.array([[2, -2], [-1, 1]], dtype=np.int64)
        effective = np.kron(np.eye(2, dtype=np.int64), f)-np.kron(m, np.eye(g, dtype=np.int64))
        self.assertTrue(np.array_equal(h @ embedding, embedding @ effective))
        # Uniform graph populations, arbitrary phase, and an entangled passive R.
        i = 3
        trees, _, _, _, h, _, _, _ = old.model(i)
        n, g, leaves = 8, len(trees), 5
        source = np.zeros((n*g, g), dtype=complex)
        graph_reference = np.diag(np.exp(2j*np.pi*np.arange(g)/g))/np.sqrt(g)
        for v in range(i, n):
            source[v*g:(v+1)*g] = graph_reference/np.sqrt(leaves)
        t = Q(1, 4096)
        final = old.evolve(h, source, t)
        reference_residual = np.linalg.norm(final.conj().T @ final-source.conj().T @ source)
        probabilities = [float(np.linalg.norm(final[v*g:(v+1)*g])**2) for v in range(n)]
        self.assertLess(reference_residual, 1e-12)
        self.assertAlmostEqual(sum(probabilities), 1, places=12)
        c = float(moments(i)[1])
        compressed = 4*c/(1+c)**2*math.sin((1+c)*float(t)/2)**2
        self.assertGreater(sum(probabilities[:i])-compressed, float(t*t/16))
        OBS['interface'] = dict(I2_full_unknown_graph_intertwining_exact=True,
            graph_reference_retained=True, passive_reference_residual=old.short(reference_residual),
            all_port_outcomes_retained=True, actual_internal_probability=old.short(sum(probabilities[:i])),
            reference_dimension=g, graph_measured_or_reset=False)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output, verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(output.getvalue())
    return dict(round=505, scientific_baseline_round=504,
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__,
        dependency_sha256={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'coherent_graph_mean_obstruction.py', 'role_covariant_graph_obstruction.py',
            'branching_tree_distance_audit.py', 'research_note_496.md')},
        scope=dict(common_graph_independent_subspaces_classified=True, all_graphs_retained_in_exact_contract=True,
            I2_exception_retained=True, two_role_fixed_time_uniform_size_gap=True,
            uniform_population_source_can_have_graph_reference_correlations=True,
            collective_source_preparation_derived=False, all_approximate_codes_excluded=False,
            dimension_three_generated=False, full_GR_goal_completed=False, phase_closure_triggered=False),
        observations=OBS)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    data = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(data, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(data, ensure_ascii=False, indent=2))
