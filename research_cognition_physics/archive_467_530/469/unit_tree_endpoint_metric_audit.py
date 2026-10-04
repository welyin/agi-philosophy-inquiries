"""Round 469: unit-edge binary-tree means on selected endpoint labels.

Baseline 467. The Euclidean target, selected endpoint labels, distribution and
common scale are preparation inputs, not generated spatial dimension.
"""
import argparse
from collections import deque
from dataclasses import dataclass
from fractions import Fraction as Q
import io
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np

TARGET = Path(__file__).with_name('unit_tree_endpoint_metric_audit_results.json')
OBS = {}


def short(x):
    return float(f'{float(x):.12g}')


@dataclass(frozen=True)
class Radical:
    """Exact linear combinations a+b*sqrt(2)+c*sqrt(3); no field division."""
    a: Q = Q(0)
    b: Q = Q(0)
    c: Q = Q(0)

    def __add__(self, other):
        if not isinstance(other, Radical):
            other = Radical(Q(other))
        return Radical(self.a+other.a, self.b+other.b, self.c+other.c)

    __radd__ = __add__

    def __neg__(self):
        return Radical(-self.a, -self.b, -self.c)

    def __sub__(self, other):
        return self + (-other if isinstance(other, Radical) else -Q(other))

    def __mul__(self, scalar):
        scalar = Q(scalar)
        return Radical(self.a*scalar, self.b*scalar, self.c*scalar)

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        return self*Q(1, scalar)

    def bounds(self):
        lo2, hi2 = Q(1414213562373095, 10**15), Q(1414213562373096, 10**15)
        lo3, hi3 = Q(1732050807568877, 10**15), Q(1732050807568878, 10**15)
        lo = self.a+self.b*(lo2 if self.b >= 0 else hi2)+self.c*(lo3 if self.c >= 0 else hi3)
        hi = self.a+self.b*(hi2 if self.b >= 0 else lo2)+self.c*(hi3 if self.c >= 0 else lo3)
        return lo, hi

    def __float__(self):
        return float(self.a)+float(self.b)*2**.5+float(self.c)*3**.5

    def exact(self):
        return dict(one=str(self.a), sqrt2=str(self.b), sqrt3=str(self.c))


ZERO, ONE = Radical(), Radical(Q(1))
R2, R3 = Radical(b=Q(1)), Radical(c=Q(1))


def edge(a, b):
    return (min(a, b), max(a, b))


def adjacent(edges, total):
    out = [set() for _ in range(total)]
    for a, b in edges:
        out[a].add(b)
        out[b].add(a)
    return out


def tree_data(n, shore, length):
    """Same labels and degree for every shore at fixed n,length."""
    assert n >= 2 and length >= 1 and shore and len(shore) < n
    shore = set(shore)
    nxt, edges = n, set()

    def rooted(leaves):
        nonlocal nxt
        if len(leaves) == 1:
            return leaves[0]
        root = nxt
        nxt += 1
        edges.add(edge(root, leaves[0]))
        edges.add(edge(root, rooted(leaves[1:])))
        return root

    a, b = rooted(sorted(shore)), rooted(sorted(set(range(n))-shore))
    edges.add(edge(a, b))
    assert nxt == 2*n-2
    base = frozenset(edges)
    edges.remove(edge(a, b))
    chain = [a]+list(range(2*n-2, 2*n+length-3))+[b]
    for x, y in zip(chain, chain[1:]):
        edges.add(edge(x, y))
    for k, x in enumerate(chain[1:-1]):
        edges.add(edge(x, 2*n+length-3+k))
    total = 2*n+2*length-4
    degrees = [1]*n+[3]*(n+length-3)+[1]*(length-1)
    assert len(degrees) == total
    return frozenset(edges), total, degrees, base, (a, b)


def distances(edges, total, n):
    adj = adjacent(edges, total)
    result = np.zeros((n, n), dtype=np.int64)
    for source in range(n):
        ds = [-1]*total
        ds[source] = 0
        todo = deque([source])
        while todo:
            x = todo.popleft()
            for y in adj[x]:
                if ds[y] < 0:
                    ds[y] = ds[x]+1
                    todo.append(y)
        assert all(x >= 0 for x in ds)
        result[source] = ds[:n]
    return result


def cut(shore, n):
    return np.array([[int((a in shore) != (b in shore)) for b in range(n)]
                     for a in range(n)], dtype=np.int64)


def cube_terms():
    weights = {1:(R2+R3-ONE)/4,
               2:(R2-R3+ONE)/4,
               3:(R3-3*R2+3*ONE)/4}
    return [(set(x for x in range(8) if (x&mask).bit_count()%2),
             weights[mask.bit_count()]) for mask in range(1, 8)]


def nni(edges, a, b, c, d):
    removed = {edge(a, b), edge(c, d)}
    added = {edge(a, c), edge(b, d)}
    assert len({a,b,c,d}) == 4 and edge(b,c) in edges
    assert removed <= edges and not (added & edges)
    return frozenset((edges-removed)|added)


class Audit(unittest.TestCase):
    def check_tree(self, edges, total, degrees):
        self.assertEqual(len(edges), total-1)
        adj = adjacent(edges, total)
        self.assertEqual([len(a) for a in adj], degrees)
        reached, todo = set(), [0]
        while todo:
            x = todo.pop()
            if x not in reached:
                reached.add(x)
                todo.extend(adj[x]-reached)
        self.assertEqual(len(reached), total)

    def test_01_exact_cube_cut_certificate(self):
        for rad, lo, hi in [(2, Q(1414213562373095,10**15), Q(1414213562373096,10**15)),
                            (3, Q(1732050807568877,10**15), Q(1732050807568878,10**15))]:
            self.assertLess(lo*lo, rad)
            self.assertGreater(hi*hi, rad)
        terms = cube_terms()
        total = sum((w for _, w in terms), ZERO)
        self.assertEqual(total, (3*ONE+3*R2+R3)/4)
        for _, w in terms:
            self.assertGreater(w.bounds()[0], 0)
        for a, b in itertools.combinations(range(8), 2):
            value = sum((w for s,w in terms if (a in s)!=(b in s)), ZERO)
            self.assertEqual(value, {1:ONE, 2:R2, 3:R3}[(a^b).bit_count()])
        OBS['cut_representation'] = dict(
            cube_points=8, target_distances=['1','sqrt(2)','sqrt(3)'],
            parity_cuts=7, all_coefficients_certified_positive=True,
            coefficients_by_mask_size={str(k):cube_terms()[[0,2,6][k-1]][1].exact() for k in [1,2,3]},
            total_weight=total.exact(), total_weight_approx=short(float(total)),
            all_28_pair_distances_exact=True,
            Euclidean_L1_cut_representation_is_mature_mathematics=True)

    def test_02_same_label_degree_unit_tree_construction(self):
        checked = 0
        for n in range(2, 7):
            for mask in range(1, 2**(n-1)):
                s = {i+1 for i in range(n-1) if mask>>i&1}
                for length in (1, 7):
                    edges, total, degrees, base, special = tree_data(n, s, length)
                    self.check_tree(edges, total, degrees)
                    self.assertEqual(len(base), 2*n-3)
                    d = distances(edges, total, n)
                    residual = d-length*cut(s,n)
                    self.assertGreaterEqual(residual.min(), 0)
                    self.assertLessEqual(residual.max(), 2*n-4)
                    self.assertTrue(set(degrees) <= {1,3})
                    checked += 1
        OBS['unit_same_budget_construction'] = dict(
            exact_constructed_trees=checked,
            observed_endpoint_counts=[2,3,4,5,6], central_path_lengths=[1,7],
            total_vertices='N=2*n+2*L-4',
            added_unobserved_leaf_labels='L-1', added_internal_labels='L-1',
            same_per_label_degrees_for_every_cut=True,
            degree_two_vertices=0, every_edge_unit=True,
            exact_distance_identity='d_T(i,j)=L*delta_S(i,j)+r_S(i,j)',
            residual_range='0 <= r_S(i,j) <= 2*n-4')

    def test_03_large_cube_resource_and_error(self):
        n, length, terms = 8, 4096, cube_terms()
        total_weight = sum((w for _,w in terms), ZERO)
        max_error, matrices, graphs = ZERO, [], []
        degrees_first = None
        for s, w in terms:
            edges, total, degrees, _, _ = tree_data(n,s,length)
            self.check_tree(edges,total,degrees)
            if degrees_first is None:
                degrees_first = degrees
            self.assertEqual(degrees,degrees_first)
            matrices.append(distances(edges,total,n))
            graphs.append(edges)
        self.assertEqual(len(set(graphs)), 7)
        bound = total_weight*(2*n-4)/length
        exact_errors = []
        for a,b in itertools.combinations(range(n),2):
            target = {1:ONE,2:R2,3:R3}[(a^b).bit_count()]
            mean = sum((w*int(mat[a,b])/length for (_,w),mat in zip(terms,matrices)),ZERO)
            error = mean-target
            self.assertGreaterEqual(error.bounds()[0], 0)
            self.assertGreaterEqual((bound-error).bounds()[0], 0)
            if float(error)>float(max_error):
                max_error=error
            exact_errors.append(error)
        # Floating values only nominate a candidate; certify every exact ordering.
        for error in exact_errors:
            self.assertGreaterEqual((max_error-error).bounds()[0], 0)
        self.assertLess(bound.bounds()[1], Q(1,100))
        OBS['cube_unit_sector'] = dict(
            N=total, observed_leaves=8, unobserved_leaves=length-1,
            internal_degree_three=8+length-3, central_length=length,
            graph_support_count=7, common_scale='W/4096',
            common_scale_approx=short(float(total_weight)/length),
            prepared_probabilities=[short(float(w)/float(total_weight)) for _,w in terms],
            maximum_absolute_distance_error=short(float(max_error)),
            maximum_error_exact=max_error.exact(),
            theorem_upper_bound=short(float(bound)), theorem_upper_bound_exact=bound.exact(),
            theorem_upper_bound_less_than='1/100',
            common_N_common_labels_common_degrees=True,
            no_weighted_edge_or_degree_two_shortcut=True,
            Hamiltonian_state_preparation_not_implemented=True)

    def test_04_NNI_sector_and_internal_label_swap(self):
        checked_flips, checked_swaps = 0, 0
        for s,_ in cube_terms():
            edges,total,degrees,_,_ = tree_data(8,s,4)
            adj=adjacent(edges,total)
            for u,v in edges:
                if len(adj[u]) != 3 or len(adj[v]) != 3:
                    continue
                left,right=sorted(adj[u]-{v}),sorted(adj[v]-{u})
                for a,c in itertools.product(left,right):
                    changed=nni(edges,a,u,v,c)
                    self.check_tree(changed,total,degrees)
                    checked_flips+=1
                first=nni(edges,left[0],u,v,right[0])
                second=nni(first,left[1],u,v,right[1])
                def swap(x):
                    return v if x==u else u if x==v else x
                relabeled=frozenset(edge(swap(a),swap(b)) for a,b in edges)
                self.assertEqual(second,relabeled)
                checked_swaps+=1
        OBS['fixed_degree_NNI_interface'] = dict(
            active_NNI_moves_checked=checked_flips,
            two_NNI_internal_label_transpositions_checked=checked_swaps,
            all_other_labels_fixed=True, entire_fixed_degree_tree_sector_closed=True,
            leaf_labeled_NNI_connectivity_attributed_to_existing_theory=True,
            connectivity_not_autonomous_distribution_preparation=True)

    def test_05_finite_budget_leaf_obstruction(self):
        n,length=3,64
        terms=[({1,2},Q(1)),({2},Q(1))]
        mats=[]
        for s,_ in terms:
            edges,total,degrees,_,_=tree_data(n,s,length)
            self.check_tree(edges,total,degrees)
            mats.append(distances(edges,total,n))
        mean=[[sum(w*Q(int(mat[a,b]),length) for (_,w),mat in zip(terms,mats))
               for b in range(n)] for a in range(n)]
        error=max(abs(mean[a][b]-abs(a-b)) for a in range(n) for b in range(n))
        lower=Q(8,3*total+4)
        self.assertGreaterEqual(error,lower)
        self.assertEqual(error,Q(3,64))
        scale=Q(2,length)
        self.assertGreaterEqual(mean[0][1]+mean[1][2]-mean[0][2],2*scale)
        for mat in mats:
            self.assertLessEqual(int(mat.max()),total//2)
        OBS['finite_budget_limit'] = dict(
            universal_leaf_triangle_slack='mean(a,b)+mean(b,c)-mean(a,c)>=2*s',
            universal_diameter_upper='mean(a,c)<=s*N/2',
            collinear_target_distances=['ell','ell','2*ell'],
            universal_max_error_lower='epsilon>=8*ell/(3*N+4)',
            no_claim_of_sharp_lower_bound=True,
            example_N=total, example_max_error=str(error),
            example_certified_lower=str(lower), exact_finite_collinear_representation_impossible=True)

    def test_06_graph_coherence_and_reference(self):
        terms,length=cube_terms(),32
        weights=np.array([float(w) for _,w in terms])
        probability=weights/weights.sum()
        rng=np.random.default_rng(469)
        refs=rng.normal(size=(7,3))+1j*rng.normal(size=(7,3))
        refs/=np.linalg.norm(refs,axis=1)[:,None]
        phase=np.exp(1j*np.arange(7)*.37)
        psi=np.sqrt(probability)[:,None]*phase[:,None]*refs
        rho_g=psi@psi.conj().T
        self.assertLess(np.linalg.norm(rho_g.diagonal().real-probability),1e-13)
        off=rho_g-np.diag(np.diag(rho_g))
        self.assertGreater(np.linalg.norm(off),.1)
        max_difference=0.
        for a,b in itertools.combinations(range(8),2):
            values=[]
            for s,_ in terms:
                edges,total,_,_,_=tree_data(8,s,length)
                values.append(float(distances(edges,total,8)[a,b])*weights.sum()/length)
            diagonal=np.diag(values)
            mean_joint=np.vdot(psi,diagonal@psi).real
            mean_prepared=probability@values
            max_difference=max(max_difference,abs(mean_joint-mean_prepared))
        self.assertLess(max_difference,1e-13)
        OBS['coherent_reference_interface'] = dict(
            graph_support_dimension=7, reference_dimension=3,
            graph_off_diagonal_Frobenius_norm=short(np.linalg.norm(off)),
            maximum_diagonal_mean_residual=short(max_difference),
            graph_reference_state_not_dephased=True,
            general_guarantee_uses_graph_diagonal_operator_identity=True,
            arbitrary_data_and_reference_allowed_at_fixed_graph_populations=True,
            actual_distance_measurement_device_or_time_average_not_derived=True)


def run():
    OBS.clear()
    stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=469,baseline_round=467,tests_run=result.testsRun,
        failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,observations=OBS,
        scope=dict(selected_endpoint_mean_representation=True,
            arbitrary_finite_Euclidean_dimension=True,
            same_unit_edge_fixed_per_label_degree_sector=True,
            auxiliary_labels_explicitly_counted=True,
            arbitrary_dimension_not_selected_as_three=True,
            source_of_endpoint_selection_not_derived=True,
            graph_distribution_and_scale_prepared_inputs=True,
            same_H_autonomous_preparation_not_proved=True,
            stable_geometry_or_ergodic_time_average_not_proved=True,
            actual_distance_readout_not_implemented=True,
            physical_spatial_dimension_derived=False,
            full_GR_goal_completed=False,phase_closure_triggered=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.check:
        assert result==json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
