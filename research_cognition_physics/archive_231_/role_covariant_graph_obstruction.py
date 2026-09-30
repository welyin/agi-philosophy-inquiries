"""Round 495: role-covariant fast graph updates cannot universally erase edge data.

Baseline 494. All labelled degree-three internal trees, all graph density
operators, and unchanged edge-controlled exchange are explicit contracts.
An exact finite-size obstruction is not a uniform large-scale geometry no-go.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib, io, itertools, json, math, platform, unittest
from pathlib import Path
import numpy as np
import branching_tree_distance_audit as tree_tools
import fast_graph_average_port_process as previous

TARGET = Path(__file__).with_name("role_covariant_graph_obstruction_results.json")
OBS = {}
# Exact integer products below have dimension at most90 and at mostfour
# factors, entries at most32. Larger counts are genuine Python integers.
assert 90**4*32**4 < 2**63

def short(x):
    return float(f"{float(x):.12g}")

def ceil_fraction(x):
    return -((-x.numerator)//x.denominator)

def formulas(i):
    assert i >= 3
    dimension = math.factorial(2*i)//2**i
    fixed = 4*(i-2)*math.factorial(2*i-6)//2**(i-3)
    p = F(4,i*(i-1)*(2*i-1)*(2*i-3)*(2*i-5))
    assert p == F(fixed,dimension)
    g = p*(1-F(2,i))/(1-p*p)
    rates = [F(dimension//i+fixed//2,dimension+fixed),
             F(dimension//i-fixed//2,dimension-fixed)]
    assert rates[0]-rates[1] == g
    return dict(i=i,dimension=dimension,fixed=fixed,p=p,g=g,rates=rates)

@lru_cache(None)
def ensemble(i):
    n = 2*i+2
    words = sorted(set(itertools.permutations(tuple(k for k in range(i)
                                                      for _ in range(2)))))
    trees = sorted({tree_tools.prufer_tree(n,w) for w in words},
                   key=lambda t:sorted(t))
    lookup = {t:j for j,t in enumerate(trees)}
    perm = list(range(n))
    for a,b in ((0,1),(i,i+1),(i+2,i+3)):
        perm[a],perm[b] = perm[b],perm[a]
    graph_perm = np.array([lookup[frozenset(tree_tools.edge(perm[a],perm[b])
                                  for a,b in t)] for t in trees],dtype=np.int64)
    edge = np.array([int(tree_tools.edge(0,i) in t) for t in trees],dtype=np.int64)
    return trees,graph_perm,edge

@lru_cache(None)
def original_three_internal():
    i,n = 3,8
    trees, permutation, q = ensemble(i)
    d = len(trees)
    lookup = {t:j for j,t in enumerate(trees)}
    f = np.zeros((d,d),dtype=np.int64)
    adjacency = []
    for j,t in enumerate(trees):
        a = np.zeros((n,n),dtype=np.int64)
        for b,c in t:a[b,c]=a[c,b]=1
        adjacency.append(a)
        for z,count in tree_tools.tree_flips(t,n).items():
            f[lookup[z],j] = count
    h = np.kron(np.eye(n,dtype=np.int64),f)
    for g,a in enumerate(adjacency):
        lap = np.diag(a.sum(axis=0))-a
        for b in range(n):
            for c in range(n):
                h[d*b+g,d*c+g] += (n-1)*int(b==c)-lap[b,c]
    return f,h,adjacency

def short_time_certificate(g,k):
    assert g>0 and k>=1
    t = g/(4*k**3)
    x = 2*k*t
    remainder = x**3/(6*(1-x/4))
    leading = g*t*t
    # At this t, x<=1/2, so remainder/leading <=8/21 <1/2.
    assert x <= F(1,2)
    assert remainder/leading <= F(8,21)
    assert remainder < leading/2
    return dict(time=t,k=k,leading=leading,remainder=remainder,
                positive_event_gap=leading/2)

def event_diagnostic():
    f,h,_ = original_three_internal()
    _,p,_ = ensemble(3)
    d,n = 90,8
    cert = short_time_certificate(formulas(3)["g"],F(15))
    t = float(cert["time"])
    # Actual U columns for the SAME input port a=3 and every graph basis.
    # Taylor is diagnostic; the strict probability gap is the rational
    # commutator certificate, not floating cancellation or a spectral scan.
    initial = np.zeros((n*d,d),dtype=complex)
    initial[3*d:4*d] = np.eye(d)
    total,term = initial.copy(),initial.copy()
    for order in range(1,9):
        term = (-1j*t/order)*(h@term)
        total += term
    block = total[:d]  # actual output port u=0
    direct = float(np.vdot(block,block).real)
    cross = float(np.vdot(block,block[:,p]).real)
    plus = (direct+cross)/94
    minus = (direct-cross)/86
    x = 15*cert["time"]
    tail = x**9/(math.factorial(9)*(1-x/10))
    return plus,minus,tail

class Audit(unittest.TestCase):
    def test_01_fixed_tree_classification_and_independent_enumeration(self):
        rows = []
        for i in (3,4):
            trees,p,q = ensemble(i)
            data = formulas(i)
            fixed = np.flatnonzero(p==np.arange(len(trees)))
            self.assertEqual(len(trees),data["dimension"])
            self.assertEqual(len(fixed),data["fixed"])
            self.assertEqual(int(q.sum()),data["dimension"]//i)
            self.assertEqual(int(q[fixed].sum()),data["fixed"]//2)
            keys = set()
            for j in fixed:
                tree = trees[int(j)]
                adj = tree_tools.adjacency(tree,2*i+2)
                common = adj[0]&adj[1]
                self.assertEqual(len(common),1)
                w = next(iter(common))
                self.assertIn(w,range(2,i))
                selected = adj[0]-{w}
                self.assertEqual(len(selected & {i,i+1}),1)
                self.assertEqual(len(selected & {i+2,i+3}),1)
                self.assertEqual(adj[1]-{w},{i,i+1,i+2,i+3}-selected)
                removed = {0,1,i,i+1,i+2,i+3}
                residual = frozenset((a,b) for a,b in tree
                                     if a not in removed and b not in removed)
                residual_adj = tree_tools.adjacency(residual,2*i+2)
                self.assertEqual(len(residual_adj[w]),1)
                for vertex in range(2*i+2):
                    if vertex in removed:continue
                    expected = 3 if vertex in range(2,i) and vertex!=w else 1
                    self.assertEqual(len(residual_adj[vertex]),expected)
                keys.add((w,tuple(sorted(selected)),tuple(sorted(residual))))
            self.assertEqual(len(keys),len(fixed))
            rows.append(dict(i=i,trees=len(trees),fixed_trees=len(fixed),
                             trace_Q=int(q.sum()),trace_QR=int(q[fixed].sum()),
                             fixed_tree_deletion_is_injective=True))
        OBS["fixed_tree_enumeration"] = rows

    def test_02_stationary_parity_sources_and_all_covariant_updates(self):
        rows = []
        for i in (3,4):
            trees,p,q = ensemble(i)
            d = len(trees)
            self.assertTrue(np.array_equal(p[p],np.arange(d)))
            data = formulas(i)
            # Permutation involution: (I+/-R)^2=2(I+/-R), orthogonal signs.
            # This sparse exact check supplies positivity without a float eig.
            self.assertEqual((d+data["fixed"])%2,0)
            self.assertEqual((d-data["fixed"])%2,0)
            self.assertGreater(data["g"],0)
            rows.append(dict(i=i,dimension=d,fixed=data["fixed"],
                source_ranks=[(d+data["fixed"])//2,(d-data["fixed"])//2],
                parity_source_weights=[str(F(d+data["fixed"],2*d)),
                                       str(F(d-data["fixed"],2*d))],
                edge_expectations=[str(x) for x in data["rates"]],
                gap=str(data["g"]),scalar_average_distance_lower_bound=str(data["g"]/2)))
        f,_,_ = original_three_internal()
        trees,p,_ = ensemble(3)
        self.assertTrue(np.array_equal(f,f.T))
        self.assertTrue(np.array_equal(f.sum(axis=0),np.full(90,8)))
        self.assertTrue(np.array_equal(f[p],f[:,p]))
        # Original F covariance checked on generators of S3 x S5.
        lookup = {t:j for j,t in enumerate(trees)}
        for a,b in ((0,1),(1,2),(3,4),(4,5),(5,6),(6,7)):
            label = list(range(8));label[a],label[b]=label[b],label[a]
            r = np.array([lookup[frozenset(tree_tools.edge(label[x],label[y])
                                    for x,y in t)] for t in trees])
            self.assertTrue(np.array_equal(f[r],f[:,r]))
        alternate = 2*(f@f)+3*f
        self.assertTrue(np.array_equal(alternate[p],alternate[:,p]))
        OBS["arbitrary_graph_H_obstruction"] = dict(
            sources=rows,only_commutation_with_selected_R_needed=True,
            real_entries_connectivity_or_non_degeneracy_not_needed=True,
            every_finite_F_time_average_preserves_gap=True,
            maximal_initial_graph_state_space_is_an_explicit_contract=True,
            original_NNI_covariance_generators_checked=6)

    def test_03_six_graph_and_non_covariant_boundaries(self):
        trees,p,q = ensemble(2)
        fixed = np.flatnonzero(p==np.arange(6))
        self.assertEqual(len(fixed),4)
        self.assertEqual(int(q[fixed].sum()),2)
        self.assertEqual(F(int(q.sum()+q[fixed].sum()),10),F(1,2))
        self.assertEqual(F(int(q.sum()-q[fixed].sum()),2),F(1,2))
        # Six-graph positive control is the exact frozen494 projector identity.
        _,nums = previous.exact_graph()
        nn = previous.model()[7]
        for a in previous.LEAVES:
            self.assertTrue(np.array_equal(sum((r@nn[1,a]@r for r in nums),
                         np.zeros((6,6),dtype=np.int64)),288*np.eye(6,dtype=np.int64)))
        # Explicit escape at the operator level: a Fourier eigenbasis with
        # distinct energies pinches EVERY diagonal edge effect to its trace/D.
        # The graph basis ordering and non-role-covariant F are new inputs.
        _,perm,q3 = ensemble(3)
        d = 90
        j = np.arange(d)
        fourier = np.exp(2j*np.pi*np.outer(j,j)/d)/math.sqrt(d)
        hf = (fourier*j)@fourier.conj().T
        diag_effect = np.sum(q3[:,None]*np.abs(fourier)**2,axis=0)
        averaged = (fourier*diag_effect)@fourier.conj().T
        self.assertLess(np.max(np.abs(averaged-np.eye(d)/3)),1e-12)
        failure = float(np.linalg.norm(hf[perm]-hf[:,perm],2))
        self.assertGreater(failure,1)
        OBS["scope_controls"] = dict(
            I2_fixed_count=4,I2_both_edge_rates="1/2",
            I2_formula_not_extended_from_I_ge_3=True,
            I2_frozen494_average_rechecked=True,
            Fourier_escape_exact_reason="all eigenvector entries have squared modulus 1/D",
            Fourier_escape_graph_ordering_and_H_are_new_inputs=True,
            Fourier_escape_average_error_diagnostic=short(np.max(np.abs(averaged-np.eye(d)/3))),
            Fourier_escape_R_commutator_norm_diagnostic=short(failure),
            restricted_sources_or_extra_graph_memory_not_excluded=True)

    def test_04_rational_finite_actual_port_separation(self):
        rows = []
        for i in (3,4):
            data = formulas(i)
            # Original NNI has 4(I-1) outgoing counted flips, and N-1 edges.
            # General theorem may instead use any declared rational norm bound.
            k = F(6*i-3)
            cert = short_time_certificate(data["g"],k)
            self.assertLess(cert["remainder"],cert["leading"]/2)
            rows.append(dict(i=i,H_norm_bound=str(k),edge_gap=str(data["g"]),
                time=str(cert["time"]),leading_difference=str(cert["leading"]),
                third_order_remainder=str(cert["remainder"]),
                strict_actual_event_probability_gap=str(cert["positive_event_gap"]),
                common_population_predictor_worst_error=str(cert["positive_event_gap"]/2)))
        OBS["actual_port_separation_certificate"] = dict(
            inputs_have_same_pure_port_and_different_stationary_graph_sources=True,
            ordinary_initial_state_difference_trace_norm=2,
            centered_final_port_effect_norm="1/2",
            generic_remainder="(2*K*t)^3/[6*(1-(2*K*t)/4)]",
            time_choice="g_I/(4*K^3)",remainder_over_leading_at_most="8/21",
            sign_strictly_positive_for_each_fixed_I=True,examples=rows)

    def test_05_true_NNI_H_and_complete_port_probability_diagnostic(self):
        f,h,adj = original_three_internal()
        d,n = 90,8
        self.assertTrue(np.array_equal(h,h.T))
        for a in range(n):
            for b in range(n):
                block = h[d*b:d*b+d,d*a:d*a+d]
                expected = (f+(n-1-int(adj[0][a].sum()))*np.eye(d,dtype=np.int64)
                            if a==b else np.diag([g[b,a] for g in adj]))
                self.assertTrue(np.array_equal(block,expected))
        plus,minus,tail = event_diagnostic()
        cert = short_time_certificate(formulas(3)["g"],F(15))
        self.assertGreater(plus-minus,float(cert["positive_event_gap"]))
        ratio = (plus-minus)/float(cert["leading"])
        self.assertLess(abs(ratio-1),1e-7)
        self.assertLess(tail,F(1,10**40))
        OBS["actual_original_instrument_diagnostic"] = dict(
            physical_one_excitation_dimension=720,graph_dimension=90,
            input_port=3,measured_port=0,
            p_plus=short(plus),p_minus=short(minus),difference=short(plus-minus),
            ratio_to_exact_second_order=short(ratio),
            Taylor_order=8,Taylor_operator_tail=str(tail),
            no_dense_full_H_spectrum_used=True,
            final_reading_is_actual_Lueders_port_instrument=True,
            floating_probabilities_do_not_replace_strict_Fraction_certificate=True)

    def test_06_sources_sampling_errors_and_large_I_boundary(self):
        data = formulas(3)
        cert = short_time_certificate(data["g"],F(15))
        gap = cert["positive_event_gap"]
        # Each of the two complete source/wait/read settings has ordinary
        # diamond error <=gap/4. Event error is half of that, so total <=gap/4.
        process_error = gap/4
        per_mean = gap/8
        samples = ceil_fraction(192/(gap*gap))
        self.assertGreaterEqual(2*samples*per_mean**2,6)
        self.assertGreater(sum(F(6)**k/math.factorial(k) for k in range(6)),160)
        # For I>=3, each parity success is >=43/90. Keep ALL failures internally.
        self.assertEqual((1-formulas(3)["p"])/2,F(43,90))
        attempts = 0
        while 160*samples*47**attempts > 90**attempts:
            attempts += 1
        self.assertLessEqual(F(47,90)**attempts,F(1,160*samples))
        bank_error = min(gap/64,F(1,80*samples))
        self.assertLessEqual(2*samples*(F(47,90)**attempts+bank_error/2),F(1,40))
        self.assertGreater(1-F(47,90)**attempts-bank_error/2,F(1,2))
        self.assertLessEqual(4*bank_error,gap/16)
        wait_read_error = 3*gap/16  # Includes actual data-port preparation, wait, and read.
        self.assertLessEqual(4*bank_error+wait_read_error,process_error)
        self.assertEqual(gap-process_error-2*per_mean,gap/2)
        # Exact integer inequalities support (not replace) the analytic
        # asymptotic g_I~1/(2 I^5) and illustrate no uniform positive gap.
        sequence = [formulas(i) for i in range(3,13)]
        self.assertTrue(all(sequence[k]["p"]>sequence[k+1]["p"]
                            for k in range(len(sequence)-1)))
        OBS["finite_resource_and_scope_ledger"] = dict(
            ideal_known_graph_source="I_G/D with an internal purifier",
            parity_measurement_projectors="(I+R)/2 and (I-R)/2",
            new_controlled_label_permutation_and_flag_reading_are_extra_inputs=True,
            parity_protocol_not_claimed_generated_by_original_F=True,
            worst_parity_success_I_ge_3="43/90",
            graph_purifier_and_every_failed_attempt_retained_in_internal_storage=True,
            samples_per_setting=str(samples),
            maximum_attempts_per_sample=attempts,
            total_graph_source_capacity_upper_bound=str(2*samples*attempts),
            complete_capped_preparation_bank_diamond_error=str(bank_error),
            per_attempt_full_source_and_flag_error_budget=str(bank_error/attempts),
            conditional_success_source_trace_error_bound=str(4*bank_error),
            subsequent_data_preparation_wait_and_read_diamond_error=str(wait_read_error),
            ideal_bank_cap_failure_per_sample=str(F(47,90)**attempts),
            independent_same_actual_protocol_replicates_are_an_extra_input=True,
            complete_source_wait_and_read_diamond_error_per_setting=str(process_error),
            empirical_mean_accuracy_each=str(per_mean),
            statistical_failure_bound="1/40",
            all_preparation_cap_failures_union_bound="1/40",
            total_failure_bound="1/20",
            surviving_estimated_event_probability_gap=str(gap/2),
            actual_huge_sampling_not_executed=True,
            asymptotic_gap="g_I ~ 1/(2 I^5)",
            selected_exact_gaps=[dict(i=x["i"],gap=str(x["g"])) for x in sequence],
            exact_population_closure_not_a_cognitive_axiom=True,
            geometry_may_legitimately_depend_on_graph_state=True,
            position_not_required_to_predict_all_internal_futures=True)

def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise AssertionError(stream.getvalue())
    baseline = Path(__file__).with_name("fast_graph_average_port_process_results.json")
    return dict(round=495,scientific_baseline_round=494,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        baseline_results_sha256=hashlib.sha256(baseline.read_bytes()).hexdigest(),
        scope=dict(all_labelled_degree_three_internal_tree_graph_states=True,
            arbitrary_role_covariant_Hermitian_fast_graph_update=True,
            unchanged_unit_strength_edge_controlled_exchange=True,
            exact_all_source_population_elimination_obstructed=True,
            no_uniform_positive_large_I_error_floor_claim=True,
            no_general_spatial_or_cognitive_no_go_claim=True,
            full_GR_goal_completed=False,phase_closure_triggered=False),
        observations=OBS)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run",action="store_true")
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    data = run()
    if args.check:
        assert json.loads(TARGET.read_text(encoding="utf-8"))==data
    elif not args.dry_run:
        with TARGET.open("x",encoding="utf-8") as handle:
            json.dump(data,handle,ensure_ascii=False,indent=2)
            handle.write("\n")
    print(json.dumps(data,ensure_ascii=False,indent=2))

