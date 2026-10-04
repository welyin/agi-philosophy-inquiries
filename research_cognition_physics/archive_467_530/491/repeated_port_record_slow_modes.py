"""Round 491: actual refreshed port records and the six-tree slow-mode quotient.

Baseline 490. Fresh independent graph sources, port instruments, route control,
internal storage and clocks are inputs. This is not closed-H dissipation,
a physical time rescaling, or a derivation of three-dimensional space.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import io,itertools,json,math,platform,unittest
from pathlib import Path
import numpy as np
import port_spectral_speed_ruler as prior

TARGET=Path(__file__).with_name('repeated_port_record_slow_modes_results.json')
OBS={}
LEAVES=(0,3,4,5)
# Core H/projector certificates: at most four factors, dimension <=384, entries <=64.
assert 384**4*64**4 < 2**63
# General I=2,3 tree check: dimension <=8, |I*lap|<=810, |v|<=5.
assert 8*810*5 < 2**63


def short(x):return float(f'{float(x):.12g}')


@lru_cache(None)
def model():
    trees,_,_,adj,_,hj,hf,_=prior.model()
    h=(hj+hf+5*np.eye(36,dtype=np.int64)).astype(np.int64)
    leaf=np.array([1,0,0,1,1,1],dtype=np.int64)
    c=np.array([0,1,-1,0,0,0],dtype=np.int64)
    w=np.array([-1,2,2,-1,-1,-1],dtype=np.int64)
    pl4=4*np.diag(leaf)-np.outer(leaf,leaf)
    pc2=np.outer(c,c);pm12=np.outer(w,w)
    lap6=sum(np.diag(a.sum(axis=1))-a for a in adj)
    return trees,h,adj,lap6,pl4,pc2,pm12


def transition(delta):
    u=prior.unitary(model()[1],delta).reshape(6,6,6,6)
    conditional=np.sum(np.abs(u)**2,axis=1) # b,a,initial graph g
    return conditional.mean(axis=2),conditional,u


def certificates():
    delta=F(1,1024);s=delta**2;n=8*1024**2;gain=3**9
    remainder=F(6*18**4,math.factorial(4))*delta**4
    assert remainder<=s/32
    slow=(1-33*s/32,1-31*s/32)
    middle=(1-97*s/32,1-95*s/32)
    fast=(1-129*s/32,1-127*s/32)
    assert 0<fast[0]<=fast[1]<middle[0]<=middle[1]<slow[0]<slow[1]<1
    assert middle[1]/slow[0]<=1-s
    assert F(33,32)/(1-33*s/32)<=F(9,8)
    assert n*s==8
    exp8_lower=sum(F(8**k,math.factorial(k)) for k in range(21))
    assert exp8_lower>2000
    return delta,s,n,gain,remainder,slow,middle,fast,exp8_lower


def labeled_trees(internals):
    # Each degree-three internal label occurs exactly twice in a Prufer word.
    i=internals;n=2*i+2
    words=sorted(set(itertools.permutations(tuple(a for a in range(i) for _ in range(2)))))
    result=[]
    for word in words:
        degree=[1]*n
        for a in word:degree[a]+=1
        edges=[]
        for a in word:
            b=next(k for k in range(n) if degree[k]==1)
            edges.append((a,b));degree[a]-=1;degree[b]-=1
        a,b=(k for k in range(n) if degree[k]==1);edges.append((a,b))
        matrix=np.zeros((n,n),dtype=np.int64)
        for a,b in edges:matrix[a,b]=matrix[b,a]=1
        result.append(matrix)
    return result


class Audit(unittest.TestCase):
    def test_01_actual_instrument_chain_and_independent_source_requirement(self):
        trees,h,adj,*_=model();m,conditional,u=transition(.1)
        flat=u.reshape(36,36)
        projectors=[np.diag(np.repeat(np.eye(6,dtype=np.int64)[b],6)) for b in range(6)]
        self.assertTrue(np.array_equal(sum(projectors),np.eye(36,dtype=np.int64)))
        kraus=[p@flat for p in projectors]
        self.assertLess(np.linalg.norm(sum(k.conj().T@k for k in kraus)-np.eye(36)),1e-12)
        for a in range(6):
            history=np.array([[m[c,b]*m[b,a] for b in range(6)] for c in range(6)])
            self.assertAlmostEqual(float(history.sum()),1,places=12)
            self.assertTrue(np.allclose(history.sum(axis=1),(m@m)[:,a],atol=1e-13))
        # Fresh graphs with identical marginal I/6 but a shared classical g.
        correlated=float(np.mean(conditional[1,0]*conditional[3,1]))
        independent=float(m[1,0]*m[3,1])
        edges01=np.array([a[0,1] for a in adj]);edges13=np.array([a[1,3] for a in adj])
        joint=F(int(edges01@edges13),6)
        product=F(int(edges01.sum()),6)*F(int(edges13.sum()),6)
        self.assertEqual(joint,F(1,6));self.assertEqual(product,F(1,4))
        self.assertEqual(joint-product,-F(1,12));self.assertLess(correlated-independent,0)
        OBS['actual_record_chain']=dict(
            transition='M_delta(b,a)=sum_hg abs(U_(b,h),(a,g))^2/6',
            ideal_history_probability='product of M_delta(a_(k+1),a_k)',
            record_instrument='K_b=(|b><b| tensor I_G) U_delta; sum K_b^dag K_b=I',
            complete_history_isometry_retains_old_networks_purifiers_and_records=True,
            Markov_claim_only_for_classical_record_marginal=True,
            fresh_graph_must_be_independent_of_all_earlier_history=True,
            fresh_graph_marginal_I_over_6_alone_is_insufficient=True,
            correlated_two_step_path_delta4_coefficient=str(joint),
            independent_two_step_path_delta4_coefficient=str(product),
            correlated_minus_independent_delta4=str(joint-product),
            correlated_minus_independent_at_delta_point1=short(correlated-independent),
            old_networks_may_continue_own_H_while_isolated=True,
            no_closed_H_dissipation_or_free_graph_reset_claim=True)

    def test_02_integer_generator_symmetry_and_even_time(self):
        trees,h,adj,lap6,pl4,pc2,pm12=model()
        generator6=np.zeros((6,6),dtype=np.int64)
        for a in range(6):
            for b in range(6):
                if a!=b:generator6[b,a]=int(np.sum(h[6*b:6*b+6,6*a:6*a+6]**2))
            generator6[a,a]=-generator6[:,a].sum()
        self.assertTrue(np.array_equal(generator6,-lap6))
        self.assertTrue(np.array_equal(2*lap6,3*pl4+3*pm12+24*pc2))
        lookup={frozenset(tuple(sorted(e)) for e in t):i for i,t in enumerate(trees)}
        count=0
        for image in itertools.permutations(LEAVES):
            for swap in (False,True):
                p=list(range(6))
                for a,b in zip(LEAVES,image):p[a]=b
                if swap:p[1],p[2]=p[2],p[1]
                gp=[lookup[frozenset(tuple(sorted((p[a],p[b]))) for a,b in t)] for t in trees]
                idx=[6*p[a]+gp[g] for a in range(6) for g in range(6)]
                self.assertTrue(np.array_equal(h[np.ix_(idx,idx)],h));count+=1
        self.assertEqual(count,48)
        m,_,_=transition(.1);minus,_,_=transition(-.1)
        self.assertTrue(np.allclose(m,minus,atol=1e-14))
        self.assertTrue(np.allclose(m,m.T,atol=1e-14))
        self.assertTrue(np.allclose(m.sum(axis=0),1,atol=1e-14))
        OBS['short_step_origin']=dict(
            exact_generator_times_six=generator6.tolist(),
            expansion='M_delta=I-J^2 delta^2 L_average+O(delta^4)',
            odd_orders_zero_from_real_symmetric_H=True,
            exact_S4_times_S2_covariance_order=count,
            kappa_F_on_during_every_quantum_interval=True,
            graph_F_does_not_enter_second_order_port_population_generator=True,
            no_classical_Laplacian_added_to_H=True,
            remainder_bound_at_unit_couplings='norm_2(R)<=26244 delta^4')

    def test_03_fixed_step_spectral_order_and_finite_iteration_certificate(self):
        delta,s,n,gain,remainder,slow,middle,fast,exp8_lower=certificates()
        _,_,_,_,pl4,pc2,pm12=model();pl=pl4/4;pc=pc2/2;pm=pm12/12;p0=np.ones((6,6))/6
        m,_,_=transition(float(delta))
        lam_l=float(np.trace(pl@m)/3);lam_m=float(np.trace(pm@m));lam_c=float(np.trace(pc@m))
        self.assertLess(np.linalg.norm(m-(p0+lam_l*pl+lam_m*pm+lam_c*pc)),1e-12)
        for x,interval in zip((lam_l,lam_m,lam_c),(slow,middle,fast)):
            self.assertGreater(x,float(interval[0]));self.assertLess(x,float(interval[1]))
        self.assertLess(max((lam_m/lam_l)**n,(lam_c/lam_l)**n),1/2000)
        self.assertLess(lam_l**(-n),gain)
        OBS['finite_slow_mode_certificate']=dict(
            delta=str(delta),steps=n,quantum_waiting_time=str(n*delta),
            full_elapsed_time_also_includes_preparation_reading_routing_and_handoff=True,
            processing_parameter_n_delta_squared=str(n*s),
            operator_remainder_bound=str(remainder),
            operator_remainder_relative_bound='<=delta^2/32',
            slow_eigenvalue_interval=[str(x) for x in slow],
            middle_eigenvalue_interval=[str(x) for x in middle],
            fast_eigenvalue_interval=[str(x) for x in fast],
            exact_operator='lambda_L^(-n) (M_delta^n-P0)',
            distance_to_leaf_difference_projector_bound='less than 1/2000',
            normalization_gain_upper_bound=gain,
            diagnostic_eigenvalues=[short(lam_l),short(lam_m),short(lam_c)],
            diagnostic_gain=short(lam_l**(-n)),
            diagnostic_operator_error=short(max((lam_m/lam_l)**n,(lam_c/lam_l)**n)),
            exponential_lower_bound_via_finite_rational_series=str(exp8_lower),
            no_arbitrary_spectral_cutoff_or_changed_H=True,
            no_identification_of_n_delta_squared_with_physical_time=True)

    def test_04_endpoint_quotient_is_five_points_not_continuous_space(self):
        _,_,_,_,pl4,pc2,pm12=model()
        for numerator,denominator,rank in ((pl4,4,3),(pc2,2,1),(pm12,12,1)):
            self.assertTrue(np.array_equal(numerator@numerator,denominator*numerator))
            self.assertEqual(int(np.trace(numerator)),denominator*rank)
        self.assertFalse(np.any(pl4@pc2));self.assertFalse(np.any(pl4@pm12))
        self.assertFalse(np.any(pc2@pm12))
        self.assertTrue(np.array_equal(3*pl4+6*pc2+pm12,12*np.eye(6,dtype=np.int64)-2*np.ones((6,6),dtype=np.int64)))
        self.assertFalse(np.any(pl4[:,1]));self.assertFalse(np.any(pl4[:,2]))
        for a,b in itertools.combinations(LEAVES,2):
            v=pl4[:,a]-pl4[:,b];self.assertEqual(int(v@v),32)
        for a in LEAVES:self.assertEqual(int(pl4[:,a]@pl4[:,a]),12)
        OBS['endpoint_quotient']=dict(
            leaf_projector_times_four=pl4.tolist(),
            limiting_leaf_leaf_squared_distance='2',
            limiting_leaf_to_common_internal_center_squared_distance='3/4',
            both_internal_labels_have_same_limit=True,
            five_statistical_points_four_tetrahedron_vertices_and_center=True,
            affine_dimension_of_these_statistical_points=3,
            finite_n_centered_record_map_still_rank_five=True,
            no_continuous_position_neighborhood_translation_or_dimension_selection=True,
            normalization_is_data_processing_not_unit_success_physical_amplification=True,
            old_quantum_correlations_and_future_process_not_encoded_by_three_statistics=True)

    def test_05_general_symmetric_tree_source_does_not_select_three(self):
        checks=[]
        for i in (2,3):
            trees=labeled_trees(i);n=2*i+2;l=i+2;m=len(trees)
            total=sum(trees);lap=m*np.diag([3]*i+[1]*l)-total
            self.assertEqual(m,math.factorial(2*i)//(2**i))
            for a in range(n):
                for b in range(a+1,n):
                    expected=F(2,i) if b<i else F(1,i) if a<i else F(0)
                    self.assertEqual(F(int(total[a,b]),m),expected)
            for b in range(i+1,n):
                v=np.zeros(n,dtype=np.int64);v[i]=1;v[b]=-1
                self.assertTrue(np.array_equal(lap@v,m*v))
            for b in range(1,i):
                v=np.zeros(n,dtype=np.int64);v[0]=1;v[b]=-1
                self.assertTrue(np.array_equal(i*lap@v,m*(3*i+2)*v))
            v=np.array([l]*i+[-i]*l,dtype=np.int64)
            self.assertTrue(np.array_equal(i*lap@v,m*(2*i+2)*v))
            self.assertFalse(np.any(lap@np.ones(n,dtype=np.int64)))
            checks.append(dict(internal=i,leaves=l,source_trees=m,slow_multiplicity=i+1,
                nonzero_spectrum=['1 repeated '+str(i+1),str(F(2*i+2,i)),str(F(3*i+2,i))+' repeated '+str(i-1)]))
        OBS['general_tree_boundary']=dict(
            proven_family='I>=2 internal degree3, L=I+2 degree1 leaves, full maximally mixed source on all labelled trees of these degrees',
            leaf_internal_edge_probability='1/I',internal_internal_edge_probability='2/I',
            spectrum='0; 1 multiplicity I+1; 2+2/I multiplicity1; 3+2/I multiplicity I-1',
            exact_small_family_checks=checks,
            each_fixed_I_requires_own_sufficiently_small_delta=True,
            no_six_point_delta_claim_uniform_in_size=True,
            three_slow_components_come_from_I_equals_two_not_cognitive_principles=True,
            no_minimal_reconfigurable_tree_selection_axiom_added=True)

    def test_06_normalization_calibration_statistics_and_history_cost(self):
        delta,s,n,a,*_=certificates();eps=F(1,100)
        dlam=eps/(24*n*a);step=eps/(12*a*n);alpha=eps/(36*a)
        self.assertLessEqual(dlam,F(1,4*n))
        self.assertLessEqual(6*n*a*dlam,eps/4)
        self.assertLessEqual(3*a*n*step,eps/4)
        self.assertLessEqual(9*a*alpha,eps/4)
        self.assertLess(F(3,4)*eps+F(1,2000),eps)
        terminal_copies=math.ceil(8/alpha**2);calibration_copies=math.ceil(64/dlam**2)
        self.assertGreaterEqual(2*terminal_copies*alpha**2,16)
        self.assertGreaterEqual(2*calibration_copies*(dlam/4)**2,8)
        self.assertLess(F(72,2000**2)+F(4,2000),F(1,100))
        # Two active unit networks during handoff: norm <=18, diamond drift <=36*tau.
        control=step/2;handoff_window=step/72
        self.assertEqual(control+36*handoff_window,step)
        OBS['resource_contract']=dict(
            normalized_endpoint_vector_error_target=str(eps),
            normalized_truncation_error='less than 1/2000',
            calibrated_lambda_identity='lambda_L=M_delta(0,0)-M_delta(3,0)',
            calibrated_lambda_absolute_error=str(dlam),
            each_calibration_probability_source_plus_instrument_error=str(dlam/4),
            each_calibration_probability_statistical_error=str(dlam/4),
            perturbed_normalization_gain_at_most=3*a,
            gain_calibration_error_bound='6 n A delta_lambda <= epsilon/4',
            complete_each_step_channel_error=str(step),
            each_step_control_source_error_budget=str(control),
            two_active_network_handoff_window=str(handoff_window),
            all_history_source_error_bound='n epsilon_step; amplified <=epsilon/4',
            terminal_each_probability_statistical_error=str(alpha),
            six_inputs_times_six_terminal_effects=36,
            each_input_six_terminal_effects_share_one_actual_history_channel=True,
            final_classical_effect_choice_does_not_change_prior_source_or_process=True,
            terminal_histories_per_effect=terminal_copies,
            total_terminal_histories=36*terminal_copies,
            fresh_network_uses_for_terminal_histories=36*n*terminal_copies,
            two_calibration_effects=2,calibration_trials_per_effect=calibration_copies,
            combined_statistical_failure_below='1/100',
            no_experimental_trials_actually_run=True,
            every_network_graph_purifier_port_record_clock_route_and_isolation_counted=True,
            no_free_erasure_and_no_reuse_as_independent_source_of_same_unknown_network=True,
            old_hidden_records_never_used_as_extra_feedback=True)


def run():
    OBS.clear();stream=io.StringIO();result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise RuntimeError(stream.getvalue())
    answer=dict(round=491,scientific_baseline_round=490,reused_frozen_rounds=[265,267,280,283,442,443,465,472,475,490],
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,observations=OBS,
        scope=dict(actual_quantum_instruments_supply_this_specific_record_chain=True,
            diffusion_map_mathematics_reused_not_claimed_new=True,
            fresh_independent_network_sources_and_control_are_extra_inputs=True,
            no_closed_autonomous_H_dissipation_or_finite_resource_infinite_limit=True,
            no_physical_time_identification_of_delta_squared_steps=True,
            slow_quotient_does_not_derive_continuous_three_dimensional_space=True,
            full_GR_goal_completed=False,phase_closure_triggered=False))
    assert json.loads(json.dumps(answer))==answer
    return answer


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--dry-run',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();answer=run()
    if args.check:assert answer==json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as f:f.write(json.dumps(answer,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(answer,ensure_ascii=False,indent=2))
