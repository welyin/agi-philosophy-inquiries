"""Round 496: graph-blind limits at fixed waiting and port-read counts.

Scientific baseline 494, independent of 495. The full reader isometries inherit
465's weighted bound. All-scale results are analytic; finite matrices check
identities and examples. No large-I Hilbert space or long experiment is run.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib, io, json, math, platform, unittest
from pathlib import Path
import numpy as np
import fast_graph_average_port_process as prior
import branching_tree_distance_audit as tree_tools
import port_spectral_speed_ruler as spectral

ROOT=Path(__file__).resolve().parent
TARGET=ROOT/'graph_blind_scaling_limit_results.json'
OBS={}
# Exact int64 products below have dimension <=384, at most four factors,
# primitive absolute entries <=64. The huge finite witness uses Python ints.
assert 384**4*64**4 < 2**63


def short(x): return float(f'{float(x):.12g}')


@lru_cache(None)
def model():
    trees,full_h,f,adj,embedding,h,ports,n,degree,h0=prior.model()
    a=h-h0
    d=-np.kron(np.diag(degree),np.eye(6,dtype=np.int64))
    hf=np.kron(np.eye(6,dtype=np.int64),f)
    assert np.array_equal(h,5*np.eye(36,dtype=np.int64)+hf+d+a)
    return trees,h,ports,a,d,hf


def pinch(x):
    p=model()[2]
    return sum((z@x@z for z in p),np.zeros_like(x))


def comb(i):
    assert i>=2
    # Internals 1..i, source leaf 0, then leaves i+1..2i+1.
    edges={(j,j+1) for j in range(1,i)}|{(0,1)}
    nxt=i+1
    for j in range(1,i+1):
        count=1 if j==1 or 1<j<i else 2
        for _ in range(count):
            edges.add(tuple(sorted((j,nxt))));nxt+=1
    assert nxt==2*i+2
    return frozenset(edges)


def coefficients():
    delta=F(1,256)
    remainder=24+F(18,1-2*delta)
    lower=(1-48*delta)**2
    upper=(F(3,2)+48*delta)**2 # sqrt(2)<3/2
    return delta,remainder,lower,upper


def finite_witness():
    delta=F(1,256);n=16384;tau=n*delta**2;t=n*delta;b=F(41,4)
    exponent=2*b*t
    assert exponent.denominator==1
    r=56*3**int(exponent);i=56*r
    return dict(delta=delta,n=n,tau=tau,t=t,b=b,exponent=int(exponent),r=r,i=i)


class Audit(unittest.TestCase):
    def test_01_full_reader_weighted_intertwining_and_true_repeated_process(self):
        trees,h,ports,*_=model();u=spectral.unitary(h,1/256)
        identities=0
        for receiver in range(6):
            distances=np.array([len(tree_tools.path_between(g,6,a,receiver)) for a in range(6) for g in trees])
            integer_d=np.diag(distances)
            for p in ports:
                self.assertTrue(np.array_equal(integer_d@p,p@integer_d));identities+=1
            weights=1/distances;v=weights[:,None]*u/weights[None,:]
            gram=sum(((p@v).conj().T@(p@v) for p in ports),np.zeros((36,36),complex))
            self.assertLess(np.max(np.abs(gram-v.conj().T@v)),1e-13)
            self.assertLessEqual(np.linalg.norm(v,2),math.exp(41/(4*256))+1e-12)
        # Exact source label 0, graph correlated with a passive two-level R.
        psi=np.zeros((36,2),complex);psi[:6,:]=np.arange(1,13).reshape(6,2)+1j*np.arange(12).reshape(6,2)/5
        psi/=np.linalg.norm(psi);rho=np.outer(psi.ravel(),psi.ravel().conj())
        ur=np.kron(u,np.eye(2));pr=[np.kron(p,np.eye(2)) for p in ports]
        reference0=np.trace(rho.reshape(36,2,36,2),axis1=0,axis2=2)
        for _ in range(4):
            evolved=ur@rho@ur.conj().T
            rho=sum((p@evolved@p for p in pr),np.zeros_like(rho))
        reference1=np.trace(rho.reshape(36,2,36,2),axis1=0,axis2=2)
        self.assertLess(np.max(np.abs(reference1-reference0)),1e-13)
        populations=[np.trace(p@rho).real for p in pr]
        self.assertAlmostEqual(sum(populations),1,places=12)
        OBS['weighted_complete_reader_bridge']=dict(
            exact_commutations_checked=identities,
            isometry_identity='sum_c (P_c W U W^-1)^dagger(P_c W U W^-1)=(W U W^-1)^dagger(W U W^-1)',
            all_reader_histories_retained_in_proof=True,
            tail_uses_initial_fixed_graph_distance=True,
            repeat_bound='exp(B*n*delta), B=41/4',
            diagnostic_steps=4,diagnostic_delta='1/256',
            diagnostic_passive_reference_error=short(np.max(np.abs(reference1-reference0))),
            diagnostic_populations=[short(x) for x in populations],
            reduced_state_computation_is_not_physical_record_erasure=True)

    def test_02_uniform_flip_commutator_and_short_step_certificate(self):
        trees,h,ports,a,d,hf=model()
        comm=hf@a-a@hf
        self.assertLessEqual(int(np.abs(comm).sum(axis=1).max()),48)
        self.assertLessEqual(int(np.abs(comm).sum(axis=0).max()),48)
        sample=[]
        for ts,size in [(trees,6),(prior.larger_graph()[0],8)]:
            touched_max=row_max=0
            for g in ts:
                adj=tree_tools.adjacency(g,size)
                flips=tree_tools.tree_flips(g,size)
                for v in range(size):
                    touched=row=0
                    for gp,mult in flips.items():
                        ap=tree_tools.adjacency(gp,size)
                        change=len(adj[v]^ap[v])
                        self.assertIn(change,(0,2))
                        if change:touched+=mult;row+=change*mult
                    touched_max=max(touched_max,touched);row_max=max(row_max,row)
                    self.assertLessEqual(touched,24);self.assertLessEqual(row,48)
            sample.append(dict(vertices=size,graphs=len(ts),maximum_touched_flips=touched_max,maximum_commutator_row_sum=row_max))
        delta,rem,lo,up=coefficients()
        self.assertLess(rem,48);self.assertGreater(lo,F(1,2));self.assertLess(up,3)
        self.assertLessEqual(np.linalg.norm(d+a,2),6+1e-12)
        u=spectral.unitary(h,1/256)
        ui=spectral.unitary(hf,-1/256)@u*np.exp(5j/256)
        error=np.linalg.norm(ui-(np.eye(36)-1j*(d+a)/256),2)
        self.assertLess(error,float(rem*delta**2))
        OBS['uniform_step_certificate']=dict(
            general_active_NNI_count_per_vertex=24,changed_adjacency_entries_per_flip=2,
            general_commutator_operator_bound=48,small_exhaustive_checks=sample,
            interaction_generator_norm_at_most=6,
            maximal_delta=str(delta),second_order_remainder_coefficient=str(rem),
            adopted_remainder_coefficient=48,
            leaf_to_internal_probability_coefficient_at_least=str(lo),
            internal_to_leaf_probability_coefficient_at_most=str(up),
            internal_to_leaf_first_order_amplitude='sqrt(2), not sqrt(3)',
            diagnostic_interaction_picture_remainder=short(error),
            no_full_F_norm_used_in_uniform_short_step_bound=True)

    def test_03_conditional_graph_uniform_population_recursion(self):
        trees,h,ports,*_=model();delta=F(1,256);u=spectral.unitary(h,float(delta))
        pint=ports[1]+ports[2];pleaf=np.eye(36)-pint
        entries=[]
        for a in range(6):
            sl=slice(6*a,6*a+6)
            effect=u.conj().T@(pleaf if a in (1,2) else pint)@u
            ev=np.linalg.eigvalsh(effect[sl,sl])
            if a in (1,2):self.assertLessEqual(ev[-1],3*float(delta**2)+1e-13)
            else:self.assertGreaterEqual(ev[0],float(delta**2/2)-1e-13)
            entries.append(dict(port=a,minimum=short(ev[0]),maximum=short(ev[-1])))
        q=F(1,2)*delta**2;c=F(7,2)*delta**2
        self.assertGreater(1-c,0);self.assertEqual(q/c,F(1,7))
        low=F(0)
        for k in range(1,9):
            low=(1-c)*low+q
            self.assertEqual(low,(1-(1-c)**k)/7)
        exp_lower=1+F(7,8)+F(7,8)**2/2
        self.assertGreater(exp_lower,2)
        OBS['population_lower_bound']=dict(
            actual_graph_and_passive_reference_arbitrary=True,
            old_readers_may_be_conditionally_correlated_with_graph=True,
            bound_is_for_unconditional_population_after_summing_new_outcomes=True,
            lower_recursion='p_(k+1)>=(1-7*delta^2/2)*p_k+delta^2/2',
            lower_after_n='(1-exp(-7*n*delta^2/2))/7',
            fixed_tau='1/4',fixed_lower_bound='strictly greater than 1/14',
            exp_7_over_8_rational_lower=str(exp_lower),
            finite_operator_diagnostic=entries)

    def test_04_caterpillar_initial_labels_and_uniform_mean_target(self):
        examples=[]
        for i in (2,3,17):
            g=comb(i);adj=tree_tools.adjacency(g,2*i+2)
            self.assertEqual(len(g),2*i+1)
            self.assertEqual([len(adj[j]) for j in range(1,i+1)],[3]*i)
            self.assertTrue(all(len(adj[j])==1 for j in [0]+list(range(i+1,2*i+2))))
            self.assertEqual([len(tree_tools.path_between(g,2*i+2,0,j))-1 for j in range(1,i+1)],list(range(1,i+1)))
            r=max(1,i//2)
            self.assertLess(sum((F(1,(j+1)**2) for j in range(r+1,i+1)),F(0)),F(1,r+1))
            if i>=3:
                rate=F(2)+F(2,i)
                self.assertGreaterEqual(1/rate,F(3,8))
                self.assertGreaterEqual(rate,2)
            examples.append(dict(internal_nodes=i,total_vertices=2*i+2,internal_distances=list(range(1,i+1))))
        # Any target column has an R-element internal set with mass <=R/I.
        target=[F(1,20),F(1,10),F(1,5),F(3,20),F(2,5),F(1,10)]
        self.assertEqual(sum(target),1)
        internal=[1,2,3,4];chosen=sorted(internal,key=lambda j:target[j])[:2]
        self.assertLessEqual(sum(target[j] for j in chosen),F(2,len(internal)))
        self.assertEqual(F(3,8)*(1-F(2,3)),F(1,8))
        OBS['initial_comb_tail_and_target']=dict(
            examples=examples,
            initial_internal_vertex_j_distance_from_source_leaf='j',
            internal_outside_first_R_probability_at_most='exp(2*B*T)/(R+1)',
            no_exponential_sphere_count_used=True,
            mean_generator_internal_mass='(1-exp(-(2+2/I)*tau))/(2+2/I)',
            mean_target_lower_I_at_least_3='3*(1-exp(-2*tau))/8',
            mean_target_lower_at_tau_quarter='strictly greater than 1/8',
            arbitrary_graph_blind_target_need_not_be_covariant=True,
            target_minimum_mass_R_internal_labels_at_most='R/I',
            relabeling_source_and_detector_preserves_original_H_and_readout=True)

    def test_05_finite_all_size_separation_without_constructing_huge_graphs(self):
        p=finite_witness();r=p['r'];i=p['i']
        self.assertEqual(p['tau'],F(1,4));self.assertEqual(p['t'],64)
        self.assertEqual(p['exponent'],1312)
        tail=F(3**p['exponent'],r+1)
        self.assertLess(tail,F(1,56));self.assertEqual(F(r,i),F(1,56))
        lower=F(1,14)-F(1,56)-F(1,56)
        self.assertEqual(lower,F(1,28))
        mean_lower=F(1,8)*(1-F(r,i))-F(1,56)
        self.assertEqual(mean_lower,F(47,448))
        # The second actual graph relabels the first R internals to the last R.
        self.assertGreater(i-r+1,r)
        pair_lower=F(1,14)-2*F(1,56)
        self.assertEqual(pair_lower,F(1,28))
        OBS['finite_worst_case_witness']=dict(
            delta=str(p['delta']),port_instrument_calls=p['n'],coarse_parameter=str(p['tau']),
            original_H_waiting_time=str(p['t']),B=str(p['b']),twice_B_T=p['exponent'],
            proof_only_near_set_cardinality=str(r),internal_nodes=str(i),total_vertices=str(2*i+2),
            exact_definitions={'R':'56*3^1312','I':'56*R'},
            far_probability_strictly_less_than='1/56',
            any_graph_blind_target_worst_TV_strictly_greater_than=str(lower),
            uniform_mean_semigroup_TV_strictly_greater_than=str(mean_lower),
            two_actual_relabelled_comb_source_event_gap_strictly_greater_than=str(pair_lower),
            R_is_only_a_proof_event_not_a_physical_threshold=True,
            huge_graph_or_Hilbert_space_not_constructed=True,
            time_cost_bound='TV>=q0-exp(2BT)/(R+1)-R/I; choose R=floor(I*q0/4)',
            no_claim_bound_or_schedule_is_optimal=True)

    def test_06_actual_event_readout_error_and_complete_resource_scope(self):
        p=finite_witness();n=p['n'];i=p['i'];ports=2*i+2
        epsilon=F(1,224);source=epsilon/2;gamma=epsilon/(2*n)
        self.assertEqual(source+n*gamma,epsilon)
        alpha=F(1,224);samples=8*224**2
        self.assertEqual(2*samples*alpha**2,16)
        exp8_lower=sum(F(8**j,math.factorial(j)) for j in range(3))
        self.assertGreater(exp8_lower,20)
        self.assertLess(F(4,20**2),F(1,50))
        self.assertLess(F(4,exp8_lower**2),F(1,100))
        observed_gap=F(1,28)-epsilon-2*alpha
        self.assertEqual(observed_gap,F(5,224))
        bits=(ports-1).bit_length()
        self.assertLessEqual(ports,2**bits);self.assertGreater(ports,2**(bits-1))
        OBS['actual_readout_and_resources']=dict(
            event='reported endpoint label belongs to S',
            two_real_source_preparations='same leaf port; two role-relabelled comb graph basis states',
            each_trial_source_trace_error=str(source),each_complete_instrument_ordinary_diamond_error=str(gamma),
            total_ordinary_process_error_per_trial=str(epsilon),
            event_probability_error_per_trial=str(epsilon/2),
            independent_repetitions_per_source=samples,statistics_each_event_accuracy=str(alpha),
            joint_failure_probability_below='1/100',certified_empirical_event_gap_above=str(observed_gap),
            independent_full_trial_sources_are_explicit_inputs=True,
            port_reader_alphabet_size=str(ports),record_qubits_per_reader=bits,
            record_encoding_capacity_qubits_per_trial=n*bits,
            total_record_encoding_capacity_qubits_two_batches=2*samples*n*bits,
            encoding_capacity_is_not_full_detector_clock_control_or_network_cost=True,
            natural_wait_per_trial='64',reading_control_routing_and_preparation_times_are_additional=True,
            finite_reading_duration_is_included_in_complete_instrument_error=True,
            diagonal_occupation_reader_H_commutes_with_W_b_but_duration_adds_to_T=True,
            total_physical_time_and_all_hardware_cost_not_claimed_fixed=True,
            all_unknown_information_kept_in_the_reader_isometry=True,
            no_hidden_reader_feedback_or_postselection=True,
            no_actual_massive_network_preparation_or_statistical_run=True)


def run():
    OBS.clear();stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise AssertionError(stream.getvalue())
    answer=dict(round=496,scientific_baseline_round=494,reused_frozen_rounds=[442,465,490,491,492],
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        dependency_results_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
            ('fast_graph_average_port_process_results.json','relational_propagation_bound_audit_results.json','repeated_port_record_slow_modes_results.json')},
        observations=OBS,scope=dict(
            fixed_original_H_waiting_and_port_instrument_call_count=True,
            arbitrary_graph_blind_target_worst_error_nonvanishing_across_scale=True,
            only_one_excitation_sector_and_declared_port_Luders_interface=True,
            pure_graph_sources_already_suffice_for_lower_bound=True,
            large_size_and_fast_measurement_limits_not_interchanged=True,
            size_dependent_waiting_or_new_controls_not_excluded=True,
            graph_dependent_geometry_not_refuted=True,
            no_spatial_dimension_or_cognitive_principle_no_go=True,
            full_GR_goal_completed=False,phase_closure_triggered=False))
    assert json.loads(json.dumps(answer))==answer
    return answer


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();answer=run()
    if args.check:assert answer==json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as f:f.write(json.dumps(answer,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(answer,ensure_ascii=False,indent=2))
