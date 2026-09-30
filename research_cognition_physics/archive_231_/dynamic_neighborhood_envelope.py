"""Round 499: a fixed initial neighborhood envelope persists approximately.

No graph measurement or postselection is implemented. The projector in the
continuation lemma is only a comparison state for a uniform CPTP bound.
"""
import argparse
from fractions import Fraction as F
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import branching_tree_distance_audit as trees_model
import neighborhood_phase_boundary_obstruction as frozen497
import quantum_neighborhood_signal_audit as frozen443

TARGET=Path(__file__).with_name('dynamic_neighborhood_envelope_results.json')
OBS={}


def constants(c,r,j=1,kappa=1):
    assert c>=2 and r>=1
    s=c*sum((c-1)**z for z in range(r))
    lam=c*abs(j)+2*c*(c-1)**2*abs(kappa)
    a=s*c**(2*r)*(r+4)
    k=12*lam*(r+4)*c**6
    return s,lam,a,k


def root_ball(tree,n,b,r):
    return frozen443.ball(tree,b,r,n)


def path_count(tree,n,b,r,area):
    adj=trees_model.adjacency(tree,n)
    total=0
    for length in range(1,r+1):
        for end in itertools.permutations([v for v in range(n) if v!=b],length):
            if end[-1] in area:continue
            path=(b,)+end
            total+=int(all(y in adj[x] for x,y in zip(path[:-1],path[1:])))
    return total


def small_tree_families():
    six=trees_model.six_vertex_sector()[0]
    eight=sorted({trees_model.prufer_tree(8,w)
        for w in set(itertools.permutations((0,0,1,1,2,2)))},key=lambda g:sorted(g))
    return [(6,six),(8,eight)]


def rational_uniform_example():
    s,lam,a,k=constants(3,1)
    time=F(1,16*k);radius=31;m=(radius-1)//3+1
    # tanh(x)<=x and pi<22/7 give w<=11*K*t/7, with rational arithmetic.
    w=11*k*time/7
    q=2*a*w**m
    assert q<F(1,10**8)
    return dict(C=3,r=1,R=radius,zero_order=m,S=s,lam=lam,A=a,K=k,
                time=time,conformal_upper=w,q_upper=q)


def six_model():
    trees,full=frozen497.six_model()
    ids=np.array([(1<<(5-v))*6+g for v in range(6) for g in range(6)])
    h=full[np.ix_(ids,ids)]
    flip=trees_model.six_vertex_sector()[3].astype(np.int64)
    # Independent exact restriction of each data SWAP to one excitation.
    expected=np.kron(np.eye(6,dtype=np.int64),flip)
    for g,tree in enumerate(trees):
        adj=np.zeros((6,6),dtype=np.int64)
        for a,b in tree:adj[a,b]=adj[b,a]=1
        block=5*np.eye(6,dtype=np.int64)-np.diag(adj.sum(axis=1))+adj
        indices=np.arange(6)*6+g
        expected[np.ix_(indices,indices)]+=block
    assert np.array_equal(h,expected)
    return trees,h,full,flip


def diagnostics_data():
    trees,h,full,flip=six_model()
    g=frozen497.comb_pair(2)[1]
    other=frozenset({(0,1),(1,2),(1,4),(2,3),(2,5)})
    source_indices=[2*6+trees.index(g),3*6+trees.index(other)]
    bad=np.array([int(not root_ball(x,6,0,1).issubset({0,1})) for x in trees],dtype=np.int64)
    psi=np.zeros((36,2),complex)
    psi[source_indices[0],0]=1/math.sqrt(2)
    psi[source_indices[1],1]=1j/math.sqrt(2)
    return trees,h,full,flip,source_indices,bad,psi


def exact_effect_series(h,indices,effect,degree,time):
    # H has row/column absolute sums <=9. Each product is bounded BEFORE
    # int64 multiplication; a diagonal two-source expectation adds two terms.
    assert np.max(np.sum(abs(h),axis=0))<=9
    assert np.max(np.sum(abs(h),axis=1))<=9
    assert 4*18**degree<2**63
    a=np.diag(effect).astype(np.int64);coef=[]
    for n in range(degree+1):
        value=int(a[indices[0],indices[0]])+int(a[indices[1],indices[1]])
        assert n%2==0 or value==0
        coef.append(F(value*((-1)**(n//2)),2*math.factorial(n)) if n%2==0 else F(0))
        if n<degree:a=h@a-a@h
    polynomial=sum(c*time**n for n,c in enumerate(coef))
    tail=(18*time)**(degree+1)/(2*math.factorial(degree+1)*(1-18*time/(degree+2)))
    return coef,polynomial,tail


def evolution(h,time):
    values,vectors=np.linalg.eigh(h.astype(float))
    return (vectors*np.exp(-1j*time*values))@vectors.T


def receiver_reference(vector):
    x=vector.reshape(6,6,2)
    excited=np.einsum('gr,gs->rs',x[0],x[0].conj())
    ground=np.einsum('vgr,vgs->rs',x[1:],x[1:].conj())
    out=np.zeros((4,4),complex)
    out[:2,:2]=ground;out[2:,2:]=excited
    return out


def trace_distance(a,b):
    return float(np.linalg.svd(a-b,compute_uv=False).sum()/2)


def ceil_fraction(x):
    return -(-x.numerator//x.denominator)


class Audit(unittest.TestCase):
    def test_01_joint_event_is_rooted_path_count(self):
        cases=0;flips=0
        for n,family in small_tree_families():
            for tree in family:
                for radius in (1,2,3):
                    for area in ({0},{0,1,3}):
                        value=path_count(tree,n,0,radius,area)
                        outside=len(root_ball(tree,n,0,radius)-area)
                        self.assertEqual(value,outside)
                        self.assertLessEqual(int(outside>0),value)
                        self.assertLessEqual(value,constants(3,radius)[0])
                        cases+=1
                degrees=[len(x) for x in trees_model.adjacency(tree,n)]
                for neighbor,multiplicity in trees_model.tree_flips(tree,n).items():
                    self.assertEqual(multiplicity,1)
                    self.assertEqual(degrees,[len(x) for x in trees_model.adjacency(neighbor,n)])
                    flips+=1
        OBS['joint_path_event']=dict(exact_graph_region_radius_cases=cases,
            original_NNI_transitions=flips,
            M_equals_number_of_A_complement_vertices_in_current_ball=True,
            Q_bad_le_M_le_Sr_identity=True,
            no_sum_of_individual_N_dependent_probability_bounds=True)

    def test_02_initial_envelope_and_rooted_word_zero_region(self):
        rows=[]
        for outside_radius in (2,5,8):
            branch=outside_radius+2
            n,initial,_,_=frozen497.comb_pair(branch)
            area=root_ball(initial,n,0,outside_radius)
            max_flips=(outside_radius-1)//3
            frontier={initial};seen=set(frontier)
            for _ in range(max_flips):
                frontier={new for g in frontier for new in trees_model.tree_flips(g,n)}
                seen.update(frontier)
            self.assertLess(len(area),n)
            self.assertTrue(all(root_ball(g,n,0,1).issubset(area) for g in seen))
            rows.append(dict(R=outside_radius,r=1,N=n,checked_flip_depth=max_flips,
                             distinct_reached_graphs=len(seen),A_size=len(area)))
        trees,h,full,flip,indices,bad,psi=diagnostics_data()
        data=full-np.kron(np.eye(64,dtype=np.int64),flip)
        m=np.diag(np.tile(bad,64))
        self.assertTrue(np.array_equal(data@m,m@data))
        # Exact full H restriction, and nonzero graph leakage from a legal good state.
        self.assertTrue(all(bad[k%6]==0 for k in indices))
        amplitude=(h@psi)*np.tile(bad,6)[:,None]
        self.assertAlmostEqual(float(np.vdot(amplitude,amplitude).real),2,places=12)
        OBS['initial_support_and_dynamic_change']=dict(reachable_graph_crosschecks=rows,
            full_data_SWAP_commutes_with_graph_event=True,
            exact_one_excitation_restriction_matches_full_384_H=True,
            initial_good_subspace_not_claimed_invariant=True,
            leading_bad_probability_coefficient='2',
            general_zero_order='ad_H^n(M) sandwiched by initial Q_R is zero for r+3n<=R',
            finite_reachability_checks_do_not_replace_all_word_proof=True)

    def test_03_uniform_constants_and_all_time_certificate(self):
        checked=0
        for c,r in itertools.product((2,3,4),(1,2,3,5)):
            s,lam,a,k=constants(c,r)
            for n in range(1,21):
                size=r+3*n+1
                unoptimized=s*size*c**(2*(r+3*n))*2**n*(lam*size)**n
                self.assertLessEqual(unoptimized,a*k**n*math.factorial(n))
                checked+=1
        example=rational_uniform_example()
        self.assertEqual((example['S'],example['A'],example['K']),(3,135,1180980))
        self.assertEqual(example['zero_order'],11)
        n,initial,_,_=frozen497.comb_pair(33)
        area=root_ball(initial,n,0,31)
        self.assertEqual(n,68);self.assertLess(len(area),n)
        self.assertEqual(constants(2,5)[0],10)
        # C<=1 has no active NNI; kappa=0 leaves every graph diagonal observable fixed.
        self.assertEqual(trees_model.tree_flips(frozenset({(0,1)}),2),{})
        OBS['uniform_all_real_time_certificate']={key:str(value) for key,value in example.items()}
        OBS['uniform_all_real_time_certificate'].update(
            integer_constant_crosschecks=checked,
            nonempty_R31_example_N=68,nonempty_R31_example_A_size=len(area),
            theorem='q<=min(1,2*A_r*tanh(pi*K_r*abs(t)/2)^m), m=floor((R-r)/3)+1',
            all_real_times_come_from_unitary_strip_translation=True,
            C_at_most_one_or_zero_graph_coupling_gives_q_zero=True,
            initial_radius_must_satisfy_R_ge_r_ge_1=True,
            initial_label_envelope_not_derived=True)

    def test_04_finite_leakage_and_real_port_certificates(self):
        _,h,_,_,indices,bad,_=diagnostics_data()
        qc,qp,qe=exact_effect_series(h,indices,np.tile(bad,6),8,F(1,64))
        pc,pp,pe=exact_effect_series(h,indices,np.repeat([1,0,0,0,0,0],6),12,F(9,64))
        self.assertEqual(qc[:3],[F(0),F(0),F(2)])
        self.assertGreater(qp-qe,F(1,4096));self.assertLess(qp+qe,F(1,1000))
        self.assertGreater(pp-pe,F(1,10000))
        OBS['finite_actual_model_certificates']=dict(
            leakage_time='1/64',leakage_coefficients=[str(x) for x in qc],
            leakage_interval=[str(qp-qe),str(qp+qe)],
            actual_port_reading_time='9/64',port_coefficients=[str(x) for x in pc],
            port_probability_interval=[str(pp-pe),str(pp+pe)],
            port_probability_lower_bound='1/10000',
            largest_prior_int64_bound=4*18**12,
            graph_event_is_computed_only_not_physically_measured=True,
            physical_reading_is_two_outcome_projector_on_D0=True)

    def test_05_joint_reference_projection_comparison_only(self):
        _,h,_,_,_,bad,initial=diagnostics_data()
        current=evolution(h,1/64)@initial
        q=float(np.sum(abs(current*np.tile(bad,6)[:,None])**2))
        projected=current*np.tile(1-bad,6)[:,None]
        projected/=math.sqrt(1-q)
        overlap=abs(np.vdot(projected,current))**2
        self.assertAlmostEqual(overlap,1-q,places=12)
        current_density=np.outer(current.ravel(),current.ravel().conj())
        projected_density=np.outer(projected.ravel(),projected.ravel().conj())
        self.assertAlmostEqual(trace_distance(current_density,projected_density),math.sqrt(q),places=12)
        self.assertLess(np.linalg.norm(current.conj().T@current-initial.conj().T@initial),1e-12)
        later=evolution(h,1/8)@current
        later_projected=evolution(h,1/8)@projected
        actual=receiver_reference(later);comparison=receiver_reference(later_projected)
        joint_difference=trace_distance(actual,comparison)
        self.assertLessEqual(joint_difference,math.sqrt(q)+1e-12)
        actual_probability=float(actual[2:,2:].trace().real)
        comparison_probability=float(comparison[2:,2:].trace().real)
        self.assertGreater(actual_probability,F(1,10000))
        self.assertGreater(abs(actual_probability-comparison_probability),1e-5)
        self.assertLessEqual(abs(actual_probability-comparison_probability),joint_difference+1e-12)
        # Both physical port projectors together form a normalized CP instrument.
        self.assertAlmostEqual(float(np.trace(actual).real),1,places=12)
        OBS['reference_and_conditional_continuation']=dict(
            initial_reference_dimension=2,
            initial_source='(|2,G5,0_R>+i|3,G4,1_R>)/sqrt(2)',
            diagnostic_q=float(f'{q:.13g}'),
            full_DG_reference_distance=float(f'{math.sqrt(q):.13g}'),
            future_receiver_reference_distance=float(f'{joint_difference:.13g}'),
            actual_future_port_probability=float(f'{actual_probability:.13g}'),
            comparison_future_port_probability=float(f'{comparison_probability:.13g}'),
            projector_used_only_to_construct_mathematical_comparison=True,
            no_graph_measurement_no_postselection_no_source_refresh=True,
            conditional_predictor_lemma='good-domain D error eps implies actual error <=2*sqrt(q)+eps',
            unknown_reference_dimension_unrestricted_in_analytic_proof=True,
            no_specific_unfrozen_predictor_claimed=True)

    def test_06_noise_floor_and_finite_reading_resource_ledger(self):
        example=rational_uniform_example();q=example['q_upper']
        source_past_error=F(1,10**8)
        self.assertLess(4*(q+source_past_error),F(3,10000)**2)
        eps=F(1,1000);future_predictor_error=F(1,10000)
        self.assertEqual(F(3,10000)+eps+future_predictor_error,F(7,5000))
        g=F(1,10000);source=g/16;dynamic=g/16;clock=g/288;read=g/16
        self.assertEqual(source+dynamic+18*clock+read,g/4)
        source_gates=14
        self.assertEqual(g/32+source_gates*g/(32*source_gates),source)
        mean=g/8;samples=ceil_fraction(128/(g*g))
        self.assertGreaterEqual(2*samples*mean**2,4)
        self.assertGreater(sum(F(4)**j/math.factorial(j) for j in range(8)),40)
        self.assertEqual(g-(g/4)/2-mean,3*g/4)
        OBS['finite_error_and_resource_scope']=dict(
            conditional_uniform_example_past_trace_distance_budget=str(source_past_error),
            conditional_good_domain_predictor_error=str(eps),
            conditional_future_and_predictor_implementation_trace_distance_budget=str(future_predictor_error),
            conditional_actual_prediction_error_bound='7/5000',
            fixed_preparation_error_does_not_vanish_by_increasing_radius=True,
            noise_bad_event_includes_any_illegal_tree_sector_leakage=True,
            path_M_bound_applied_only_on_ideal_legal_tree_sector=True,
            noisy_comparison_requires_full_carrier_CPTP_implementations=True,
            demonstration_signal_lower_bound=str(g),
            demonstration_source_trace_error=str(source),
            demonstration_source_gates=source_gates,
            source_blank_trace_error=str(g/32),per_source_gate_diamond_error=str(g/(32*source_gates)),
            demonstration_dynamics_diamond_error=str(dynamic),
            demonstration_clock_tolerance=str(clock),
            clock_norm_9_used_after_ideal_source_replacement=True,
            demonstration_final_reader_diamond_error=str(read),
            complete_demonstration_process_trace_norm_error=str(g/4),
            samples=str(samples),empirical_mean_accuracy=str(mean),statistical_failure_bound='1/20',
            positive_empirical_port_frequency_lower_bound=str(3*g/4),
            listed_qubits_per_trial=23,total_listed_trial_qubits=str(23*samples),
            total_original_H_wait=str(samples*F(9,64)),
            phase_reference_control_isolation_clock_routing_and_blank_banks_are_extra_inputs=True,
            preparation_may_leave_tree_sector_before_legal_endpoint=True,
            retained_reference_and_all_measurement_records_not_erased=True,
            independent_same_actual_protocol_sampling_is_additional_contract=True,
            no_direct_sampling_of_Q_bad_and_no_gentle_projection_implementation=True,
            huge_sampling_and_control_hardware_not_executed=True)


def run():
    OBS.clear();stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise AssertionError(stream.getvalue())
    names=['neighborhood_phase_boundary_obstruction_results.json',
           'unknown_background_propagation_audit_results.json']
    return dict(round=499,scientific_baseline_round=497,reused_frozen_rounds=[442,443,466,497],
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        dependency_results_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in names},
        scope=dict(full_unknown_data_graph_coherence_and_passive_reference=True,
            fixed_initial_label_envelope_is_an_input=True,
            joint_escape_probability_bound_independent_of_N=True,
            approximate_envelope_not_exact_dynamic_support=True,
            projector_is_only_mathematical_comparison=True,
            specific_future_predictor_not_assumed_or_constructed=True,
            no_unfrozen_498_scientific_dependency=True,
            no_three_dimensional_geometry_or_autonomous_envelope_generation_claim=True,
            full_GR_goal_completed=False,phase_closure_triggered=False),observations=OBS)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--check',action='store_true');args=parser.parse_args();data=run()
    if args.check:assert json.loads(TARGET.read_text(encoding='utf-8'))==data
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8') as f:
            json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps(data,ensure_ascii=False,indent=2))
