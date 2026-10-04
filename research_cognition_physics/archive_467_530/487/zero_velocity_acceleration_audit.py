"""Round 487: zero initial distance velocity is not low-action reversibility.

Science uses frozen 474/479/480/485. A general two-source preparation family has
a control-independent initial acceleration. The all-control result uses
H+u(t) SWAP_12 and the declared budget 9*T+integral|u| only.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import correlated_reference_local_access as access
import positive_time_local_motion_cost as motion
import rigid_leaf_reference_audit as frame
import symmetric_control_reachability_audit as words
import sequential_tree_distance_readout_audit as reader

TARGET=Path(__file__).with_name('zero_velocity_acceleration_audit_results.json')
OBS={}
HNORM=9
BNORM=324
ALPHA=F(2,3)
CUTOFF=F(1,1944)
SCHEDULE=motion.SCHEDULE

def short(x):
    return float(f'{float(x):.13g}')

def gkron(a,b):
    ar,ai=a;br,bi=b
    return np.kron(ar,br)-np.kron(ai,bi),np.kron(ar,bi)+np.kron(ai,br)

def gadd(a,b):
    return a[0]+b[0],a[1]+b[1]

@lru_cache(None)
def data_coefficients():
    i=np.eye(2,dtype=np.int64);z=np.diag([1,-1]).astype(np.int64)
    x=np.array([[0,1],[1,0]],dtype=np.int64)
    y=np.array([[0,-1],[1,0]],dtype=np.int64)
    zero=np.zeros((2,2),dtype=np.int64)
    p=((i,zero),(x,zero),(zero,y),(z,zero))
    output=[]
    for source in p:
        factors=[p[0],source,p[0],gadd(p[0],p[1]),gadd(p[0],p[2]),gadd(p[0],p[3])]
        out=(np.ones((1,1),dtype=np.int64),np.zeros((1,1),dtype=np.int64))
        for factor in factors:out=gkron(out,factor)
        assert max(np.max(abs(v)) for v in out)<=2
        assert np.array_equal(out[0],out[0].T)
        assert np.array_equal(out[1],-out[1].T)
        output.append(out)
    return output

@lru_cache(None)
def operators():
    raw_h=frame.system()[1]; raw_k=words.s12(); raw_q=np.tile(words.qgraph(),(1,64))
    assert np.issubdtype(raw_h.dtype,np.integer)
    h=raw_h.astype(np.int64);k=raw_k.astype(np.int64);q=raw_q.astype(np.int64)
    assert np.array_equal(raw_k,k) and np.array_equal(raw_q,q)
    assert int(np.max(np.sum(abs(h),axis=1)))<=9 and int(np.max(abs(q)))<=1
    # c row norm <=18, second row norm <=324; all ensuing integer products are safe.
    assert 384**2*8*324<2**63
    first=[];second=[]
    for diagonal in q:
        c=h*diagonal[None,:]-diagonal[:,None]*h
        first.append(c)
        second.append(-(h@c-c@h))
    return h,k,q,first,second

def pair_density_coeff(source,uniform):
    graph=np.ones((6,6),dtype=np.int64) if uniform else np.eye(6,dtype=np.int64)
    dr,di=data_coefficients()[source]
    return np.kron(dr,graph),np.kron(di,graph)

def exact_moments(pair,denominator):
    rr,ri=pair;_,_,q,first,second=operators()
    # Hermitian Gaussian-integer densities (or affine coefficients).
    # i[H,Q] is imaginary, -[H,[H,Q]] real.
    assert rr.dtype==ri.dtype==np.dtype('int64')
    limit=384**2*int(max(np.max(abs(rr)),np.max(abs(ri)),1))*BNORM
    assert limit<2**63
    zeroth=[];velocity=[];acceleration=[]
    for diagonal,c,b in zip(q,first,second):
        self_imag=int(np.sum(np.diag(ri)*diagonal));assert self_imag==0
        zeroth.append(F(int(np.sum(np.diag(rr)*diagonal)),denominator))
        assert int(np.sum(rr*c.T))==0
        velocity.append(F(-int(np.sum(ri*c.T)),denominator))
        assert int(np.sum(ri*b.T))==0
        acceleration.append(F(int(np.sum(rr*b.T)),denominator))
    return zeroth,velocity,acceleration

@lru_cache(None)
def numerical_continuation():
    h,k,_,_,_=operators()
    initial=access.pure_columns()
    output=initial.copy()
    unitary=np.eye(384,dtype=complex)
    points=[motion.q_from_columns(output)]
    for dt,area in SCHEDULE:
        u=words.unitary(float(dt)*h+float(area)*k)
        output=u@output;unitary=u@unitary
        points.append(motion.q_from_columns(output))
    return initial,output,unitary,points

def ceil_fraction(x):
    return (x.numerator+x.denominator-1)//x.denominator

class Audit(unittest.TestCase):
    def test_01_exact_relative_degree_two_identities(self):
        h,k,q,c,b=operators()
        self.assertEqual(int(np.max(np.sum(abs(h),axis=1))),HNORM)
        self.assertTrue(np.array_equal(k@k,np.eye(384,dtype=np.int64)))
        graph_f=frame.system()[2]
        commutators=[]
        for j in range(3):
            self.assertLessEqual(int(np.max(abs(q[j]))),1)
            self.assertTrue(np.array_equal(k*q[j][None,:],q[j][:,None]*k))
            self.assertTrue(np.array_equal(k@c[j],c[j]@k))
            graph_q=words.qgraph()[j]
            expected=np.kron(np.eye(64,dtype=np.int64),graph_f*graph_q[None,:]-graph_q[:,None]*graph_f)
            self.assertTrue(np.array_equal(c[j],expected))
            self.assertTrue(np.array_equal(c[j],-c[j].T))
            self.assertTrue(np.array_equal(b[j].reshape(64,6,64,6).trace(axis1=1,axis2=3),np.zeros((64,64),dtype=np.int64)))
            self.assertTrue(np.array_equal(b[j],b[j].T))
            self.assertLessEqual(int(np.max(np.sum(abs(b[j]),axis=1))),BNORM)
            commutators.append(int(np.max(abs(k@b[j]-b[j]@k))))
        self.assertTrue(all(v>0 for v in commutators))
        OBS['exact_operator_identities']=dict(active_dimension=384,H_norm_bound=9,
            Q_order=['D13-D10','D14-D10','D15-D10'],
            K_commutes_with_Q_and_first_Heisenberg_derivative=True,
            second_derivative_operator='B_j=-[H,[H,Q_j]]',
            second_derivative_norm_bound=BNORM,
            K_second_derivative_commutator_max_entries=commutators,
            control_not_claimed_invisible_at_all_orders=True,
            Gaussian_integer_trace_safety_example_upper=384**2*8*BNORM,
            no_float_used_for_integer_certificates=True)

    def test_02_whole_two_source_family_and_independent_seed_certificate(self):
        records=[]
        for uniform in (True,False):
            branch=[]
            for source in range(4):
                q,v,b=exact_moments(pair_density_coeff(source,uniform),384)
                self.assertEqual(q,[F(0)]*3);self.assertEqual(v,[F(0)]*3)
                expected=[F(-4,3) if uniform and source==j+1 else F(0) for j in range(3)]
                self.assertEqual(b,expected)
                branch.append([str(x) for x in b])
            records.append(dict(graph='uniform_pure' if uniform else 'maximally_mixed',
                source_affine_basis=['I','X','Y','Z'],acceleration_coefficients=branch))
        parts=[pair_density_coeff(0,True),pair_density_coeff(1,True),
               pair_density_coeff(0,False),pair_density_coeff(2,False)]
        pair=(sum(v[0] for v in parts),sum(v[1] for v in parts))
        cr,ci=access.integer_columns()
        cr=cr.astype(np.int64);ci=ci.astype(np.int64)
        self.assertTrue(np.array_equal(pair[0],2*(cr@cr.T+ci@ci.T)))
        self.assertTrue(np.array_equal(pair[1],2*(ci@cr.T-cr@ci.T)))
        q,v,b=exact_moments(pair,768)
        self.assertEqual(b,[F(-2,3),F(0),F(0)])
        OBS['preparation_family']=dict(
            state='w eta(r) tensor |s><s| + (1-w) eta(s) tensor I6/6',
            domain='0<=w<=1, |r|<=1, |s|<=1',
            initial_Q_and_velocity_zero_for_entire_family=True,
            acceleration='-(4*w/3)*r', exact_affine_branch_records=records,
            seed_weight='1/2',seed_r=[1,0,0],seed_s=[0,1,0],
            seed_acceleration=[str(x) for x in b],
            seed_is_frozen_480_original_preparation_not_its_post_pulse_center=True,
            density_identity_with_frozen_480_columns_exact=True,
            full_arbitrary_source_family_not_proved_by_sampling=True)

    def test_03_uniform_all_control_cone_and_noise_boundary(self):
        self.assertEqual(2*BNORM*CUTOFF,ALPHA/2)
        self.assertEqual(ALPHA/4,F(1,6))
        examples=[]
        for w,rj in ((F(1,2),F(1)),(F(3,4),F(-1,2)),(F(1,7),F(2,3))):
            alpha=4*w*abs(rj)/3
            cutoff=alpha/(4*BNORM)
            self.assertEqual(-alpha+2*BNORM*cutoff,-alpha/2)
            self.assertEqual(cutoff,w*abs(rj)/972)
            examples.append(dict(w=str(w),r_j=str(rj),alpha=str(alpha),cutoff=str(cutoff),
                                 negative_quadratic_coefficient=str(alpha/4)))
        t=sum((v for v,_ in SCHEDULE),F(0))
        action=motion.action_budget(SCHEDULE)
        self.assertLess(action,CUTOFF)
        eta=F(1,10**16)
        robust=18*eta*t+(-ALPHA+BNORM*eta+2*BNORM*action)*t*t/2
        self.assertLess(robust,-t*t/6)
        OBS['all_control_acceleration_cone']=dict(
            controls='positive time, real L1 u, Hamiltonian H+u(t)K; fixed initial family member',
            declared_action='A=9*T+integral |u|',
            prefix_trace_distance_upper='2*A',
            exact_acceleration='q_j_double_dot(t)=Tr rho(t) B_j almost everywhere',
            aligned_coordinate='sign(r_j)*Q_j',
            family_if_w_abs_rj_positive='A<=w*abs(r_j)/972 implies signed_Delta_Q_j<=-w*abs(r_j)*T^2/3',
            family_exact_rational_examples=examples,
            seed_cutoff=str(CUTOFF),seed_bound='Delta_Qx<=-T^2/6<0',
            allowed_control_strength_not_bounded_pointwise=True,
            finite_words_and_L1_controls_covered_analytically=True,
            state_error_eta_bound='signed_Delta_Q_j <= 18*eta*T+(-alpha+324*eta+648*A)*T^2/2',
            fixed_nonzero_eta_does_not_give_uniform_sign_for_arbitrarily_small_T=True,
            robust_finite_example_eta=str(eta),robust_finite_example_upper=str(robust),
            no_arbitrary_generator_noise_theorem=True,
            changing_r_between_runs_is_changing_the_initial_source_not_continuation=True)

    def test_04_finite_signed_control_and_preserved_reference(self):
        initial,output,u,points=numerical_continuation()
        t=sum((dt for dt,_ in SCHEDULE),F(0));action=motion.action_budget(SCHEDULE)
        gap=t*t/6
        delta=points[-1][0]-points[0][0]
        self.assertLess(delta,-float(gap)+3e-15)
        stronger=(ALPHA-2*BNORM*action)*t*t/2
        self.assertLess(delta,-float(stronger)+3e-15)
        errors=[]
        for columns in (initial,output):
            for subset in (slice(0,4),slice(4,28)):
                g=motion.graph_from_columns(columns[:,subset],192)
                errors.append(max(abs(np.diag(g).real@frame.distance(a,b)-8/3)
                                  for a,b in itertools.combinations(frame.LEAVES,2)))
        self.assertLess(max(errors),5e-12)
        rng=np.random.default_rng(487)
        unknown=rng.normal(size=(384,3))+1j*rng.normal(size=(384,3))
        unknown/=np.linalg.norm(unknown)
        recovery=np.linalg.norm(u.conj().T@(u@unknown)-unknown)
        self.assertLess(recovery,4e-12)
        OBS['finite_control_diagnostic']=dict(schedule=[dict(time=str(dt),area=str(a)) for dt,a in SCHEDULE],
            time=str(t),action=str(action),strict_negative_gap=str(gap),
            stronger_negative_gap=str(stronger),actual_delta_Qx=short(delta),
            Q_at_segment_endpoints=[[short(x) for x in point] for point in points],
            max_conditional_leaf_mean_error=short(max(errors)),
            reference_matching_protection_follows_475_477_symmetry_not_round486=True,
            unknown_R_recovery_residual=short(recovery),
            mathematical_inverse_not_claimed_as_allowed_negative_time=True,
            arbitrary_R_information_preservation_not_arbitrary_R_leaf_moment_protection=True,
            original_data_graph_correlations_retained=True)

    def test_05_degenerate_family_and_scope_boundary(self):
        # w=0, s=e_x has zero initial acceleration, but is not H-stationary.
        pair=gadd(pair_density_coeff(0,False),pair_density_coeff(1,False))
        h,_,_,_,_=operators()
        commr=h@pair[0]-pair[0]@h
        commi=h@pair[1]-pair[1]@h
        self.assertGreater(int(np.max(abs(commr)))+int(np.max(abs(commi))),0)
        q,v,b=exact_moments(pair,384)
        self.assertEqual(q,[F(0)]*3);self.assertEqual(v,[F(0)]*3);self.assertEqual(b,[F(0)]*3)
        frozen=json.loads(words.TARGET.read_text(encoding='utf-8'))
        self.assertEqual((frozen['round'],frozen['tests_run'],frozen['failures'],frozen['errors']),(479,6,0,0))
        self.assertTrue(frozen['scope']['exact_two_coefficient_covariant_response_classification'])
        skew=words.K
        self.assertTrue(np.array_equal(skew.T,-skew))
        self.assertTrue(np.array_equal(skew@np.ones(3,dtype=np.int64),np.zeros(3,dtype=np.int64)))
        # Independent original-H diagnostic of the degenerate mixed-graph branch.
        # Its nonstationarity does not guarantee 3D access; the old theorem fixes a plane.
        state=(pair[0]+1j*pair[1])/384
        unitary=np.eye(384,dtype=complex)
        for dt,area in ((.13,.27),(.21,-.19),(.08,.35)):
            unitary=words.unitary(dt*h+area*operators()[1])@unitary
        after=unitary@state@unitary.conj().T
        final_q=np.tile(words.qgraph(),(1,64))@np.diag(after).real
        self.assertLess(abs(final_q[1]+final_q[2]),4e-14)
        self.assertGreater(np.linalg.norm(final_q),1e-8)
        OBS['scope_boundary']=dict(
            family_degenerate_case='w=0 or r=0 has zero initial acceleration',
            zero_acceleration_example='w=0, s=e_x',
            exact_nonstationary_H_commutator_max=[int(np.max(abs(commr))),int(np.max(abs(commi)))],
            frozen_round479_two_coefficient_classification_reused=True,
            degenerate_family_all_control_outputs='(1-w)*(a*s+b*K_skew*s)',
            degenerate_family_fixed_plane='span{s,K_skew*s}, dimension at most 2',
            degenerate_s_parallel_ones_at_most_one_dimension=True,
            degenerate_s_zero_or_w_one_gives_zero_output=True,
            exact_full_plane_reachability_not_claimed=True,
            mixed_graph_finite_control_Q=[short(x) for x in final_q],
            entire_family_dichotomy='w*r nonzero: low-action halfspace; w*r=0: all-action fixed plane',
            all_fixed_two_source_family_members_fail_continuous_low_action_3D_neighborhood_under_declared_contract=True,
            no_claim_all_instantaneously_stationary_states_fail=True,
            no_claim_all_control_classes_fail=True,
            no_claim_no_future_information_without_position_summary=True,
            later_endpoint_access_from_seed_not_proved_by_481_current_chart_theorem=True,
            remaining_open_case='zero acceleration states outside this family or broader sourced relational controls')

    def test_06_true_CP_readout_and_complete_finite_budget(self):
        initial,output,_,_=numerical_continuation()
        probe=F(1,65536)
        true=[];measured=[];records=[]
        for columns in (initial,output):
            graph=motion.graph_from_columns(columns)
            true.append(float(words.qgraph()[0]@np.diag(graph).real))
            local=[];means=[]
            for leaf in (0,3):
                minus,plus,_=reader.instrument((1,leaf),probe)
                probabilities=[float(np.trace(np.einsum('ghij,ij->gh',op,graph)).real)
                               for op in (minus,plus)]
                self.assertLess(abs(sum(probabilities)-1),5e-12)
                self.assertGreaterEqual(min(probabilities),-5e-12)
                means.append((probabilities[1]-probabilities[0])/float(probe))
                local.append(dict(ordered_edge=[1,leaf],measured_port=1,
                                  probabilities=[short(x) for x in probabilities]))
            measured.append(means[0]-means[1]);records.append(local)
        x=18*probe;bias=x*x/(2*(1-x/3)*probe)
        self.assertLess(max(abs(a-b) for a,b in zip(true,measured)),float(2*bias)+3e-9)
        t=sum((v for v,_ in SCHEDULE),F(0))
        gap=t*t/6;source=gap/8;a=gap/16;h=a/(16*81);gamma=a*h/4
        x=18*h
        self.assertLessEqual((x*x/(2*(1-x/3))+gamma)/h,a/2)
        copies=ceil_fraction(80/(h*h*a*a))
        self.assertGreaterEqual(copies*h*h*a*a/8,10)
        self.assertGreater(sum(F(10**j,math.factorial(j)) for j in range(10)),800)
        remaining=gap-2*source-4*a
        self.assertEqual(remaining,gap/2)
        OBS['actual_readout_and_resources']=dict(
            actual_CP_instrument_round=472,diagnostic_probe_wait=str(probe),
            diagnostic_records=records,true_Qx_before_after=[short(x) for x in true],
            instrument_Qx_before_after=[short(x) for x in measured],
            diagnostic_each_Qx_bias_upper=str(2*bias),
            diagnostic_probe_not_claimed_to_certify_tiny_quadratic_gap=True,
            rigorous_negative_gap=str(gap),each_setting_complete_source_diamond_error=str(source),
            same_actual_prepared_state_for_both_edge_readouts_in_each_setting=True,
            edge_selected_after_common_source_preparation_and_without_old_record_feedback=True,
            differing_edge_specific_sources_need_separate_error_accounting=True,
            each_adjacency_mean_error=str(a),certified_probe_wait=str(h),
            complete_storage_handoff_and_probe_instrument_diamond_error=str(gamma),
            independent_means=4,copies_per_mean=str(copies),total_known_history_preparations=str(4*copies),
            simultaneous_failure_probability_less_than='1/100',
            negative_measured_difference_magnitude_lower=str(remaining),
            old_complete_data_moved_to_isolated_internal_storage_and_G_remains_active=True,
            old_storage_graph_record_reference_correlations_retained=True,
            final_CP_readout_can_disturb_G=True,
            no_clone_of_unknown_state_or_enormous_sampling_run=True,
            source_control_clock_storage_and_repetition_are_explicit_extra_inputs=True,
            selected_action_not_identified_with_unique_cognitive_cost_or_dissipation=True)

def run():
    OBS.clear();stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise RuntimeError(stream.getvalue())
    output=dict(round=487,baseline_round=485,scientific_baselines=[474,479,480,485],
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,observations=OBS,
        scope=dict(zero_initial_velocity_not_sufficient_for_low_action_two_sided_motion=True,
            explicit_two_source_family_only=True,
            no_complete_cognitive_or_spatial_no_go=True,
            independent_of_parallel_round486=True,
            full_GR_goal_completed=False,phase_closure_triggered=False))
    assert json.loads(json.dumps(output))==output
    return output

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();r=run()
    if args.check:assert r==json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(r,ensure_ascii=False,indent=2))
