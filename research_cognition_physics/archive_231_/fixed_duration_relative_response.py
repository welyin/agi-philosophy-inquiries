"""Round 488: fixed-duration, signed relative distance response about a real baseline.

Baseline 487; actual 480 finite-pulse current state, unchanged source parameters,
original H and K, common clocks and independently replayed preparations.  No
inverse baseline unitary, free drift subtraction or position-group claim.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import io,json,math,platform,unittest
from pathlib import Path
import numpy as np
import rigid_leaf_reference_audit as frame
import symmetric_control_reachability_audit as words
import correlated_reference_local_access as prior
import singlet_supply_capacity_audit as fixed
import sequential_tree_distance_readout_audit as reader

TARGET=Path(__file__).with_name('fixed_duration_relative_response_results.json')
TAU=F(1,10**10)
BINV=22000
RADIUS=F(1,144*BINV)
COVER=RADIUS/(2*BINV)
OBS={}


def short(x):return float(f'{float(x):.12g}')


def rational_upper(x,denominator=10**20):
    return str(F((x.numerator*denominator+x.denominator-1)//x.denominator,denominator))


@lru_cache(None)
def certificate():
    h=frame.system()[1];cr,ci=prior.integer_columns();q=np.tile(words.qgraph(),(1,64)).astype(np.int64)
    exchange=np.argmax(words.s12(),axis=1);out=[[0]*3 for _ in range(3)];scale=1<<80;error=F(0);blocks=[]
    for excitation in range(7):
        idx=np.array([6*d+g for d in range(64) if d.bit_count()==excitation for g in range(6)])
        local={v:j for j,v in enumerate(idx)};swap=np.array([local[exchange[v]] for v in idx]);subh=h[np.ix_(idx,idx)]
        ur,ui,base,err,x=fixed.dyadic_unitary_certificate(subh,time=1,squarings=2,degree=14,bits=80)
        assert base==scale;assert all(isinstance(v,int) for v in ur.flat)
        error=max(error,err)
        def U(v):return ur@v[0]-ui@v[1],ur@v[1]+ui@v[0]
        def K(v):return v[0][swap],v[1][swap]
        def power(v,n):
            for _ in range(n):v=U(v)
            return v
        # Ignore the common -i phase of the old pi/2 exchange. Seven U factors.
        start=K(U((cr[idx],ci[idx])))
        value=power(start,6)
        for col,m in enumerate((1,2,3)):
            derivative=power(K(power(start,m+2)),4-m)
            derivative=(derivative[1],-derivative[0])
            for row in range(3):
                weights=q[row,idx,None].astype(object)
                numerator=np.sum(weights*(value[0]*derivative[0]+value[1]*derivative[1]))
                assert isinstance(numerator,int);out[row][col]+=2*numerator
        blocks.append(dict(excitation=excitation,dimension=len(idx),unitary_error=str(err)))
    jac=[[F(z,384*scale**14) for z in row] for row in out]
    inverse=prior.inverse3(jac);inverse_norm=max(sum(abs(z) for z in row) for row in inverse)
    word_error=(1+error)**7-1
    integer_entry_error=2*word_error*(2+word_error)
    pulse_entry_error=144*TAU
    beta=inverse_norm*3*(integer_entry_error+pulse_entry_error)
    return dict(jacobian=jac,inverse=inverse,inverse_norm=inverse_norm,
        determinant=prior.determinant(jac),unitary_error=error,word_error=word_error,
        integer_entry_error=integer_entry_error,pulse_entry_error=pulse_entry_error,beta=beta,blocks=blocks)


@lru_cache(None)
def numerics():
    h=frame.system()[1].astype(float);k=words.s12();q=np.tile(words.qgraph(),(1,64));u=words.unitary(h/5)
    c=prior.numerical_control()[0]@prior.pure_columns()/math.sqrt(384)
    return h,k,q,u,c


def output(theta):
    h,k,q,u,c=numerics();v=u@c
    for angle in theta:v=u@words.unitary(float(TAU)*h+float(angle)*k)@v
    values=q@np.sum(abs(v)**2,axis=1)
    return values,v


def finite_jacobian():
    h,k,q,u,c=numerics();e,v=np.linalg.eigh(h);ep=np.exp(-1j*float(TAU)*e)
    pulse=(v*ep)@v.conj().T
    # Exact Frechet integral at angle zero, evaluated stably by sinc.
    a=float(TAU)*(e[:,None]-e[None,:])/2
    kernel=np.exp(-1j*float(TAU)*(e[:,None]+e[None,:])/2)*np.sinc(a/np.pi)
    derivative=v@(-1j*kernel*(v.conj().T@k@v))@v.conj().T
    out=u@c
    for _ in range(3):out=u@pulse@out
    columns=[]
    for chosen in range(3):
        d=u@c
        for j in range(3):d=u@(derivative if j==chosen else pulse)@d
        columns.append(2*q@np.sum((out.conj()*d).real,axis=1))
    return np.column_stack(columns)


class Audit(unittest.TestCase):
    def test_01_same_current_state_and_fixed_clock_contract(self):
        h=frame.system()[1].astype(np.int64);k=words.s12().astype(np.int64);c=prior.pure_columns()
        self.assertEqual(int(np.sum(abs(c)**2)),384)
        self.assertEqual(max(sum(abs(int(z)) for z in row) for row in h),9)
        self.assertTrue(np.array_equal(k@k,np.eye(384)))
        self.assertLessEqual(np.max(abs(words.qgraph())),1)
        total=F(4,5)+3*TAU
        self.assertEqual(9*total,F(36,5)+27*TAU)
        OBS['fixed_process_contract']=dict(
            current_state='480 actual finite-pulse center, same fixed correlated source',
            new_free_waits=['1/5']*4,new_control_window_duration=str(TAU),
            duration_same_for_all_theta=str(total),baseline_continuation_u=0,
            baseline_continuation_cost=str(9*total),
            old_preparation_cost='9(3/5+tau)+pi/2 plus correlated source, clock, records',
            signed_additional_control_area='sum_j abs(theta_j)',
            baseline_and_perturbed_processes_use_same_source_and_clock_contract=True,
            baseline_is_a_real_independently_replayed_experiment=True,
            no_inverse_baseline_unitary_implemented=True,
            three_source_parameters_not_varied=True,
            original_H_remains_on_inside_all_four_old_and_new_pulses=True)

    def test_02_integer_full_rank_certificate_and_finite_pulses(self):
        cert=certificate();det=cert['determinant'];inv=cert['inverse_norm']
        self.assertLess(F(-1,10**8),det);self.assertLess(det,F(-9,10**9))
        self.assertLess(inv,11000)
        self.assertLess(cert['integer_entry_error'],F(1,10**13))
        self.assertLess(cert['beta'],F(1,2000))
        self.assertLess(inv/(1-cert['beta']),BINV)
        numerical=finite_jacobian();integer=np.array([[float(z) for z in row] for row in cert['jacobian']])
        self.assertLess(np.max(abs(numerical-integer)),float(cert['integer_entry_error']+cert['pulse_entry_error']))
        OBS['exact_response_certificate']=dict(
            rational_jacobian=[[str(z) for z in row] for row in cert['jacobian']],
            rational_determinant=str(cert['determinant']),
            rational_inverse_inf_norm=str(cert['inverse_norm']),
            rational_jacobian_numeric_display=[[short(z) for z in row] for row in cert['jacobian']],
            rational_determinant_interval=['-1/100000000','-9/1000000000'],
            rational_inverse_inf_norm_bound=11000,
            integer_unitary_error_rational_upper=rational_upper(cert['unitary_error']),
            seven_factor_word_error_rational_upper=rational_upper(cert['word_error']),
            integer_Jacobian_entry_error_rational_upper=rational_upper(cert['integer_entry_error']),
            finite_pulse_Jacobian_entry_error=str(cert['pulse_entry_error']),
            finite_pulse_and_integer_Banach_perturbation_beta_upper='1/2000',
            actual_finite_pulse_inverse_inf_norm_bound=BINV,
            blocks=cert['blocks'],
            all_certificate_multiplications_use_Python_integers=True,
            finite_Jacobian_crosscheck=[[short(z) for z in row] for row in numerical],
            finite_Jacobian_singular_values=[short(z) for z in np.linalg.svd(numerical,compute_uv=False)])

    def test_03_exact_neighborhood_and_vanishing_perturbation_cost(self):
        hessian=4;jac_lipschitz=3*3*hessian
        self.assertEqual(BINV*jac_lipschitz*RADIUS,F(1,4))
        self.assertEqual(BINV*COVER,F(1,2)*RADIUS)
        self.assertEqual(BINV*COVER+RADIUS/4,3*RADIUS/4)
        preparation_radius=F(1,12*BINV)
        self.assertEqual(6*BINV*preparation_radius,F(1,2))
        self.assertEqual((2*BINV)*36*(RADIUS/2),F(1,4))
        OBS['relative_endpoint_neighborhood']=dict(
            scalar_mixed_second_derivative_bound=hessian,Jacobian_inf_Lipschitz_constant=jac_lipschitz,
            parameter_cube_radius=str(RADIUS),relative_Q_cube_radius=str(COVER),
            contraction_constant_at_most='1/4',fixed_point_map_image_radius=str(3*RADIUS/4),
            actual_finite_pulse_Q_minus_actual_zero_control_baseline_covers_cube=True,
            parameter_inverse_Lipschitz_bound=f'4*{BINV}/3',
            incremental_control_area_bound=f'{4*BINV} * norm(relative_Q,infinity)',
            incremental_control_area_tends_to_zero=True,
            absolute_duration_and_maintenance_cost_do_not_tend_to_zero=True,
            relative_center_is_baseline_endpoint_not_original_current_Q=True,
            nearby_current_state_trace_norm_radius=str(preparation_radius),
            nearby_state_inverse_bound=2*BINV,nearby_state_parameter_radius=str(RADIUS/2),
            nearby_state_relative_Q_radius=str(COVER/4),
            nearby_state_targets_are_relative_to_its_own_baseline=True,
            nearby_state_leaf_protection_requires_separate_old_selection_contract=True,
            different_actual_states_can_require_different_calibrated_angle_solutions=True,
            no_physical_displacement_group_or_frame_change_law_inferred=True)

    def test_04_signed_actual_processes_and_common_leaf_reference(self):
        baseline,v0=output([0,0,0]);rho0=v0@v0.conj().T
        jac=finite_jacobian();inverse=np.linalg.inv(jac)
        step=float(COVER)/4;records=[];leaf_error=0.;finite_residual=0.
        for axis in range(3):
            for sign in (-1,1):
                target=np.zeros(3);target[axis]=sign*step
                theta=inverse@target;values,v=output(theta);difference=values-baseline
                residual=float(np.max(abs(difference-target)));finite_residual=max(finite_residual,residual)
                self.assertLess(residual,2e-14)
                self.assertLess(np.max(abs(theta)),float(RADIUS))
                # Original 480 source branches have 4 and 24 Gaussian columns.
                for col in (slice(0,4),slice(4,28)):
                    w=v[:,col];rho=2*w@w.conj().T;graph=rho.reshape(64,6,64,6).trace(axis1=0,axis2=2)
                    for ia,a in enumerate((0,3,4,5)):
                        for b in (0,3,4,5)[ia+1:]:leaf_error=max(leaf_error,abs(frame.distance(a,b)@np.diag(graph).real-8/3))
                records.append(dict(axis=axis,sign=sign,theta=[short(x) for x in theta],
                    relative_Q=[short(x) for x in difference],incremental_area=short(sum(abs(theta)))))
        self.assertLess(leaf_error,2e-12)
        # Unknown input/reference is preserved by the full pre-readout unitary.
        h,k,q,u,c=numerics();theta=np.array([.001,-.002,.003]);word=u.copy()
        for angle in theta:word=u@words.unitary(float(TAU)*h+angle*k)@word
        self.assertLess(np.linalg.norm(word.conj().T@word-np.eye(384)),2e-11)
        OBS['actual_signed_processes']=dict(baseline_Q=[short(x) for x in baseline],
            target_magnitude=str(COVER/4),records=records,
            maximum_linear_prediction_residual=short(finite_residual),
            maximum_leaf_reference_error_in_each_old_source_branch=short(leaf_error),
            all_pre_readout_evolution_unitary_for_arbitrary_unknown_reference=True,
            no_claim_arbitrary_reference_conditioned_leaf_moments_equal_8_over_3=True,
            tiny_numeric_endpoint_signs_not_used_as_proof=True)

    def test_05_two_real_processes_and_actual_CP_distance_readout(self):
        probe=F(1,65536);settings=([0,0,0],[.01,-.01,.02]);records=[];largest=0.
        for theta in settings:
            true,v=output(theta);rho=v@v.conj().T;graph=rho.reshape(64,6,64,6).trace(axis1=0,axis2=2)
            means=[]
            for leaf in (0,3,4,5):
                minus,plus,kraus=reader.instrument((1,leaf),probe)
                means.append(np.trace(reader.apply_map(plus-minus,graph,g=6,r=1)).real/float(probe))
                complete=np.einsum('bdgi,bdgj->ij',kraus.conj(),kraus,optimize=True)
                self.assertLess(np.max(abs(complete-np.eye(6))),3e-12)
            measured=means[0]-np.array(means[1:]);largest=max(largest,float(np.max(abs(measured-true))))
            records.append(dict(theta=list(theta),true_Q=[short(z) for z in true],measured_Q=[short(z) for z in measured]))
        bound=(math.expm1(18*float(probe))-18*float(probe))/float(probe)
        self.assertLess(largest,2*bound)
        OBS['two_real_instruments']=dict(probe_wait=str(probe),settings=records,
            total_actual_CP_instruments=8,largest_absolute_Q_readout_bias=short(largest),
            readout_difference_compares_two_measured_outputs=True,
            baseline_drift_is_not_a_free_computed_subtraction=True,
            old_D_transferred_to_noninteracting_S_G_remains_active_all_correlations_preserved=True,
            final_probe_readout_may_disturb_G=True,
            diagnostic_perturbation_outside_certified_small_cube_and_not_used_for_neighborhood_proof=True)

    def test_06_complete_source_control_readout_and_statistical_costs(self):
        epsilon=COVER/8;a=epsilon/8;source=epsilon/8;probe=a/(16*81);gamma=a*probe/8
        x=18*probe;tail=x*x/(2*(1-x/3));bias=(tail+gamma)/probe
        self.assertLess(bias,a/2)
        # Eight means, error a/2 for statistics; log(1600)<8.
        self.assertGreater(sum(F(8)**n/math.factorial(n) for n in range(40)),1600)
        copies=math.ceil(64/(a*a*probe*probe))
        self.assertEqual(2*source+4*a,3*epsilon/4)
        ctrl=gamma/2;duration=gamma/36;self.assertEqual(ctrl+18*duration,gamma)
        OBS['complete_resource_certificate']=dict(
            desired_relative_Q_infinity_error=str(epsilon),per_edge_total_error=str(a),
            per_process_common_source_and_all_controls_error_budget=str(source),
            four_edge_choices_share_one_actual_channel_in_each_process=True,
            two_process_channels_have_separate_errors_no_cancellation_assumed=True,
            signed_theta_sensitivity_to_total_angle_error='2 sum(abs(delta_theta))',
            certified_probe_wait=str(probe),instrument_error_budget=str(gamma),
            instrument_bias_bound=str(bias),instrument_control_only_error=str(ctrl),
            total_instrument_control_duration=str(duration),copies_per_mean=str(copies),
            total_independent_known_history_replays=str(8*copies),failure_probability_at_most='1/100',
            complete_relative_Q_error_bound=str(3*epsilon/4),
            source_clock_pulse_and_reference_costs_are_additional_explicit_contracts=True,
            no_cloning_of_one_unknown_current_system=True,
            arbitrary_target_precision_requires_increasing_source_control_and_readout_resources=True,
            giant_sampling_budget_not_executed=True)


def run():
    OBS.clear();stream=io.StringIO();result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise RuntimeError(stream.getvalue())
    output=dict(round=488,scientific_baseline_round=487,reused_frozen_rounds=[480,485,475,472],
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,observations=OBS,
        scope=dict(fixed_duration_actual_relative_three_direction_readout_neighborhood=True,
            baseline_maintenance_is_separately_costed=True,
            incremental_control_cost_not_total_physical_cost=True,
            not_a_rebuttal_of_485_or_487_absolute_small_action_bounds=True,
            single_unknown_copy_not_replayed_or_cloned=True,
            no_position_group_or_cross_basepoint_identification=True,
            no_necessary_three_dimensional_physical_space_derived=True,
            full_GR_goal_completed=False,phase_closure_triggered=False))
    assert json.loads(json.dumps(output))==output
    return output


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--dry-run',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();result=run()
    if args.check:assert result==json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
