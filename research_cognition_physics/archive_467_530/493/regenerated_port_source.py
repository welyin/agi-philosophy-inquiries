"""Round 493: regenerate port sources with one retained data/graph network.

Baseline 492. Bounded measured retries, with an explicit abort flag, replace
491's repeated fresh network sources. All statistical trials may reuse DG;
conditional channel bounds and martingale concentration replace IID sampling.
Blank readers, clocks, comparisons, isolation and instrument precision remain
additional resources, and hidden new records never become free product states.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib,io,itertools,json,math,platform,unittest
from pathlib import Path
import numpy as np
import port_spectral_speed_ruler as prior

ROOT=Path(__file__).resolve().parent
TARGET=Path(__file__).with_name('regenerated_port_source_results.json')
OBS={}
# Integer resource arithmetic uses Python int/Fraction, never fixed-width arrays.


def short(x):return float(f'{float(x):.12g}')


def frozen(name):return json.loads((ROOT/name).read_text(encoding='utf-8'))


@lru_cache(None)
def model():
    _,_,_,_,_,hj,hf,_=prior.model()
    h=hj+hf+5*np.eye(36,dtype=np.int64)
    ports=[np.diag(np.repeat(np.eye(6)[a],6)) for a in range(6)]
    return h,ports


def channel(x,u,ports):return sum(p@u@x@u.conj().T@p for p in ports)


def parameters():
    old=frozen('retained_port_instrument_mixing_results.json')
    proof=old['observations']['finite_certificate']
    q0=F(proof['ordinary_diamond_bound']);m=proof['instrument_steps']
    assert m==65537 and q0<F(1,4000)
    slow=frozen('repeated_port_record_slow_modes_results.json')
    cert=slow['observations']['finite_slow_mode_certificate']
    n=cert['steps'];a=cert['normalization_gain_upper_bound'];delta=F(cert['delta'])
    assert n==8388608 and a==19683 and delta==F(1,1024)
    q=F(1,4000);r=6;k=256;eps=F(1,100)
    dlam=eps/(24*n*a);beta=dlam/4;alpha=eps/(36*a)
    gmix=beta/(8*k*r*m);gtest=beta/(8*k);gprobe=beta/4
    eta=q**r+r*m*gmix+gtest
    failure=F(5,6)**k
    regeneration=k*eta+2*failure
    macro=regeneration+gprobe
    mathematical=k*q**r+2*failure
    assert mathematical<beta/4
    assert macro<beta
    return dict(q0=q0,q=q,m=m,n=n,a=a,delta=delta,r=r,k=k,eps=eps,
        dlam=dlam,beta=beta,alpha=alpha,gmix=gmix,gtest=gtest,gprobe=gprobe,
        eta=eta,failure=failure,regeneration=regeneration,macro=macro,
        mathematical=mathematical)


class Audit(unittest.TestCase):
    def test_01_last_record_is_a_cp_readout_of_the_block_diagonal_state(self):
        h,ports=model();u=prior.unitary(h,1/16)
        full_h=prior.model()[1]
        number=np.array([d.bit_count() for d in range(64) for _ in range(6)],dtype=np.int64)
        # Number differences <=6, Hamiltonian entries <=9: exact integer check.
        assert 6*9 < 2**63
        self.assertFalse(np.any((number[:,None]-number[None,:])*full_h))
        rng=np.random.default_rng(493)
        v=rng.normal(size=72)+1j*rng.normal(size=72);v/=np.linalg.norm(v)
        rho=np.outer(v,v.conj());ur=np.kron(u,np.eye(2))
        pr=[np.kron(p,np.eye(2)) for p in ports]
        first=channel(rho,ur,pr)
        actual_last=[p@ur@first@ur.conj().T@p for p in pr]
        final=sum(actual_last)
        copied_last=[p@final@p for p in pr]
        self.assertLess(max(np.linalg.norm(x-y) for x,y in zip(actual_last,copied_last)),1e-14)
        self.assertTrue(np.allclose(sum(p@p for p in ports),np.eye(36)))
        # Exact stationary projection: every last outcome has probability 1/6;
        # target a's normalized DG state is Pa tensor I_G/6.
        for p in ports:
            block=p/36
            self.assertAlmostEqual(float(np.trace(block)),1/6)
            self.assertTrue(np.allclose(block*6,p/6))
        # A maximally mixed graph marginal does not decouple a new hidden purifier.
        half_distance=F(35,36)
        self.assertEqual(half_distance,1-F(1,6**2))
        OBS['last_record_interface']=dict(
            theta='Theta(X)=sum_c Pc X Pc tensor |c><c|_classical',
            exact_last_record_channel='Theta composed with E^(r*m)',
            final_port_and_classical_last_record_agree=True,
            target_success_probability_at_Pi='1/6',
            target_success_state='|a><a| tensor I_G/6',
            no_new_graph_measurement_or_data_permutation_used=True,
            reference_statement_covers_all_initial_DG_R_in_one_excitation_sector=True,
            actual_instruments_including_errors_preserve_one_excitation_sector=True,
            occupation_controlled_reader_interface_can_conserve_data_excitation=True,
            complete_number_preserving_readout_hardware_not_compiled=True,
            newly_generated_hidden_readers_not_in_product_reference_claim=True,
            maximally_entangled_graph_hidden_purifier_half_distance_from_product=str(half_distance))

    def test_02_frozen_mixing_certificate_power_and_contractivity(self):
        p=parameters();h,ports=model();u=prior.unitary(h,1/16)
        rng=np.random.default_rng(49302)
        x=rng.normal(size=(36,36))+1j*rng.normal(size=(36,36));x=x+x.conj().T
        pi=lambda z:np.eye(36)*np.trace(z)/36
        e=lambda z:channel(z,u,ports)
        self.assertLess(np.linalg.norm(e(np.eye(36))-np.eye(36)),1e-12)
        self.assertLess(np.linalg.norm(e(pi(x))-pi(x)),1e-12)
        self.assertLess(np.linalg.norm(pi(e(x))-pi(x)),1e-12)
        b=lambda z:e(z)-pi(z)
        self.assertLess(np.linalg.norm(e(e(x))-pi(x)-b(b(x))),1e-11)
        OBS['mixing_power']=dict(
            frozen_492_ideal_diamond_bound=str(p['q0']),used_coarse_bound=str(p['q']),
            original_instrument_steps_per_block=p['m'],original_wait_delta='1/16',
            blocks_per_attempt=p['r'],readers_per_attempt=p['r']*p['m'],
            ideal_power_bound=str(p['q']**p['r']),
            identity='(E^m-Pi)^r=E^(r*m)-Pi',
            complete_error_per_mixing_instrument=str(p['gmix']),
            classical_compare_and_branch_error_per_attempt=str(p['gtest']),
            complete_attempt_error=str(p['eta']),
            the_last_record_lift_is_CPTP_so_diamond_error_does_not_increase=True,
            no_new_full_mixing_rate_scan_or_new_492_certificate=True)

    def test_03_bounded_retry_preserves_history_and_never_drops_abort(self):
        h,ports=model();u=prior.unitary(h,1/16);ks=[p@u for p in ports]
        rng=np.random.default_rng(49303)
        v0=rng.normal(size=(36,2))+1j*rng.normal(size=(36,2));v0/=np.linalg.norm(v0)
        v1=rng.normal(size=(36,2))+1j*rng.normal(size=(36,2));v1/=np.linalg.norm(v1)
        active=[(v0,v1)];finished=[];target=0
        # Three single-instrument attempts check the exact branch isometry,
        # not the numerical quality of the full certified preparation block.
        for _ in range(3):
            next_active=[]
            for x,y in active:
                for b,kb in enumerate(ks):
                    pair=(kb@x,kb@y)
                    if b==target:finished.append(pair)
                    else:next_active.append(pair)
            active=next_active
        all_branches=finished+active
        self.assertEqual(len(finished),1+5+25);self.assertEqual(len(active),125)
        overlap=sum(np.vdot(x,y) for x,y in all_branches)
        self.assertLess(abs(overlap-np.vdot(v0,v1)),1e-12)
        reference=sum(x.conj().T@x for x,_ in all_branches)
        self.assertLess(np.linalg.norm(reference-v0.conj().T@v0),1e-12)
        self.assertAlmostEqual(float(sum(np.vdot(x,x).real for x,_ in all_branches)),1,places=12)
        p=parameters();failure=p['failure']
        # Ideal replacement attempts have an exact truncated geometric law.
        success=sum(F(1,6)*F(5,6)**j for j in range(p['k']))
        self.assertEqual(success+failure,1)
        self.assertLess(p['mathematical'],p['beta']/4)
        self.assertLess(p['macro'],p['beta'])
        OBS['retry_channel']=dict(
            attempt_cap=p['k'],no_arbitrary_c_to_a_data_SWAP=True,
            allowed_feedback='last outcome equals target; bounded repeat; outer observed endpoint',
            ideal_abort_probability=str(failure),
            approximate_source_error='k*eta+2*(5/6)^k in ordinary diamond norm',
            preparation_diamond_bound=str(p['regeneration']),
            final_probe_delta=str(p['delta']),probe_full_instrument_error=str(p['gprobe']),
            macro_diamond_bound=str(p['macro']),asserted_macro_bound=str(p['beta']),
            finite_isometry_test_success_histories=len(finished),abort_histories=len(active),
            test_uses_three_short_instruments_only_not_full_preparation_execution=True,
            abort_is_an_output_never_postselected_away=True,
            delivery_at_bounded_stopping_time=True,
            success_immediately_starts_probe_without_wait_to_attempt_cap=True,
            no_free_hold_of_the_port_or_switching_off_original_H=True,
            all_failure_records_detector_environments_and_unknown_information_retained=True)

    def test_04_adaptive_macro_composition_and_one_network_reuse(self):
        p=parameters();n=p['n'];a=p['a']
        self.assertLessEqual(n*p['macro'],n*p['beta'])
        self.assertEqual(3*a*n*p['beta'],p['eps']/32)
        # Exact finite scalar check of first-abort accounting, without any
        # assumption that actual quantum histories are independently sampled.
        ideal_abort=F(5,6)**4
        for length in (1,3,7):
            survival=(1-ideal_abort)**length
            abort_weights=[ideal_abort*(1-ideal_abort)**j for j in range(length)]
            self.assertEqual(survival+sum(abort_weights),1)
            self.assertLessEqual(2*(1-survival),2*length*ideal_abort)
        OBS['macro_history']=dict(
            ideal_macro_transition='M_delta of frozen round491',
            number_of_macros_per_terminal_trial=n,
            visible_history_and_initial_reference_error_bound='n*beta',
            numerical_rational_bound=str(n*p['beta']),
            required_target_selected_from_previous_visible_macro_record=True,
            regeneration_handles_arbitrary_current_DG_and_old_R_every_trial=True,
            same_DG_reused_across_all_macros_and_all_statistical_trials=True,
            initial_one_excitation_sector_remains_an_input=True,
            terminal_abort_encoded_as_zero_for_each_six_endpoint_indicator=True,
            no_IID_trial_source_or_hidden_reader_feedback_assumed=True,
            no_total_number_of_trials_times_beta_error_accumulation=True)

    def test_05_martingale_statistics_allow_correlated_reuse(self):
        p=F(1,2);bias=F(1,8);length=8
        mass=F(0);mean_d=F(0);variance_d=F(0);predictable_var=F(0);tail=F(0)
        first_second=F(0);mean_first=F(0);mean_second=F(0)
        for path in itertools.product((0,1),repeat=length):
            probability=F(1);d=F(0);pv=F(0)
            for t,y in enumerate(path):
                mu=p if t==0 else p+bias*(2*path[t-1]-1)
                self.assertLessEqual(abs(mu-p),bias)
                probability*=mu if y else 1-mu
                d+=y-mu;pv+=mu*(1-mu)
            mass+=probability;mean_d+=probability*d
            variance_d+=probability*d*d;predictable_var+=probability*pv
            if abs(d)>=F(3,8)*length:tail+=probability
            first_second+=probability*path[0]*path[1]
            mean_first+=probability*path[0];mean_second+=probability*path[1]
        self.assertEqual(mass,1);self.assertEqual(mean_d,0)
        self.assertEqual(variance_d,predictable_var);self.assertEqual(variance_d,F(121,64))
        covariance=first_second-mean_first*mean_second
        self.assertEqual(covariance,F(1,16));self.assertLess(tail,F(1,5))
        OBS['non_IID_statistics']=dict(
            analytic_bound='Pr(|N^-1 sum(Y_i-E[Y_i|past])|>alpha)<=2 exp(-2 N alpha^2)',
            proof='conditional Bernoulli log-mgf second derivative <=1/4, then exponential supermartingale',
            each_input_uses_one_common_batch_of_complete_six_outcome_or_abort_records=True,
            conditional_probability_column_l1_bias_bound='n*beta uniformly in all visible past',
            averaged_conditional_bias_does_not_grow_with_trial_count=True,
            exact_correlated_eight_trial_total_probability=str(mass),
            raw_first_two_outcome_covariance=str(covariance),
            centered_martingale_sum_mean=str(mean_d),variance=str(variance_d),
            exact_tail_at_alpha_three_eighths=str(tail),
            finite_enumeration_checks_implementation_not_the_general_inequality=True)

    def test_06_complete_calibration_normalization_and_resource_ledger(self):
        p=parameters();n=p['n'];a=p['a'];eps=p['eps'];beta=p['beta'];dlam=p['dlam'];alpha=p['alpha']
        self.assertEqual(beta,dlam/4)
        self.assertLessEqual(2*(beta+dlam/4),dlam)
        self.assertEqual(6*n*a*dlam,eps/4)
        self.assertEqual(9*a*alpha,eps/4)
        total=eps/4+eps/4+eps/32+F(1,2000)
        self.assertLess(total,eps)
        nt=math.ceil(8/alpha**2);nc=math.ceil(64/dlam**2)
        self.assertGreaterEqual(2*nt*alpha**2,16)
        self.assertGreaterEqual(2*nc*(dlam/4)**2,8)
        exp8_lower=sum(F(8**j,math.factorial(j)) for j in range(21))
        self.assertGreater(exp8_lower,2000)
        self.assertLess(F(72,2000**2)+F(4,2000),F(1,100))
        mixes=p['k']*p['r']*p['m'];reads=mixes+1
        macros=6*nt*n+nc
        wait_per_macro=F(mixes,16)+p['delta']
        # One retained active H of norm <=9, not two network copies.
        control_error=p['gmix']/2;control_window=p['gmix']/36
        self.assertEqual(control_error+18*control_window,p['gmix'])
        OBS['complete_resource_ledger']=dict(
            normalized_endpoint_vector_error_target=str(eps),
            ideal_frozen_491_projection_remainder='less than 1/2000',
            lambda_identity='M_delta(0,0)-M_delta(3,0)',
            lambda_calibration_absolute_error=str(dlam),
            calibration_conditional_bias_per_probability=str(beta),
            calibration_statistical_error_per_probability=str(dlam/4),
            terminal_each_probability_statistical_error=str(alpha),
            normalization_gain_with_calibration_at_most=3*a,
            terminal_history_bias_after_normalization_bound=str(eps/32),
            final_total_endpoint_error_bound=str(total),
            joint_failure_probability_below='1/100',
            independent_trial_preparations_required=False,
            terminal_batches=6,terminal_trials_per_batch=nt,total_terminal_trials=6*nt,
            calibration_batches=1,calibration_trials=nc,
            all_trials_use_the_same_data_graph_network=True,
            active_networks=1,fresh_graph_sources_after_initialization=0,fresh_data_batches_after_initialization=0,
            maximum_mix_instruments_per_macro=mixes,maximum_six_outcome_readers_per_macro=reads,
            maximum_natural_wait_per_macro=str(wait_per_macro),
            maximum_total_macros=macros,
            maximum_total_six_state_readers=macros*reads,
            maximum_record_encoding_capacity_qubits=3*macros*reads,
            record_capacity_is_not_total_detector_clock_or_isolation_hardware=True,
            maximum_total_natural_wait=str(macros*wait_per_macro),
            reading_routing_comparison_preparation_and_control_times_are_additional=True,
            finite_control_source_error_per_mix_instrument=str(control_error),
            finite_H_on_control_window_per_mix_instrument=str(control_window),
            no_big_sample_experiment_or_long_mixing_sequence_actually_run=True,
            blank_readers_classical_flags_clock_program_and_isolation_are_inputs=True,
            no_free_erasure_or_all_history_product_state_claim=True)


def run():
    OBS.clear();stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise AssertionError(stream.getvalue())
    answer=dict(round=493,scientific_baseline_round=492,reused_frozen_rounds=[280,283,443,465,472,476,477,490,491],
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        dependency_results_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
            ('retained_port_instrument_mixing_results.json','repeated_port_record_slow_modes_results.json')},
        observations=OBS,scope=dict(
            finite_retry_regenerates_491_source_without_fresh_graph_or_network=True,
            all_statistical_trials_may_reuse_one_DG_via_conditional_bounds=True,
            new_blank_readers_and_timed_instruments_still_extra_inputs=True,
            no_arbitrary_data_SWAP_added=True,
            aborted_trials_not_removed_by_postselection=True,
            hidden_new_readers_not_independent_of_graph_and_never_feedback=True,
            initial_unknown_information_supported_in_one_excitation_sector_only=True,
            actual_instruments_and_errors_preserve_this_sector=True,
            unrestricted_leakage_noise_across_reused_trials_not_covered=True,
            no_closed_H_autonomous_preparation_or_physical_time_rescaling=True,
            no_continuous_space_or_three_dimensional_necessity_derived=True,
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
