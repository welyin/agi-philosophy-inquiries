"""Independent 1206 all-raw-trial pointer audit; no new author imports.
Default recomputes read-only and emits the complete saved JSON. --write creates
independent_results.json once. Frozen 1204 kernels are reused as certificates;
this script does not pretend to rerun their matrix-exponential proof.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from decimal import Decimal,localcontext
from fractions import Fraction as F
from pathlib import Path
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent.parent

def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def matvec(a,b):return [dot(row,b) for row in a]
def decimal(x,places=24):
    with localcontext() as c:
        c.prec=places+20
        return str((Decimal(x.numerator)/Decimal(x.denominator)).quantize(Decimal(10)**-places))
def lo(x,places=24):
    scale=10**places
    y=x.numerator*scale//x.denominator
    return str(Decimal(y).scaleb(-places))
def hi(x,places=24):return lo(-x,places).replace('-','',1) if x<=0 else decimal_ceiling(x,places)
def decimal_ceiling(x,places=24):
    scale=10**places;y=-((-x.numerator*scale)//x.denominator)
    return str(Decimal(y).scaleb(-places))
def exp_interval(x,n=60):
    # Positive series. For the remainder, every subsequent ratio is at most
    # x/(n+2), so a geometric bound proves the enclosure.
    assert 0<=x<n+2
    s=sum((x**k/F(math.factorial(k)) for k in range(n+1)),F(0))
    first=x**(n+1)/F(math.factorial(n+1))
    return s,s+first/(1-x/F(n+2))
def run():
    old=load(STAGE/'1204/results.json')
    prior=load(STAGE/'1204/independent_results.json')
    prep=load(STAGE/'1205/results.json')
    dims=list(map(int,prep['fixed_program']['dimensions']))
    probs=[1-(1-F(d,2**k))**10 for d,k in zip(dims,(173,195))]
    p=probs[0]*probs[1]
    assert p==probs[1]*probs[0]
    assert F('.9990864202096577')<p<F('.9990864202096579')
    rC=list(map(F,prior['cold_shell']['one_rational']))
    rH=list(map(F,prior['hot_shell']['one_rational']))
    assert sum(rC)==sum(rH)==1
    den=prior['root_denominator']
    boxes=[[F(a,den),F(b,den)] for a,b in prior['independent_root_boxes']]
    e=[(a+b)/2 for a,b in boxes]
    Ec=dot(rC,e);Eh=dot(rH,e);w=F('1.094')
    assert boxes[-1][1]-boxes[0][0]<w
    weights={'CH':[a*b for a in rC for b in rH],
             'HC':[a*b for a in rH for b in rC]}
    eps=F(prep['normal_clock']['epsilon_clock_contract'])
    assert eps==F(4,10**15)
    coeff=F(1,10**10)
    sys=F(1,1000);stat=F(4,1000);m=2*w*(sys+stat)
    assert m==F('.01094')
    contacts={}
    for t in (600,300):
        row=old['full_native_contacts'][str(t)]
        K=[list(map(F,r)) for r in row['pointer_K']]
        P=[list(map(F,r)) for r in row['full_transition_P']]
        acoef=list(map(F,row['delta_EA_coefficients']))
        bcoef=list(map(F,row['delta_EB_coefficients']))
        assert F(row['each_saved_entry_error_upper'])==coeff
        # Pure source-exchange symmetry applies to the rectangular comparison.
        swap=lambda q:4*(q%4)+q//4
        swaperr=max(abs(P[a][b]-P[swap(a)][swap(b)]) for a in range(16) for b in range(16))
        assert swaperr<=2*coeff
        margerr=max(abs(K[a][i]-sum(P[4*a+b][i] for b in range(4))) for a in range(4) for i in range(16))
        assert margerr<=5*coeff
        raw={};means={};normerrors={}
        for direction in ('CH','HC'):
            success=matvec(K,weights[direction])
            raw[direction]=[p*x+(1-p)*int(j==0) for j,x in enumerate(success)]
            means[direction]=dot(e,raw[direction])
            normerrors[direction]=abs(sum(raw[direction])-1)
            assert normerrors[direction]<=4*coeff
            assert all(0<=x<=1 for x in raw[direction])
        proxy=(means['CH']-means['HC']-p*(Ec-Eh))/2
        q_cond=dot(weights['CH'],[(a-b)/2 for a,b in zip(acoef,bcoef)])
        q=p*q_cond
        # Two independently represented centers: K histogram versus certified
        # energy-change coefficients. This tolerance is a consistency check,
        # while the analytic identity transports the existing 1e-10 certificate.
        assert abs(proxy-q)<5*coeff
        qsave=F(old['actual_flow_predictions'][f't{t}_lambda1']['Q_dir_center'])
        assert q_cond==qsave
        reverseq=p*dot(weights['HC'],[(a-b)/2 for a,b in zip(acoef,bcoef)])
        assert abs(q+reverseq)<=2*coeff
        shape=F(1,10000*(600//t))
        actual_radius=w*(2*eps+shape)+coeff
        proxy_radius=w*(eps+shape)+coeff
        coarse_radius=actual_radius
        inference=actual_radius+proxy_radius+m
        coarse_inference=2*coarse_radius+m
        all_data_lower=q-actual_radius-2*proxy_radius-2*m
        coarse_all_data_lower=q-3*coarse_radius-2*m
        assert q-actual_radius>0 and all_data_lower>0 and coarse_all_data_lower>0
        contacts[str(t)]={
          'raw_K_nominal_CH':[decimal(x) for x in raw['CH']],
          'raw_K_nominal_HC':[decimal(x) for x in raw['HC']],
          'raw_K_each_probability_total_model_radius_upper':decimal_ceiling(eps+shape+coeff),
          'raw_K_nominal_normalization_residual_max':float(max(normerrors.values())),
          'pointer_numeric_mean_CH':decimal(means['CH']),
          'pointer_numeric_mean_HC':decimal(means['HC']),
          'proxy_from_K_nominal':decimal(proxy),
          'center_from_certified_delta_energy_coefficients':decimal(q),
          'K_vs_delta_energy_center_discrepancy':float(abs(proxy-q)),
          'reverse_center_antisymmetry_residual':float(abs(q+reverseq)),
          'saved_rectangular_swap_residual':float(swaperr),
          'K_full_P_marginal_residual':float(margerr),
          'actual_raw_Q_model_radius_upper':decimal_ceiling(actual_radius),
          'proxy_model_radius_upper':decimal_ceiling(proxy_radius),
          'actual_raw_Q_lower':lo(q-actual_radius),
          'actual_raw_Q_upper':decimal_ceiling(q+actual_radius),
          'data_to_actual_Q_inference_radius_upper':decimal_ceiling(inference),
          'coarse_inference_radius_upper':decimal_ceiling(coarse_inference),
          'all_budget_compliant_data_actual_CI_lower_floor':lo(all_data_lower),
          'coarse_all_data_CI_lower_floor':lo(coarse_all_data_lower),
          'shape_error_exact':str(shape),
          'success_conditioning_denominator_used':False,
          'actual_parent_AB_swap_symmetry_assumed':False}
    elo,ehi=exp_interval(F(8))
    failure_upper=32/elo
    failure_lower=32/ehi
    confidence_lower=1-failure_upper
    assert confidence_lower>F('.9892651959071196')
    # Removal controls: common failure output and equal p are essential.
    p1,p2=F(9,10),F(8,10)
    c,h,blank=F(1),F(3),F(0)
    no_flow_false_proxy=(p1*c+(1-p1)*blank-p2*h-(1-p2)*blank-p1*(c-h))/2
    assert no_flow_false_proxy==F(3,20)
    mismatch_blank_offset=(1-p)*F(1)/2
    assert mismatch_blank_offset>0
    sourcepaths=[STAGE/'1204/results.json',STAGE/'1204/independent_results.json',
      STAGE/'1204/acceptance.json',STAGE/'1205/proof.md',STAGE/'1205/results.json',
      STAGE/'1205/acceptance.json']
    return {
      'schema':'independent_1206_all_raw_pointer_v1',
      'status':'INDEPENDENT_CONDITIONAL_AUDIT_NOT_STAGE_ACCEPTANCE',
      'all_checks_passed':True,
      'reference_counts':{'dimensions':[str(x) for x in dims],
        'success_C':decimal(probs[0]),'success_H':decimal(probs[1]),
        'joint_success_exact':str(p),'joint_success':decimal(p),
        'same_probability_in_both_directions':True,
        'cold_mean_energy_from_root_midpoints':decimal(Ec),
        'hot_mean_energy_from_root_midpoints':decimal(Eh)},
      'analytic_scope':{
        'ideal_failure_control_QQ_and_copy_both_zero':True,
        'ideal_failure_free_material_evolution_preserves_hA_hB':True,
        'ideal_failure_pointer_state':'|0><0|',
        'failure_pointer_statistic_is_e0_not_zero':True,
        'pointer_bare_H_is_zero_but_statistic_encodes_material_h':True,
        'identity':'(mu_CH-mu_HC-p*(E_C-E_H))/2 = p*Q_success = Q_raw',
        'full_parent_Duhamel_used_without_success_denominator':True,
        'two_sections_in_actual_energy_change':True,
        'separate_proxy_and_actual_comparison_errors':True,
        'success_readout_instrument_required':False,
        'no_success_postselection_or_batch_fill':True},
      'contacts':contacts,
      'sampling':{'fixed_raw_trials_per_mode':250000,'modes':4,
        'raw_trials_total':1000000,'observed_K_categories_per_mode':4,
        'each_probability_sampling_error':str(stat),
        'each_probability_terminal_readout_error_contract':str(sys),
        'energy_proxy_measurement_radius':str(m),
        'Hoeffding_exponent_exact':str(2*250000*stat**2),
        'union_tail_count':32,
        'joint_statistical_failure_interval':[lo(failure_lower),decimal_ceiling(failure_upper)],
        'joint_statistical_confidence_lower':lo(confidence_lower),
        'rational_positive_series_terms_for_exp8':61,
        'all_budget_compliant_data_have_positive_actual_Q_CI':True,
        'this_is_not_already_observed_experimental_data':True,
        'fresh_independent_full_resources_per_raw_trial':True,
        'success_readout_error_added':False,
        'systematic_point001_is_assumed_K_readout_contract_not_Hoeffding_result':True,
        'material_resource_occurrences':1000*1000000,
        'seed_bit_resource_occurrences':3900*1000000},
      'removal_controls':{
        'unequal_success_probabilities_with_no_energy_exchange_fake_proxy':str(no_flow_false_proxy),
        'different_failure_pointer_statistics_offset':decimal(mismatch_blank_offset),
        'correlated_reused_pools_not_admitted_to_iid_Hoeffding':True},
      'source_sha256':{str(p.relative_to(STAGE)).replace('\\','/'):sha(p) for p in sourcepaths},
      'scope':{'new_author_code_imported':False,'historical_code_imported':False,
        'frozen_native_numeric_certificates_reused_not_rederived':True,
        'raw_histogram_reconstructed_independently':True,
        'continuous_clock_simulated':False,'new_physical_device_built':False,
        'new_B_energy_pointer_built':False,'universal_heat_flow_theorem_claimed':False}}
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();r=run();out=HERE/'independent_results.json'
    if args.write:
        if out.exists():raise FileExistsError('--write refuses overwrite')
        out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==load(out)
    print(json.dumps(r,ensure_ascii=False,indent=2))
