"""1212: native finite-reservoir work-value bounds; default replay is read-only.

All signed claims use exact Fraction bounds loaded from accepted sources.
NumPy tests the reused quantum inequality only; it does not certify the signs.
No author/old script is imported. --write creates a new result, never overwrites.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
RESEARCH = STAGE.parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def fs(x):
    return str(x.numerator) if x.denominator == 1 else str(x)

def interval_product(a, b):
    v = [x*y for x in a for y in b]
    return min(v), max(v)

def sum_intervals(xs):
    return sum((x[0] for x in xs), F(0)), sum((x[1] for x in xs), F(0))

def initial_drift(r, tau, tau_error, roots):
    er = (sum(p*e[0] for p,e in zip(r,roots)),
          sum(p*e[1] for p,e in zip(r,roots)))
    et = sum_intervals([interval_product(e, (p-tau_error,p+tau_error))
                        for p,e in zip(tau,roots)])
    return er[0]-et[1], er[1]-et[0]

def down(x, places=12):
    scale = 10**places
    return F(x.numerator*scale//x.denominator, scale)

def record(x):
    if isinstance(x, F): return fs(x)
    if isinstance(x, dict): return {k:record(v) for k,v in x.items()}
    if isinstance(x, (tuple,list)): return [record(v) for v in x]
    return x

def quantum_diagnostic(energies):
    rng = np.random.default_rng(1212003)
    width = float(max(energies)-min(energies))
    smallest = 1.
    cases = 0
    for beta in (1/3, 6/5):
        tau = np.exp(-beta*np.array(energies)); tau /= tau.sum()
        et = float(np.dot(energies,tau))
        for _ in range(24):
            z = rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
            rho = z@z.conj().T; rho /= np.trace(rho)
            ev = np.linalg.eigvalsh(rho)
            drel = float(sum(x*math.log(x) for x in ev if x>0)
                         -np.dot(np.diag(rho).real,np.log(tau)))
            de = float(np.dot(energies,np.diag(rho).real))-et
            quadratic = 2*(de/width)**2
            assert drel+1e-12 >= quadratic
            smallest = min(smallest,drel-quadratic); cases += 1
        assert abs(float(np.dot(tau,np.log(tau)-np.log(tau)))) < 1e-15
    return dict(cases=cases, minimum_slack=smallest,
                coherent_states_included=True,
                numerical_diagnostics_not_sign_certificate=True)

def run():
    p2=json.loads((STAGE/'1202/results.json').read_text(encoding='utf-8'))
    p4=json.loads((STAGE/'1204/results.json').read_text(encoding='utf-8'))
    p5=json.loads((STAGE/'1205/results.json').read_text(encoding='utf-8'))
    den=int(p2['preparation']['root_denominator'])
    roots=[(F(a,den),F(b,den)) for a,b in p2['preparation']['root_boxes']]
    width=F('1.094')
    assert roots[-1][1]-roots[0][0] < width
    cold=list(map(F,p2['one_marginal_rational']))
    hot=list(map(F,p4['new_reference']['one_marginal_rational']))
    assert sum(cold)==sum(hot)==1 and min(cold+hot)>0
    drift={}
    for name,r in [('cold',cold),('hot',hot)]:
        d=p4['temperature_comparisons'][name]
        tau=list(map(F,d['gibbs_population_rational']))
        assert abs(sum(tau)-1) <= 4*F(d['tau_entry_error_upper'])
        drift[name]=initial_drift(r,tau,F(d['tau_entry_error_upper']),roots)
    assert F('-0.003124') < drift['hot'][0] < drift['hot'][1] < F('-0.003123')
    assert F('0.000016') < drift['cold'][0] < drift['cold'][1] < F('0.000017')
    beta_c,beta_h=F(6,5),F(1,3)
    carnot=1-beta_h/beta_c
    assert carnot==F(13,18)
    d0=F(p4['multi_reference_thermodynamics']['cold_hot_initial_relative_entropy_defect_upper'])
    eta=F(9,10**15)
    assert F(p5['normal_clock']['conditional_joint_trace_error_upper']) < eta
    # Per initial 4-state marginal: |Delta S| <= eta*(1+47+2),
    # because 2^-47 <= eta, h2(eta)<=eta*(1+47), and log(3)<2.
    assert F(1,2**47) <= eta
    d0_transport=(beta_c+beta_h)*width*eta+100*eta
    d0_extra=F(2,10**12)
    assert d0_transport < d0_extra

    rows=[]
    for family in ('1204_independent_active_contact','1205_conserved_S_contact'):
        for t in (300,600):
            old=p4['actual_flow_predictions'][f't{t}_lambda1']
            if family.startswith('1204'):
                radius=F(old['single_energy_error_upper'])
                ac,ah=drift['cold'],drift['hot']; din=d0
                independence=True
            else:
                radius=F(p5['native_flow_transport'][f't{t}_lambda1']['new_actual_energy_radius_upper'])
                ac=tuple(x+s*width*eta for x,s in zip(drift['cold'],(-1,1)))
                ah=tuple(x+s*width*eta for x,s in zip(drift['hot'],(-1,1)))
                din=d0+d0_extra
                independence=False
            qc=F(old['delta_EA_center'])
            qh=-F(old['delta_EB_center'])
            qci=(qc-radius,qc+radius); qhi=(qh-radius,qh+radius)
            # In all four certified boxes the hot side lies below tau_h,
            # and the cold side lies above tau_c.
            hot_distance=qhi[0]-ah[1]
            cold_distance=qci[0]+ac[0]
            assert hot_distance>0 and cold_distance>0
            jh=2*(hot_distance/width)**2
            jc=2*(cold_distance/width)**2
            gap_hot=(jh-din)/beta_c
            gap_pair=(jh+jc-din)/beta_c
            assert gap_pair>gap_hot>0
            if t==300: assert gap_hot>F('0.00088') and gap_pair>F('0.0019')
            else: assert gap_hot>F('0.0076') and gap_pair>F('0.015')
            rows.append(dict(family=family,t=t,
                energy_error_radius=radius,hot_heat_interval=qhi,cold_gain_interval=qci,
                initial_hot_drift_interval=ah,initial_cold_drift_interval=ac,
                initial_relative_entropy_sum_upper=din,
                hot_relative_entropy_lower=jh,cold_relative_entropy_lower=jc,
                hot_only_work_value_gap_lower=gap_hot,
                hot_only_gap_display_down=down(gap_hot),
                two_endpoint_work_value_gap_lower=gap_pair,
                two_endpoint_gap_display_down=down(gap_pair),
                formal_Carnot_work_value_upper=carnot*qhi[1],
                cyclic_work_value_upper=carnot*qhi[1]-gap_hot,
                resource_credit_needed_to_reach_formal_Carnot=gap_hot,
                independence_of_full_machine_assumed=independence,
                finite_heat_window_is_required_for_any_completion=True))

    cycle=[]
    p_success_lower=F(499,500)  # 1205 proof section 2, strict lower bound.
    for t in (300,600):
        old=p4['actual_flow_predictions'][f't{t}_lambda1']
        new=p5['native_flow_transport'][f't{t}_lambda1']
        k0=F(old['A_pointer_probabilities_center'][0])
        kerr=F(new['pointer_probability_error_upper'])
        conditioned_nonblank=1-k0-kerr
        raw_nonblank=p_success_lower*conditioned_nonblank
        assert raw_nonblank>F('0.625') if t==300 else raw_nonblank>F('0.653')
        shape=F(new['rectangular_shape_error'])
        # Full 1000-material energy width, at both sections; no AB-only bound
        # for the preparation/stock tail contribution.
        radius_clock=2*F(1094)*eta+2*width*shape+F(1,10**10)
        tail=F(p5['normal_clock']['endpoint_total_interaction_upper'])
        radius_clock+=2*tail/p_success_lower
        material_work=F(old['bare_work_center'])
        clock_interval=(-material_work-radius_clock,-material_work+radius_clock)
        entropy_bound=F(p5['normal_clock']['conditional_clock_entropy_upper'])
        f_interval=(clock_interval[0]-entropy_bound/beta_c,
                    clock_interval[1]+entropy_bound/beta_c)
        if t==300: assert f_interval[1]<F('-0.000033')
        else: assert f_interval[0]<0<f_interval[1]
        cycle.append(dict(t=t,success_probability_strict_lower=p_success_lower,
            conditioned_K_nonblank_lower=conditioned_nonblank,
            raw_K_nonblank_lower=raw_nonblank,
            raw_machine_trace_distance_from_initial_strict_lower=raw_nonblank,
            raw_machine_return_distance_display_down=down(raw_nonblank),
            no_success_read_or_postselection_required=True,
            failure_contribution_kept_nonnegative=True,
            contact_clock_energy_change_interval=clock_interval,
            contact_clock_work_value_change_interval=f_interval,
            clock_gain_strictly_excluded=(t==300),
            full_material_width_used=1094,
            clock_return_and_future_uncopy_not_inferred=True))

    sources=[
        'archive_1201_/1202/results.json','archive_1201_/1204/results.json',
        'archive_1201_/1204/proof.md','archive_1201_/1205/results.json',
        'archive_1201_/1205/proof.md','archive_1201_/1209/proof.md',
        'archive_1201_/1210/acceptance.json','archive_1201_/1211/acceptance.json',
        'archive_629_652/research_note_638.md',
        'archive_702_741/research_note_714.md',
        'archive_1009_1043/research_note_1020.md']
    return record(dict(schema='native_finite_reservoir_cycle_1212_v1',
        status='AUTHOR_CHECKS_PASS_NOT_ACCEPTANCE',
        source_base='research_cognition_physics',
        source_sha256={s:sha(RESEARCH/s) for s in sources},
        input_certificates=dict(beta_c=beta_c,beta_h=beta_h,
            formal_Carnot_coefficient=carnot,energy_width_strict_upper=width,
            initial_mean_energy_drift_intervals=drift,
            initial_defect_upper=d0,
            prepared_marginal_trace_budget=eta,
            initial_defect_transport_derived_upper=d0_transport,
            initial_defect_transport_used_upper=d0_extra),
        finite_flow_bounds=rows,
        actual_factory_cycle_obstruction=cycle,
        reused_quantum_bound_diagnostic=quantum_diagnostic(
            [float((a+b)/2) for a,b in roots]),
        scope=dict(new_axiom=False,
            mature_entropy_and_finite_reservoir_tools_reused=True,
            full_energy_and_entropy_ledger_required=True,
            clock_Gibbs_partition_function_required=False,
            clock_free_energy_is_work_value_not_achievable_work=True,
            original_QQ_heat_engine_implemented=False,
            all_cycle_completions_excluded=False,
            reverse_Kelvin_to_Clausius_equivalence_proved=False,
            unconditional_four_laws_derived=False,
            stage_complete=False),
        check_groups=6))

def compare(a,b,path='root'):
    if isinstance(a,float) and isinstance(b,(int,float)):
        assert math.isclose(a,b,rel_tol=1e-10,abs_tol=2e-12),(path,a,b)
    elif isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a: compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+f'[{i}]')
    else: assert a==b,(path,a,b)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run(); target=HERE/'results.json'
    if args.write:
        assert not target.exists(),'Do not overwrite saved results'
        target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        compare(result,json.loads(target.read_text(encoding='utf-8')))
    print(json.dumps(dict(status='PASS',all_checks_passed=True,results_match=True,
                          check_groups=result['check_groups'],
                          mode='author_write' if args.write else 'read_only_replay')))
if __name__=='__main__': main()
