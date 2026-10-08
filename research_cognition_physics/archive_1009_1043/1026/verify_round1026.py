"""Round1026: reproducibility and frozen-scope checks, no live navigation."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import re
from urllib.parse import unquote

import numpy as np
import limiting_speed_selection as science

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
ROOT=BASE.parent
NOTE=HERE.parent/'research_note_1026.md'
OUT=HERE/'research_round_1026_checks.json'
OWN=['limiting_speed_selection.py','limiting_speed_selection_results.json','review.md',
     'selection_audit.md','input_dependency_update_v0_15.md','NEXT.md','verify_round1026.py']
HISTORICAL={
    'archive_1009_/research_note_1025.md','archive_1009_/1025/input_dependency_update_v0_14.md',
    'archive_1009_/1025/NEXT.md','archive_1009_/1009/input_dependency_ledger_v0_1.md',
    *[f'archive_301_341/research_note_{n}.md' for n in (334,335,339,341)],
    'archive_342_369/research_note_352.md','archive_923_934/research_note_924.md',
    'archive_956_989/research_note_957.md'}
LINK=re.compile(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)|^\[[^\]\n]+\]:\s*(\S+)\s*$',re.M)
EXCLUDED=re.compile(r'\\\[.*?\\\]|\$\$.*?\$\$|^```.*?^```\s*$',re.S|re.M)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def local_links(text):
    for match in LINK.finditer(EXCLUDED.sub(lambda match:' '*len(match[0]),text)):
        value=(match[1] if match[1] is not None else match[2]).strip()
        target=value[1:value.index('>')] if value.startswith('<') and '>' in value else re.split(r'\s+[\"\']',value,1)[0]
        local=unquote(target.split('#',1)[0])
        if local and not re.match(r'^[a-zA-Z]+:',local) and not local.startswith('//'):
            yield target,local


def independent_constants(a0,u0,eta_a=0,eta_u=0):
    m,M=min(a0,1),max(a0,1)
    def ff(a):return (3*(a+1)**3+4*a)/(3*a*(a+1)**2)
    def bb(a):return 2+4*(2*a+1)/(a*(a+1)**2)
    return ff(M)-eta_a*u0,ff(m)+eta_a*u0,bb(M)-eta_u*u0,bb(m)+eta_u*u0


def check_rows(case,eta_a=0,eta_u=0):
    a0,u0=case['a0'],case['u0']
    fm,fp,bm,bp=independent_constants(a0,u0,eta_a,eta_u)
    assert min(fm,bm)>0
    assert [row['ell'] for row in case['rows']]==[0,1,10,100,1000]
    for row in case['rows']:
        ell=row['ell'];delta=abs(row['a']-1)
        ulo,uhi=u0/(1+bp*u0*ell),u0/(1+bm*u0*ell)
        dlo=abs(a0-1)*math.exp(-fp/bm*math.log1p(bm*u0*ell))
        dhi=abs(a0-1)*math.exp(-fm/bp*math.log1p(bp*u0*ell))
        assert ulo-1e-12<=row['u']<=uhi+1e-12
        assert dlo-1e-10<=delta<=dhi+1e-10
        for key,value in [('u_lower',ulo),('u_upper',uhi),('delta_lower',dlo),('delta_upper',dhi)]:
            assert math.isclose(row['bounds'][key],value,rel_tol=1e-12,abs_tol=1e-15)
        if a0!=1:
            assert (row['a']-1)*(a0-1)>0
        if u0==0:
            assert row['a']==a0 and row['u']==0


def verify(prospective=False):
    fresh=science.run()
    science.compare(fresh,json.loads(science.OUT.read_text('utf8')))
    assert fresh['round']==1026 and fresh['all_scientific_calibrations_passed']
    assert fresh['new_calibration_groups']==1 and fresh['cumulative_test_groups']==3803
    assert fresh['new_cognitive_axioms']==0
    for key in ('one_loop_beta_functions_rederived_from_loop_integrals','numerical_error_rigorously_certified',
                'actual_higher_loop_envelope_certified','physical_common_metric_or_all_Lorentz_operators_derived',
                'actual_instruments_or_initial_conditions_generated','fixed_line_linearization_is_exact_off_line',
                'exact_finite_scale_equalization_claimed','infinite_IR_extrapolation_is_required_for_current_effective_goal'):
        assert fresh[key] is False,key
    assert fresh['one_loop_flow_is_adopted'] and fresh['mathematical_bounds_proved_in_note_not_by_sampling']
    assert fresh['primary_source']=='https://arxiv.org/html/1102.0789' and fresh['primary_equations']==[18,19]
    assert fresh['near_line_exponent']=='7/15'
    assert len(fresh['exact_chain_rule_checks'])==4
    for row in fresh['exact_chain_rule_checks']:
        cf,cb,v=map(F,(row['cf'],row['cb'],row['v']))
        a,u=cb/cf,v/cf**3
        f=(3*(a+1)**3+4*a)/(3*a*(a+1)**2)
        b=2+4*(2*a+1)/(a*(a+1)**2)
        dotcf=4*v*(cb-cf)/(3*cb*(cf+cb)**2)
        dotcb=v*(cf*cf-cb*cb)/(cb*cf**3)
        dotv=-2*v*v*(cb**3+2*cb**2*cf+3*cb*cf**2+4*cf**3)/(cb*cf**3*(cf+cb)**2)
        assert dotcb/cf-a*dotcf/cf == -u*f*(a-1)
        assert dotv/cf**3-3*u*dotcf/cf == -b*u*u
        assert row['residuals']==['0','0','0']
    cases=fresh['finite_trajectories']
    assert len(cases)==12 and sum(c['u0']==0 for c in cases)==2
    assert sum(c['a0']==1 for c in cases)==2
    for case in cases:
        check_rows(case)
        assert case['original_reduced_max_error']<2e-9
        assert case['refinement_max_differences'][-1]<3e-9
        assert case['speed_monotonicity_checked'] and not case['refinement_difference_is_rigorous_error_bound']
        sign=np.sign(case['a0']-1)
        cf=np.array([r['cf'] for r in case['rows']]);cb=np.array([r['cb'] for r in case['rows']])
        assert np.all(sign*np.diff(cf)>=-1e-12) and np.all(sign*np.diff(cb)<=1e-12)
        for row in case['rows']:
            u=row['g']**2/(8*math.pi**2*row['cf']**3)
            assert abs(u-row['u'])<2e-10
            assert row['speed_integral_error']<1e-9 and row['gauss_refinement_error']<2e-12
            if case['a0']==1:
                assert row['a']==1 and abs(row['u']-case['u0']/(1+5*case['u0']*row['ell']))<2e-10
    witness=next(c for c in cases if c['a0']==1.1 and c['u0']==.02)
    assert .0114<witness['rows'][-1]['absolute_ratio_difference']<.0115
    near=fresh['near_line_checks'];assert len(near)==4
    assert [r['delta0'] for r in near]==[-.001,-.0005,.0005,.001]
    for row in near:
        assert .06<row['max_error_over_delta_squared']<.08
        baseline=row['delta0']*(1+.1*np.array(row['ell']))**(-7/15)
        assert np.max(np.abs(baseline-np.array(row['linearized_delta'])))<1e-15
        assert row['max_error']>1e-8
    for full,half in ((near[0],near[1]),(near[3],near[2])):
        assert 3.95<full['max_error']/half['max_error']<4.05
    for case in fresh['hypothetical_remainder_checks']:
        check_rows(case,1,1)
        assert case['actual_QFT_remainder_certified'] is False
    assert len(fresh['hypothetical_remainder_checks'])==2
    threshold=fresh['threshold_checks'];assert len(threshold)==3
    for row in threshold:
        fm,fp,bm,bp=independent_constants(row['a0'],row['u0'])
        ratio=abs(row['a0']-1)/row['epsilon']
        low=(ratio**(bm/fp)-1)/(bm*row['u0'])
        high=(ratio**(bp/fm)-1)/(bp*row['u0'])
        assert math.isclose(low,row['necessary_ell'],rel_tol=1e-12)
        assert math.isclose(high,row['sufficient_ell'],rel_tol=1e-12)
        assert 0<low<=high and row['this_is_not_a_certified_physical_RG_range']
    assert fresh['maximum_original_reduced_error']<2e-9
    assert fresh['maximum_speed_integral_error']<1e-9
    assert set(fresh['historical_source_sha256'])==HISTORICAL
    for name,digest in fresh['historical_source_sha256'].items():assert sha(BASE/name)==digest,name
    assets=[NOTE]+[HERE/name for name in OWN]
    assert len(assets)==8
    assert not any(p.name in {'README.md','research_direction.md','RESEARCH_STATE.md'} for p in assets)
    links=0
    for path in assets:
        assert path.is_file(),path
        if path.suffix=='.md':
            text=path.read_text('utf-8-sig');assert text.count('$$')%2==0,path
            for target,local in local_links(text):
                resolved=(path.parent/local.replace('\\','/')).resolve()
                assert resolved.exists() or (prospective and resolved==OUT),(path,target)
                links+=1
    assert all(f'## {n}.' in NOTE.read_text('utf8') for n in range(1,11))
    assert '独立科学签审通过' in (HERE/'review.md').read_text('utf8')
    return dict(round=1026,date='2026-10-08',all_delivery_checks_passed=True,scientific_result_reproduced=True,
                new_calibration_groups=1,cumulative_research_groups=3803,new_cognitive_axioms=0,
                calibration_kind='adopted one-loop ODE and exact comparison-bound calibration',
                exact_chain_cases=4,finite_trajectories=12,finite_scale_rows=60,near_line_controls=4,
                hypothetical_remainder_controls=2,threshold_controls=3,local_links_checked=links,
                maximum_original_reduced_difference=fresh['maximum_original_reduced_error'],
                maximum_refinement_difference=fresh['maximum_refinement_difference'],
                maximum_speed_integral_difference=fresh['maximum_speed_integral_error'],
                numerical_ODE_error_certified=False,physical_high_loop_bound_certified=False,
                actual_finite_domain_or_all_Lorentz_structure_generated=False,
                live_navigation_frozen=False,neighboring_round_frozen=False,goal_completed=False,
                visual_checks_performed=False,frozen_current_files=8,historical_input_files=len(HISTORICAL),
                historical_source_sha256=fresh['historical_source_sha256'],
                source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in assets})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write','--write-receipt',dest='write',action='store_true')
    args=parser.parse_args();result=verify(prospective=args.write)
    if args.write:
        with OUT.open('x',encoding='utf8') as stream:stream.write(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    else:assert result==json.loads(OUT.read_text('utf8'))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith('sha256')},ensure_ascii=False,indent=2))
