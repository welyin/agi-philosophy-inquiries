"""Reuse frozen 1250 enclosures with exact Fraction arithmetic; no new science group.
Default is read-only comparison. --write exclusively creates replay_results.json.
The original 1250 default replay is recorded separately in original_check_execution.json.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import argparse, hashlib, json
HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/'research_cognition_physics').is_dir() and (p/'scripts').is_dir())
ARCH=ROOT/'research_cognition_physics/archive_1234_'
OLD=ARCH/'1250'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def interval(v):
    a,b=map(F,v)
    assert a<=b
    return a,b
def add(x,y):return x[0]+y[0],x[1]+y[1]
def sub(x,y):return x[0]-y[1],x[1]-y[0]
def amp(x):return max(abs(x[0]),abs(x[1]))
def exact(x):return [str(t) for t in x]
def dec(x):
    with localcontext() as c:
        c.prec=52
        return str(Decimal(x.numerator)/Decimal(x.denominator))

def calculate():
    saved=json.loads((OLD/'results.json').read_text('utf-8'))
    independent=json.loads((OLD/'independent_results.json').read_text('utf-8'))
    execution=json.loads((HERE/'original_check_execution.json').read_text('utf-8'))
    assert execution['exit_code']==0 and execution['status']=='PASS_READ_ONLY_REPLAY'
    prior=json.loads(execution['stdout'])
    assert prior['mode']=='READ_ONLY_MATCH' and prior['status']=='PASS'
    for rel,digest in execution['frozen_files_sha256_before_and_after'].items():
        assert sha(ROOT/rel)==digest,rel
    for rel,digest in saved['source_sha256'].items():
        assert sha(ARCH/rel)==digest,rel
    for rel,digest in independent['bindings'].items():
        assert sha(OLD/rel)==digest,rel
    eta=F(saved['fixed']['error_per_menu'])
    gamma=F(saved['group_finite_balance']['raw_common_acceptance'])
    rows=[]
    for branch in saved['group_same_history_null']['branches']:
        q=branch['screen']; c=branch['read_contract']
        a=interval(q['geometric_area_delta']); B=interval(q['boost_full'])
        eg=add(interval(q['weight_correction']),interval(q['expansion_correction']))
        d=interval(q['positive_memory_gap'])
        assert eg==interval(q['geometry_correction'])
        assert a==sub((F(1),F(1)),interval(q['area_ratio_past']))
        raw=sub(a,B); corrected=sub(raw,eg)
        raw_omitted=sub(a,sub(B,d)); omitted=sub(raw_omitted,eg)
        assert raw==interval(q['raw_balance'])
        assert raw_omitted==interval(q['raw_omitted_balance'])
        assert corrected==interval(q['corrected_balance'])
        assert omitted==interval(q['corrected_omitted_balance'])
        assert sub(a,interval(q['area_from_metric_optics']))==interval(q['optical_residual'])
        assert corrected[0]<=0<=corrected[1]
        assert F(c['full_corrected_error'])==3*eta
        assert F(c['omitted_corrected_error'])==4*eta
        assert d[0]>7*eta
        legacy_margin=sub(d,(7*eta,7*eta))
        assert legacy_margin==interval(c['corrected_separation_margin'])
        # Full and omitted prediction enclosures expanded by their read contracts.
        full_with_contract=(corrected[0]-3*eta,corrected[1]+3*eta)
        omitted_with_contract=(omitted[0]-4*eta,omitted[1]+4*eta)
        gap=omitted_with_contract[0]-full_with_contract[1]
        assert gap>0
        raw_slack=gamma-(amp(raw_omitted)+3*eta)
        assert raw_slack>0 # Old uncorrected tolerance accepts the omitted ledger.
        assert interval(q['comoving_radius_tube'])[1]<5
        assert interval(q['affine_length'])[0]>0
        assert q['terminal_expansion_exact']=='0' and q['terminal_kt_exact']=='1'
        rows.append({'alpha':branch['alpha'],
            'reconstructed_corrected_balance':exact(corrected),
            'corrected_balance_absolute_upper':str(amp(corrected)),
            'positive_memory_gap':exact(d),
            'legacy_margin_after_7eta':exact(legacy_margin),
            'full_interval_with_original_contract':exact(full_with_contract),
            'omitted_interval_with_original_contract':exact(omitted_with_contract),
            'strict_disjoint_interval_gap':str(gap),
            'strict_disjoint_interval_gap_decimal':dec(gap),
            'additional_equal_error_threshold':str(gap/2),
            'additional_equal_error_threshold_decimal':dec(gap/2),
            'uncorrected_omission_acceptance_slack':str(raw_slack)})
    bounds=independent['symbolic_group']['independent_positive_bounds']
    old_lower=3*F(1,100)*F(1,27)*F(bounds['area_ratio_lower'])*F(bounds['affine_length_lower'])**2/2
    assert old_lower==F(bounds['normalized_memory_deletion_lower'])
    assert old_lower-7*eta==F(bounds['independent_corrected_read_margin_lower'])>0
    global_gap=min(F(r['strict_disjoint_interval_gap']) for r in rows)
    threshold=global_gap/2
    conservative=F('0.000000076464')
    assert conservative<threshold
    paths=[OLD/'proof.md',OLD/'check.py',OLD/'results.json',OLD/'independent_review.md',OLD/'independent_check.py',OLD/'independent_results.json',
           ARCH/'1242/check.py',ARCH/'1247/results.json',
           ARCH/'_admission/jacobson_common_realization_screen_after1278.md',
           ARCH/'_shared/ct_hg_adoption_after1278.md',HERE/'original_check_execution.json',Path(__file__)]
    return {'schema':'jacobson_1250_frozen_replay_v1','status':'PASS_WITH_SCOPE',
        'new_scientific_rounds':0,'new_scientific_groups':0,'inherited_round':1250,
        'method':'Exact Fraction reconstruction from frozen enclosures; original checker actually replayed read-only; no new model, ODE or instrument.',
        'fixed':saved['fixed'],'branches':rows,
        'all_branch_separation_gap':str(global_gap),
        'additional_equal_error_threshold':str(threshold),
        'conservative_sufficient_equal_error_bound':str(conservative),
        'hg_matching_extra_error':{'assigned':False,'value':None,
            'requirement':'e_full + e_omitted < all_branch_separation_gap; each e < threshold is sufficient.',
            'scope':'All additional HG/reference/clock/source/endpoint/instrument errors in the same normalized residual; not assumed zero.'},
        'old_independent_analytic_lower':str(old_lower),
        'old_independent_analytic_margin':str(old_lower-7*eta),
        'raw_omitted_ledger_remains_accepted':True,
        'effective_read_contract_realized_here':False,
        'source_bindings':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in paths}}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=calculate();out=HERE/'replay_results.json'
    if args.write:
        with out.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        mode='CREATED_ONCE'
    else:
        assert json.loads(out.read_text('utf-8'))==result,'saved replay differs'
        mode='READ_ONLY_MATCH'
    print(json.dumps({'status':result['status'],'mode':mode,'new_scientific_groups':0,
       'min_separation_gap':result['all_branch_separation_gap'],
       'equal_extra_error_threshold':result['additional_equal_error_threshold'],
       'HG_error_assigned':False},ensure_ascii=False))
