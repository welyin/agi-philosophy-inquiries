"""Independent bounded York measure audit; default read-only, first --write exclusive.
No author-code imports. Fraction proof gates and Decimal composite Simpson diagnostics.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
INPUTS = ('README.md', 'proof.md', 'check.py', 'results.json', 'source_baseline.json')

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def atan_small(x):
    power=x
    total=x
    n=1
    while True:
        power *= -x*x
        term=power/Decimal(2*n+1)
        total += term
        if abs(term) < Decimal('1e-85'):
            return total
        n += 1

def simpson(fun, left, right, n):
    assert n > 0 and n % 2 == 0
    h=(right-left)/Decimal(n)
    acc=fun(left)+fun(right)
    for j in range(1,n):
        acc += (4 if j%2 else 2)*fun(left+h*j)
    return h*acc/3

def baseline_audit():
    stage=HERE.parent.parent
    root=next(p for p in HERE.parents if (p/'research_cognition_physics').is_dir())
    data=json.loads((HERE/'source_baseline.json').read_text(encoding='utf-8'))
    count=0
    assert sha(stage/'_shared/goal_prompt.txt') == data['goal_sha256']
    count += 1
    for name,digest in data['source_sha256'].items():
        assert sha(stage/name) == digest, name
        count += 1
    for name,digest in data['project_sources_sha256'].items():
        assert sha(root/name) == digest, name
        count += 1
    for name,digest in data['external_source_sha256'].items():
        assert sha(Path(name)) == digest, name
        count += 1
    return count

def calculate():
    # r=16/(5 sqrt(5) t^(3/2)); monotonicity fixes both endpoints.
    tlo,thi=F(99,500),F(101,500)
    r2lo=F(256,125)/thi**3
    r2hi=F(256,125)/tlo**3
    assert r2lo == F(256000000,101**3)
    assert r2hi == F(256000000,99**3)
    assert F(8)**2 < r2lo < r2hi < F(33,2)**2
    # f''=pi(8-r), hence -27<f''<0 without evaluating a square root.
    abs_f2=F(27)
    assert (F(33,2)-8)*F(22,7) < abs_f2
    half=F(1,10)
    exponent_drop=abs_f2*half**2/2
    assert exponent_drop == F(27,200)
    gmin=1-exponent_drop # exp(-x)>1-x for x>0
    assert gmin == F(173,200) and gmin > F(6,7)
    plo,phi=gmin/(1+gmin),1/(1+gmin)
    assert 1-2*plo-plo**2 < 0
    sharp_gap=phi*(1-phi)/(1+phi)
    assert sharp_gap == F(34600,213729)
    coarse_gap=F(21,130)
    assert sharp_gap > coarse_gap > F(4,25)
    conditional=coarse_gap-2*F(1,100)
    assert conditional == F(46,325) and conditional > F(7,50)
    # Independent Lipschitz-midpoint error proof, used only to compare author diagnostic.
    lipschitz=abs_f2*half
    eta_times_n=lipschitz*half**2/4
    denominator=2*half*gmin
    normalized_error=eta_times_n/denominator
    assert normalized_error == F(27,692)
    assert normalized_error < F(99,2404) < F(297,4808)
    # For any positive equal-count halves with gmin<g<=1, the finite lower bound
    # follows from the same rational inequalities; no continuum limit is required.
    exact={
        'r_squared_lower':str(r2lo),'r_squared_upper':str(r2hi),
        'r_upper': '33/2', 'f_second_derivative_strict_lower':'-27',
        'f_second_derivative_strict_upper':'0',
        'relative_exponent_strict_lower':str(-exponent_drop),
        'g_strict_lower':str(gmin), 'g_upper':'1',
        'p0_strict_lower':str(plo),'p0_strict_upper':str(phi),
        'gap_strict_lower':str(sharp_gap),
        'simpler_gap_strict_lower':str(coarse_gap),
        'gap_after_two_1pct_errors_strict_lower':str(conditional),
        'midpoint_probability_error_times_N':str(normalized_error),
        'reweighted_probability_error_times_N':str(2*normalized_error),
        'finite_equal_count_half_grids_obey_same_bound':True,
    }
    with localcontext() as ctx:
        ctx.prec=80
        pi=16*atan_small(Decimal(1)/5)-4*atan_small(Decimal(1)/239)
        root0=(Decimal(1)/5).sqrt()
        # Rationalize twice; this expression is distinct from the author's
        # linear-in-d subtraction and avoids large O(10^4) exponent terms.
        def weight(d):
            rt=(Decimal(1)/5-d/50).sqrt()
            exponent=pi*d*d*(4-(Decimal(32)/5)/(rt+root0)**2)
            return exponent.exp()
        rows=[]
        raw=[]
        for n in (64,128,256):
            low=simpson(weight,Decimal('-0.1'),Decimal(0),n)
            high=simpson(weight,Decimal(0),Decimal('0.1'),n)
            p0=high/(low+high)
            p1=2*high/(low+2*high)
            gap=p1-p0
            assert Decimal(sharp_gap.numerator)/sharp_gap.denominator < gap
            assert abs(p1-2*p0/(1+p0)) < Decimal('1e-70')
            rows.append({'panels_per_half':n,'p0':format(p0,'.35f'),
                         'p1':format(p1,'.35f'),'gap':format(gap,'.35f')})
            raw.append((p0,p1,gap))
        successive=max(abs(a-b) for a,b in zip(raw[-2],raw[-1]))
        assert successive < Decimal('1e-10')
        authors=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
        assert authors['status']=='PASS_BOUNDED_REFERENCE_MEASURE_SCREEN'
        last=authors['diagnostics_not_proof'][-1]
        assert last['N_per_half']==512
        diagnostic_delta=abs(raw[-1][0]-Decimal(str(last['p0'])))
        assert diagnostic_delta < Decimal('1e-8')
        assert diagnostic_delta < Decimal(normalized_error.numerator)/(normalized_error.denominator*512)
        diagnostics={
            'method':'80-digit Decimal composite Simpson; Machin pi; double-rationalized exponent',
            'not_used_as_continuous_proof':True,
            'rows':rows,
            'last_two_resolution_max_difference':format(successive,'.8E'),
            'vs_author_512_midpoints_p0_difference':format(diagnostic_delta,'.8E'),
        }
    return {
        'verdict':'PASS_WITH_SCOPE','independent_bounded_suites':1,
        'new_formal_rounds':0,'new_formal_scientific_groups':0,'new_adopted_axioms':0,
        'exact':exact,'diagnostics_not_proof':diagnostics,
        'baseline_hashes_matched':baseline_audit(),
        'source_sha256':{name:sha(HERE/name) for name in INPUTS},
        'independent_code_sha256':sha(Path(__file__).resolve()),
        'scope':{
            'fixed_mass_window_not_full_ensemble_support':True,
            'two_explicit_reference_counting_inputs':True,
            'passive_unitary_transport_not_compared_as_new_physics':True,
            'area_factor_not_certified_microscopic_degeneracy':True,
            'EB_only_Gibbs_identity_certified':False,
            'actual_instrument_or_preparation_certified':False,
            'non_trace_class_implies_finite_prediction_failure':False,
            'original_fixed_measure_ambiguous':False,
            'GR_generated_or_CT_adopted':False,
            'full_goal_complete':False,
        },
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=calculate()
    path=HERE/'independent_results.json'
    if args.write:
        with path.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
        print('EXCLUSIVE_FIRST_WRITE PASS_WITH_SCOPE')
    else:
        prior=json.loads(path.read_text(encoding='utf-8'))
        assert prior == result, 'Independent read-only replay differs'
        print('READ_ONLY_MATCH PASS_WITH_SCOPE')
    print(json.dumps({'gap_lower':result['exact']['simpler_gap_strict_lower'],
                      'last_diagnostic':result['diagnostics_not_proof']['rows'][-1]},ensure_ascii=False))

if __name__=='__main__':
    main()
