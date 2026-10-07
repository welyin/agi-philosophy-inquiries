"""848: exact causal equation-test calibration, not the original spacetime PDE.

The original code is used for content differences. A finite mixed wave system
checks dual tests, causal integration, free-quotient order and source counting.
The original continuum statement uses 785/846 Green identities analytically.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,itertools,json,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent; TARGET=HERE/'causal_equation_observable_results.json'
sys.path.insert(0,str(HERE.parent/'829'))
import majorana_code_source_bridge as code

def run():
    gamma,_,_,_,comp,sx,sz=code.code_data()
    sy=[code.mul((0,0,3),code.mul(sz[b],sx[b])) for b in range(5)]
    neg=lambda a:code.mul((0,0,2),a)
    sources=[neg(sz[0]),neg(sz[1]),neg(sy[3])]
    logical=code.mul((0,0,1),code.product(gamma[i] for i in (0,1,4,5,12,14)))
    delta=F(27,50)
    def content(ops):
        kind,phase=comp(code.mul(logical,code.product(ops)))
        if kind!='scalar':return F(0)
        assert phase in (0,2)
        return delta*(1 if phase==0 else -1)
    assert all(content([a])==0 for a in sources)
    assert all(content([a,b])==0 for a in sources for b in sources)
    assert content(sources)==delta
    # Three coupled physical components, seven interior time equations.
    n=9; count=3
    c=[[F(1,5),F(1,9),F(0)],[F(1,9),F(1,4),F(1,11)],
       [F(0),F(1,11),F(1,6)]]
    profile=[[F(0)]*count for _ in range(n)]
    for t in range(1,n-1):
        for a in range(count):profile[t][a]=F((t+a)%4+1,5+a)
    zero=lambda:[[F(0) for _ in range(count)] for _ in range(count)]
    x=[zero() for _ in range(n)]
    for t in range(1,n-1):
        for a in range(count):
            for j in range(count):
                forcing=profile[t][a] if a==j else F(0)
                x[t+1][a][j]=2*x[t][a][j]-x[t-1][a][j]-sum(c[a][b]*x[t][b][j] for b in range(count))-forcing
    residuals=[]
    for t in range(1,n-1):
        for a in range(count):
            for j in range(count):
                residuals.append(x[t+1][a][j]-2*x[t][a][j]+x[t-1][a][j]
                    +sum(c[a][b]*x[t][b][j] for b in range(count))
                    +(profile[t][a] if a==j else F(0)))
    assert not any(residuals)
    tests=[]
    for a in range(count):
        h=[F(0)]*n
        weights=[F(1),F(2+a),F(1)]
        norm=sum(weights[t-3]*profile[t][a] for t in range(3,6))
        for t in range(3,6):h[t]=weights[t-3]/norm
        tests.append(h)
    def equation_tests(include_mixing):
        outputs=[]
        for a,h in enumerate(tests):
            coeff=[[F(0) for _ in range(count)] for _ in range(n)]
            for t in range(1,n-1):
                coeff[t-1][a]+=h[t];coeff[t][a]-=2*h[t];coeff[t+1][a]+=h[t]
                if include_mixing:
                    for b in range(count):coeff[t][b]+=h[t]*c[a][b]
            # This is the transposed equation acting on the compact test.
            outputs.append([sum(coeff[t][b]*x[t][b][j] for t in range(n) for b in range(count)) for j in range(count)])
        return outputs
    good=equation_tests(True);bad=equation_tests(False)
    assert good==[[-F(a==b) for b in range(count)] for a in range(count)]
    assert bad!=good
    def triple(rows):
        value=F(0)
        for indices in itertools.product(range(count),repeat=3):
            weight=F(1)
            for a,j in enumerate(indices):weight*=rows[a][j]
            value+=weight*content([sources[j] for j in indices])
        return value
    assert triple(good)==-delta and triple(bad)!=-delta
    # Nonorthogonal slice: physical field v*x, source cov*j + ell*extension.
    v=[F(2),F(-1)];cov=[F(-1),F(-3)];ell=[F(1),F(2)]
    assert sum(a*b for a,b in zip(v,cov))==1
    assert sum(a*b for a,b in zip(v,ell))==0
    extensions_checked=0
    for t in range(1,n-1):
        for a in range(count):
            extension=F((t+1)*(a+2),13)
            raw=[cov[k]*profile[t][a]+ell[k]*extension for k in range(2)]
            assert sum(v[k]*raw[k] for k in range(2))==profile[t][a]
            extensions_checked+=1
    # At hbar^3, with three roots of minimum field degree two, only (2,2,2)
    # and zero quantum/contraction order can retain six code legs.
    allowed=[]
    for degrees in itertools.product(range(2,7),repeat=3):
        for q in range(4):
            if sum(degrees)+2*q==6:allowed.append((degrees,q))
    assert allowed==[((2,2,2),0)]
    native=json.loads((HERE.parent/'847/native_triple_source_probe_results.json').read_text('utf-8'))
    return dict(round=848,all_checks_passed=True,fresh_test_groups=1,
        calibration_is_original_spacetime_PDE=False,exact_rational_causal_integration=True,
        time_sites=n,physical_components=count,equation_residual_nonzero_coefficients=0,
        field_equation_test_coefficients=[[str(a) for a in row] for row in good],
        mixed_equation_terms_cannot_be_deleted=True,
        triple_content_difference=str(triple(good)),
        omitted_mixed_equation_triple_difference=str(triple(bad)),
        premature_free_quotient_prediction='0',
        premature_free_quotient_error=str(delta),
        source_extension_tests=extensions_checked,
        leading_three_root_six_code_leg_choices=1,
        native_847_coefficients_reused_not_new_PDE_values=[
            dict(momentum=row['momentum'],joint_equation_test_cumulant=-row['normalized_triple_content_difference'])
            for row in native['coefficient_rows']],
        pure_metric_observable_claimed=False,one_point_geometry_difference_proven=False,
        finite_coupling_or_autonomous_detector_proven=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();out=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert out==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(out,ensure_ascii=False,indent=2))
