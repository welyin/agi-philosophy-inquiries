"""845: exact original-code coefficients and an explicitly auxiliary response.

The code is the original ten-CAR code. The scalar equation q+a*q^3+h*B=0
is a diagnostic of the order bound, not the original Einstein/matter PDE.
No historical certificate is recomputed or overwritten here.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,itertools,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'record_source_hbar_filtration_results.json'
sys.path.insert(0,str(HERE.parent/'829'))
import majorana_code_source_bridge as code

def run():
    gamma,_,_,_,comp,_,_=code.code_data()
    ident=code.I
    logical=code.mul((0,0,1),code.product(gamma[i] for i in (0,1,4,5,12,14)))
    bilinears=[code.mul((0,0,p),code.mul(gamma[i],gamma[j]))
               for p,i,j in ((1,0,1),(1,4,5),(3,12,14))]
    assert code.product(bilinears)==logical
    assert all(code.adj(b)==b for b in bilinears)
    assert all(code.commute(a,b) for a in bilinears for b in bilinears)

    def scalar_trace(a):
        kind,phase=comp(a)
        if kind!='scalar':return F(0)
        assert phase in (0,2)
        return F(1 if phase==0 else -1)
    def logical_trace(a):return scalar_trace(code.mul(logical,a))
    sixth_axis=[];sixth_logical=0;sixth_scalar=0
    for ids in itertools.combinations(range(20),6):
        word=code.mul((0,0,1),code.product(gamma[i] for i in ids))
        kind,_=comp(word);sixth_logical+=kind=='logical';sixth_scalar+=kind=='scalar'
        coefficient=logical_trace(word)
        if coefficient:sixth_axis.append(dict(indices=list(ids),sign=int(coefficient)))
    assert sixth_logical==240 and sixth_scalar==0 and len(sixth_axis)==80
    def add(*polys):
        out={}
        for p in polys:
            for a,c in p.items():out[a]=out.get(a,F(0))+c
        return {a:c for a,c in out.items() if c}
    def scale(p,c):return {a:b*c for a,b in p.items() if b*c}
    def mul(p,q):
        out={}
        for a,c in p.items():
            for b,d in q.items():
                ab=code.mul(a,b);out[ab]=out.get(ab,F(0))+c*d
        return {a:c for a,c in out.items() if c}
    def logical_coefficient(p):return sum((c*logical_trace(a) for a,c in p.items()),F(0))
    def scalar_coefficient(p):return sum((c*scalar_trace(a) for a,c in p.items()),F(0))
    b={a:F(i+1) for i,a in enumerate(bilinears)}
    powers=[{ident:F(1)}]
    for _ in range(7):powers.append(mul(powers[-1],b))
    coefficients=[logical_coefficient(p) for p in powers]
    assert coefficients[:4]==[0,0,0,36]
    assert scalar_coefficient(powers[2])==14
    assert coefficients[5]==1680
    outside=F(1,10);bias=F(3,5);delta=(1-outside)*bias
    alpha=F(2,7)
    # Formal solution of q + alpha*q^3 + h*B = 0, through h^7.
    series=[{},scale(b,-1)]
    for n in range(2,8):
        terms=[mul(mul(series[i],series[j]),series[n-i-j])
               for i in range(1,n) for j in range(1,n-i)]
        series.append(scale(add(*terms),-alpha))
    residual=[]
    for n in range(1,8):
        triple=add(*(mul(mul(series[i],series[j]),series[n-i-j])
                     for i in range(1,n) for j in range(1,n-i)))
        p=add(series[n],scale(triple,alpha),b if n==1 else {})
        # Pauli phases can encode equal matrices with opposite signs.
        reduced={}
        for (x,z,phase),c in p.items():
            key=(x,z,phase%2)
            reduced[key]=reduced.get(key,F(0))+c*(1 if phase<2 else -1)
        residual.append(all(c==0 for c in reduced.values()))
    assert all(residual)
    difference=[delta*logical_coefficient(p) for p in series]
    assert difference[:3]==[0,0,0]
    assert difference[3]==F(972,175)
    assert all(difference[n]==0 for n in (0,2,4,6))
    rows=[]
    for h in (1/16,1/32,1/64,1/128,1/256):
        value=0.
        for signs in itertools.product((-1,1),repeat=3):
            eigen=sum((i+1)*s for i,s in enumerate(signs))
            q=-h*eigen
            for _ in range(12):q-=(q+float(alpha)*q**3+h*eigen)/(1+3*float(alpha)*q*q)
            assert abs(q+float(alpha)*q**3+h*eigen)<1e-14
            value+=float(delta)*np.prod(signs)/8*q
        trunc=sum(float(c)*h**n for n,c in enumerate(difference))
        rows.append(dict(hbar=h,exact_auxiliary_response_difference=value,
                         difference_over_hbar_cubed=value/h**3,
                         formal_hbar_seven_error=abs(value-trunc)))
    assert abs(rows[-1]['difference_over_hbar_cubed']-float(difference[3]))<.005
    assert rows[-1]['formal_hbar_seven_error']<2e-15
    # Changing the code-complement weight spoils fourth-moment equality.
    parity_four=code.product(gamma[:4]);assert scalar_trace(parity_four)==1
    def fourth(e):return 1-e-F(2,1022)*e
    varying_weight_delta=fourth(F(1,5))-fourth(F(1,10))
    assert varying_weight_delta==F(-256,2555)
    return dict(round=845,all_checks_passed=True,fresh_test_groups=1,
        original_code_CAR_modes=10,exact_rational_and_symplectic_arithmetic=True,
        six_Majorana_subsets_checked=38760,sixth_degree_logical_words=sixth_logical,
        sixth_degree_scalar_words=sixth_scalar,
        chosen_logical_axis_sixth_moment_support=sixth_axis,
        diagnostic_bilinears_Majorana_pairs=[[0,1],[4,5],[12,14]],
        B_coefficients=[1,2,3],logical_coefficients_B_powers=[str(c) for c in coefficients],
        B_squared_code_scalar=str(scalar_coefficient(powers[2])),
        formal_auxiliary_response_difference=[str(c) for c in difference],
        auxiliary_equation='q+(2/7)*q^3+hbar*B=0',
        auxiliary_equation_residual_through_hbar_seven_exactly_zero=all(residual),
        independent_numerical_roots=rows,
        changing_complement_weight_creates_hbar_squared_difference=str(varying_weight_delta),
        normalized_logical_observable_has_hbar_zero_difference=str(delta),
        original_continuum_actual_six_leg_response_computed=False,
        auxiliary_scalar_feedback_is_original_action=False,
        hbar_filtration_is_real_time_or_coupling_expansion=False,
        arbitrary_all_order_content_independence_proven=False,
        finite_coupling_continuum_control_proven=False)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
