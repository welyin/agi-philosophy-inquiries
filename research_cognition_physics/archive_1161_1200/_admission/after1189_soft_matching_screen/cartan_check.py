"""Exact coefficient checks only; no numerical evolution or new physical claim."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json

def calculate():
    fixed=[]
    for c in (F(1,2),F(1),F(3,2)):
        a,b,v=F(0),F(1),F(1)
        kr=b*v-2*a/v
        kv=c*v-b/v+a/v**3
        # Expand a*v*(r-1/v²)²+b*v*(r-1/v²)+c*v.
        for r in (F(-2,3),F(0),F(4,7)):
            assert a*v*(r-1/v**2)**2+b*v*(r-1/v**2)+c*v == a*v*r*r+kr*r+kv
        fixed.append(dict(c=str(c),K_R=str(kr),K_V=str(kv),Lambda_algebraic=str(-6*kv/kr)))
    # Additional nonzero a and v values check all three coefficients.
    for a,b,c,v in ((F(2,3),F(5,4),F(-7,6),F(3,2)),(F(-1),F(2),F(3),F(2))):
        kr,kv=b*v-2*a/v,c*v-b/v+a/v**3
        for r in (F(-2,3),F(0),F(4,7)):
            assert a*v*(r-1/v**2)**2+b*v*(r-1/v**2)+c*v == a*v*r*r+kr*r+kv
    stable=[]
    for z in (F(1),F(2)):
        b,d,v=F(1),F(3),F(1)
        kappa,K=2*b*v,3*v*z*z/(4*b)
        U=24*v*v*(b/v**3+d*v)
        first=24*v*v*(-3*b/v**4+d)
        second=288*b/v**3
        mass2=second/(2*K)
        assert first==0 and second>0 and K>0
        assert mass2==192*b*b/(v**4*z*z)
        stable.append(dict(z=str(z),v=str(v),kappa=str(kappa),Lambda=str(U/(2*kappa)),K=str(K),U_prime=str(first),U_second=str(second),scalar_mass_squared=str(mass2)))
    return dict(status='PASS_ALGEBRA_ONLY',fixed_norm_coefficient_examples=fixed,known_stable_potential_examples=stable,coefficient_inputs_selected=False,actual_matter_clock_certified=False,quantum_completion_certified=False,new_numbered_round=False)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=calculate()
    out=Path(__file__).with_name('cartan_results.json')
    if args.write:
        with out.open('x',encoding='utf-8') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write('\n')
    else:
        assert json.loads(out.read_text(encoding='utf-8'))==result
    print(json.dumps(dict(status=result['status'],stable_examples=2,new_numbered_round=False)))

if __name__=='__main__':
    main()
