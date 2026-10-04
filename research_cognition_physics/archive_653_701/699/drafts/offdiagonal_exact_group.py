"""Exact original Haar integral from the complete two-sphere polynomial."""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import clifford_moment_probe as alg
import joint_static_haar_majorant as previous
TARGET=HERE/'offdiagonal_exact_group_results.json'


def prime_stream():
    p=2000000000-(2000000000-1)%120
    while True:
        if previous.isprime(p):yield p
        p-=120


def primitive_root(n,p):
    factors=[q for q in range(2,n+1) if n%q==0 and previous.isprime(q)]
    for a in range(2,1000):
        r=pow(a,(p-1)//n,p)
        if pow(r,n,p)==1 and all(pow(r,n//q,p)!=1 for q in factors):return r
    raise AssertionError('root')


def setup(data):
    kappa=(1,0)
    for _ in range(16):kappa=alg.mulc(kappa,(2,1))
    coeff={tuple(key):alg.mulc(tuple(value),kappa) for key,value in data['auxiliary_character_coefficients']}
    common=data['denominator']*2**16*16**16*2**8*12
    return coeff,common


def modular(p,n,coeff,auxden):
    coords=previous.chamber(n);r=primitive_root(n,p);roots=np.array([pow(r,k,p) for k in range(n)],dtype=np.int64)
    count=len(coords);indices=[coords[:,j] for j in range(4)]+[(-coords.sum(axis=1))%n]
    cos=[(roots[j]+roots[(-j)%n])*pow(2,-1,p)%p for j in indices]
    hs=[]
    for start,end in ((0,3),(3,5)):
        h=[np.ones(count,dtype=np.int64)]+[np.zeros(count,dtype=np.int64) for _ in range(8)]
        for c in cos[start:end]:
            for k in range(1,9):h[k]=(h[k]+c*h[k-1])%p
        hs.append(h)
    values=[np.zeros(count,dtype=np.int64),np.zeros(count,dtype=np.int64)]
    for (b,c),pair in coeff.items():
        basis=hs[0][b]*hs[1][c]%p
        for j in range(2):values[j]=(values[j]+basis*(pair[j]%p))%p
    physical=np.ones(count,dtype=np.int64)
    for w,_ in previous.moment.weights():
        exponent=(coords@np.array(w,dtype=np.int64))%n
        physical=physical*((2+roots[exponent]+roots[(-exponent)%n])%p)%p
    density=np.ones(count,dtype=np.int64)
    for i,j in ((0,1),(0,2),(1,2),(3,4)):
        d=(indices[i]-indices[j])%n
        density=density*((2-roots[d]-roots[(-d)%n])%p)%p
    factor=pow((auxden*2**16*16**16*n**4)%p,-1,p)
    answer=[int(np.sum(v*physical%p*density%p,dtype=np.int64))%p*factor%p for v in values]
    assert 2*p*p<2**63 and count*p<2**63
    return answer,r,count


def run():
    data=json.loads((HERE/'offdiagonal_exact_sphere_results.json').read_text('utf8'))
    coeff,den=setup(data);primegen=prime_stream();modulus=1;integers=[0,0];rows=[]
    while modulus<=2*den:
        p=next(primegen);values,r,count=modular(p,20,coeff,data['denominator'])
        residues=[v*(den%p)%p for v in values]
        for j in range(2):integers[j]+=modulus*((residues[j]-integers[j])%p*pow(modulus,-1,p)%p)
        modulus*=p;rows.append(dict(prime=p,root=r,residues=residues))
    signed=[v if v<=modulus//2 else v-modulus for v in integers]
    assert all(abs(v)<=den for v in signed)
    result=[F(v,den) for v in signed]
    p=next(primegen);values,r,check_count=modular(p,24,coeff,data['denominator'])
    assert values==[previous.rational_mod(x,p) for x in result]
    lo,hi=F(1414213562373095,10**15),F(1414213562373096,10**15)
    assert lo*lo<2<hi*hi
    a,b=result;lower=a+b*(lo if b>=0 else hi);upper=a+b*(hi if b>=0 else lo)
    old=json.loads(previous.TARGET.read_text('utf8'))['exact_integral_upper']
    target=F(old['unproved_offdiagonal_sufficient_lower']);b00=F(old['inherited_B00']);u=F(old['exact_haar_majorant'])
    alpha=F(3,10000);negative=b00-2*alpha*lower+alpha*alpha*u
    return dict(exact_B10_Qsqrt2=[str(x) for x in result],exact_B10_bounds=[str(lower),str(upper)],
        B10_decimal=float(a)+float(b)*np.sqrt(2),certified_lower_decimal=float(lower),
        target_lower=str(target),target_met=lower>=target,
        inherited_B00=str(b00),inherited_B11_upper=str(u),test_vector=[1,str(-alpha)],
        quadratic_upper=str(negative),quadratic_upper_decimal=float(negative),
        strict_full_group_and_sphere_B_negative=negative<0,
        denominator=str(den),integer_pair=[str(x) for x in signed],CRT_modulus=str(modulus),
        primes=rows,source_degree_bounds=[19,19,19,18],grid_order=20,chamber_points=count,
        check_order=24,check_chamber_points=check_count,check_prime=p,check_root=r,
        auxiliary_polynomial_terms=len(data['auxiliary_character_coefficients']),
        not_a_finite_tau_Hb_or_original_HF_statement=True)


if __name__=='__main__':
    result=run()
    with TARGET.open('x',encoding='utf8') as f:json.dump(result,f,indent=2)
    print(json.dumps({k:result[k] for k in ('B10_decimal','target_met','quadratic_upper_decimal','strict_full_group_and_sphere_B_negative')}))
