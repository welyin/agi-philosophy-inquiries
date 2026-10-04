"""Exact two-sphere moments of the verified Clifford determinant polynomial."""
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import clifford_moment_probe as alg
import joint_full_holonomy_reduction as original
TARGET=HERE/'offdiagonal_exact_sphere_results.json'


def rising(a,n):return math.prod(range(a,a+n))


@lru_cache(None)
def wick(a,b):
    out=[]
    for k in range(min(a,b)+1):
        if (a-k)%2 or (b-k)%2:continue
        i=(a-k)//2;j=(b-k)//2
        coef=math.factorial(a)*math.factorial(b)//(math.factorial(k)*math.factorial(i)*math.factorial(j)*2**(i+j))
        out.append((i+j,k,coef))
    return out


def integrate(data):
    e4={m:tuple(v) for m,v in data['elementary'][4]}
    squared=alg.pmul(e4,e4);fourth=alg.pmul(squared,squared)
    radial=2**16*rising(5,16);h={};kept=0
    for mon,value in fourth.items():
        a,s,xe,xg,ye,yg,z,w=[(mon>>(5*j))&31 for j in range(8)]
        pc,pw=xe+xg,ye+yg
        if pc%2 or pw%2:continue
        n=a+(pc+pw)//2;assert n<=16
        pref=radial//(2**n*rising(5,n))*2**a*rising(3+pc//2,a)
        for sc,zc,cc in wick(xe,xg):
            for sw,wc,cw in wick(ye,yg):
                for j in range(sw+1):
                    key=(s+sc+j,z+zc,w+wc)
                    assert sum(key)<=16
                    factor=pref*cc*cw*math.comb(sw,j)*(-1)**j
                    alg.padd(h,key,alg.scale(value,factor))
        kept+=1
    hden=data['e4_denominator']**4*radial
    coeff={};fden=rising(5,16)
    for (a,b,c),value in h.items():
        factor=fden//rising(5,a+b+c)*rising(3+b,a)*math.factorial(b)*math.factorial(c)
        alg.padd(coeff,(b,c),alg.scale(value,factor))
    den=hden*fden
    gcd=math.gcd(den,*(abs(x) for v in coeff.values() for x in v))
    coeff={key:(v[0]//gcd,v[1]//gcd) for key,v in coeff.items()};den//=gcd
    return dict(e4_terms=len(e4),e4_squared_terms=len(squared),fourth_terms=len(fourth),
        nonzero_even_E_monomials=kept,after_first_sphere_terms=len(h),
        auxiliary_character_coefficients=[[list(key),list(v)] for key,v in sorted(coeff.items())],
        denominator=den,max_total_character_degree=max(sum(k) for k in coeff),
        coefficient_field='Q(sqrt(2))',full_two_S9_average=True)


def peval(poly,variables,den):
    return sum((a+b*np.sqrt(2))*math.prod(variables[j]**((m>>(5*j))&31) for j in range(8)) for m,(a,b) in poly.items())/den


def matrix_check(data):
    rng=np.random.default_rng(69921);e4={m:tuple(v) for m,v in data['elementary'][4]};rows=[]
    for k in range(6):
        angles=rng.uniform(-1,1,(1,4));h,_=original.probe.torus(angles)
        ef=rng.normal(size=(2,10));ef[1]+=.8*ef[0];ef/=np.linalg.norm(ef,axis=1)[:,None]
        e,f=ef;g=original.b.prior.vector_rotation(np.diag(h[0])).T@f
        r=e[:6]@e[:6];s=f[:6]@f[:6]
        variables=[r,s,e[:6]@f[:6],e[:6]@g[:6],e[6:]@f[6:],e[6:]@g[6:],f[:6]@g[:6],f[6:]@g[6:]]
        d=peval(e4,variables,data['e4_denominator'])
        _,v=original.polar(h[0],True,False)
        t=np.zeros((32,32),complex)
        for site in range(2):t[16*site:16*site+16,16*site:16*site+16]=sum(x*y for x,y in zip(ef[site],original.b.internal.T))
        aux=(t+v.conj()@t@v.conj().T)/2
        ratio=float(abs(np.linalg.det(aux)-d**4)/max(abs(d**4),1e-30));assert ratio<3e-10,(k,ratio,d,np.linalg.det(aux))
        whole,gap=original.probe.entry.weight(True,False,np.eye(16),np.diag(h[0]),ef)
        calculated=np.linalg.det((np.eye(32)+v)/2)**2*d**4
        err=float(abs(whole-calculated)/max(abs(whole),1e-30));assert err<3e-10
        rows.append(dict(sample=k,determinant_polynomial=d,compressed_relative_error=ratio,
            original256_relative_error=err,gap=gap))
    return rows


def run():
    data=json.loads((HERE/'clifford_moment_probe_results.json').read_text('utf8'))
    assert data['characteristic_identity_zero']
    checks=matrix_check(data);result=integrate(data);result['original_matrix_checks']=checks
    return result


if __name__=='__main__':
    result=run()
    with TARGET.open('x',encoding='utf8') as f:json.dump(result,f,indent=2)
    print(json.dumps({k:result[k] for k in ('fourth_terms','after_first_sphere_terms','max_total_character_degree')}),flush=True)
