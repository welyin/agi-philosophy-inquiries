"""Exact Clifford reduction for the original offdiagonal centre source.
Variables: r,s,x,xp,y,yp,z,w; six vectors Ec,Fc,Gc,Ew,Fw,Gw.
"""
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
HERE=Path(__file__).resolve().parent
SHIFT=5;VAR=[1<<(SHIFT*j) for j in range(8)]


def addc(a,b):return a[0]+b[0],a[1]+b[1]
def mulc(a,b):return a[0]*b[0]+2*a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def scale(a,n):return a[0]*n,a[1]*n


def padd(out,key,value):
    if not any(value):return
    value=addc(out.get(key,(0,0)),value)
    if any(value):out[key]=value
    elif key in out:del out[key]


def pmul(a,b):
    out={}
    for u,x in a.items():
        for v,y in b.items():padd(out,u+v,mulc(x,y))
    return out


def s_terms():
    zero=(F(0),F(0));half=(F(1,2),F(0));q=(F(0),F(1,4))
    us=[[[q,addc(scale(half,-1),scale(q,-1))],[addc(half,q),q]],
        [[scale(q,-1),addc(scale(half,-1),q)],[addc(half,scale(q,-1)),scale(q,-1)]]]
    terms={}
    for t in range(2):
        for s in range(2):
            for k in range(2):
                for weak in range(2):
                    j=k+3*weak
                    for a in range(2):
                        for b in range(2):
                            value=scale(mulc(us[a][t][k],us[b][s][k]),F((-1)**(a*weak),2))
                            key=(t,s,j,a^b);terms[key]=addc(terms.get(key,zero),value)
            if t==s:
                for j in (0 if t==0 else 2,3 if t==0 else 5):
                    key=(t,s,j,0);terms[key]=addc(terms.get(key,zero),half)
    terms={k:v for k,v in terms.items() if any(v)}
    den=math.lcm(*(x.denominator for v in terms.values() for x in v))
    return den,[(t,s,j,r,(int(v[0]*den),int(v[1]*den))) for (t,s,j,r),v in terms.items()]


def gram(i,j):
    if (i<3)!=(j<3):return []
    if i>j:i,j=j,i
    table={(0,0):[(VAR[0],1)],(1,1):[(VAR[1],1)],(2,2):[(VAR[1],1)],
        (0,1):[(VAR[2],1)],(0,2):[(VAR[3],1)],(1,2):[(VAR[6],1)],
        (3,3):[(0,1),(VAR[0],-1)],(4,4):[(0,1),(VAR[1],-1)],(5,5):[(0,1),(VAR[1],-1)],
        (3,4):[(VAR[4],1)],(3,5):[(VAR[5],1)],(4,5):[(VAR[7],1)]}
    return table[i,j]


def right_vector(mask,j):
    out=[];indices=[i for i in range(6) if mask>>i&1]
    if not mask>>j&1:out.append((mask^(1<<j),0,(-1)**sum(i>j for i in indices)))
    for k,i in enumerate(indices):
        for mon,coef in gram(i,j):out.append((mask^(1<<i),mon,coef*(-1)**(len(indices)-1-k)))
    return out


def advance(power,terms):
    out={}
    for (t,k,mask,r),poly in power.items():
        for a,s,j,b,coef in terms:
            if k!=a:continue
            coef=scale(coef,(-1)**(r*(j>=3)))
            for mask2,mon,sgn in right_vector(mask,j):
                dest=out.setdefault((t,s,mask2,r^b),{})
                for old,value in poly.items():padd(dest,old+mon,scale(mulc(value,coef),sgn))
    return {k:p for k,p in out.items() if p}


def trace(power):
    out={}
    for t in range(2):
        for mon,coef in power.get((t,t,0,0),{}).items():padd(out,mon,scale(coef,2))
    return out


def run():
    den,terms=s_terms();powers=[{(t,t,0,0):{0:(1,0)} for t in range(2)}]
    counts=[]
    for n in range(1,9):
        powers.append(advance(powers[-1],terms))
        counts.append(sum(len(p) for p in powers[-1].values()))
        print('power',n,'terms',counts[-1],flush=True)
    traces=[None]+[trace(powers[2*k]) for k in range(1,5)]
    elementary=[{0:(1,0)}]
    for k in range(1,5):
        poly={}
        for j in range(1,k+1):
            for mon,c in pmul(elementary[k-j],traces[j]).items():padd(poly,mon,scale(c,(-1)**(j-1)))
        assert all(c[0]%k==0 and c[1]%k==0 for c in poly.values())
        elementary.append({m:(c[0]//k,c[1]//k) for m,c in poly.items()})
    remainder={}
    for k in range(5):
        for key,p in powers[8-2*k].items():
            dest=remainder.setdefault(key,{})
            for mon,c in pmul(elementary[k],p).items():padd(dest,mon,scale(c,(-1)**k))
    remainder={k:p for k,p in remainder.items() if p}
    print('quartic identity remainder',len(remainder),'e4 terms',len(elementary[4]),flush=True)
    result=dict(denominator=den,linear_terms=terms,power_monomial_counts=counts,
        characteristic_identity_zero=not remainder,
        elementary=[[[mon,list(c)] for mon,c in sorted(p.items())] for p in elementary],
        e4_denominator=den**8,variables=['r','s','x','xp','y','yp','z','w'])
    with (HERE/'clifford_moment_probe_results.json').open('x',encoding='utf8') as f:json.dump(result,f,indent=2)
    return result


if __name__=='__main__':run()
