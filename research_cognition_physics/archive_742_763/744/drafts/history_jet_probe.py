"""Exploratory full onsite H jets on normal Gaussian H5 packets.

Sparse polynomials multiply exp(-alpha |x|^2/2). All five target derivatives
are kept before evaluation on a Higgs radial representative. This is a native
Hamiltonian Taylor check, not a frozen classical-background time simulation.
"""
from functools import lru_cache
from math import comb
from pathlib import Path
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parent.parent;sys.path.insert(0,str(ROOT))
import joint_record_mass_feedback as old
ZERO=(0,)*5
ALPHA=1.
W=1.

def add(out,term,factor=1):
    for key,v in term.items():
        out[key]=out.get(key,0)+factor*v
    return {k:v for k,v in out.items() if abs(v)>1e-12}

def shift(poly,axis,count=1,factor=1):
    out={}
    for (state,power),v in poly.items():
        pp=list(power);pp[axis]+=count;out[state,tuple(pp)]=factor*v
    return out

def radial(poly):
    out={}
    for i in range(5):out=add(out,shift(poly,i,2))
    return out

def euler(poly):return {key:sum(key[1])*v for key,v in poly.items()}

def laplace(poly):
    out={}
    for (state,power),v in poly.items():
        for i,d in enumerate(power):
            if d>=2:
                pp=list(power);pp[i]-=2;key=(state,tuple(pp));out[key]=out.get(key,0)+d*(d-1)*v
    return out

def D(poly):return add(euler(poly),radial(poly),-ALPHA)

def kinetic(poly):
    out=add(laplace(poly),euler(poly),-2*ALPHA)
    out=add(out,radial(poly),ALPHA*ALPHA);out=add(out,poly,-5*ALPHA)
    out=add(out,D(D(poly)),1/6);out=add(out,D(poly),4/6)
    return {k:-v/(2*W) for k,v in out.items()}

def variable(poly,indices):
    out={}
    for i in indices:out=add(out,shift(poly,i,2))
    return out

def potential(poly):
    # Exact623 polynomial in a=|xH|^2,b=x5^2, independently from old parameters.
    L,u,_=old.geom.lattice.scalar.parameters();P=np.eye(2)-np.outer(u,np.ones(2))/(6*old.geom.M)
    M=P.T@L@P/4;linear=-(P.T@L@u)/(2*old.geom.M);constant=u@L@u/(4*old.geom.M**2)
    a=variable(poly,range(4));b=variable(poly,[4])
    out={k:constant*v for k,v in poly.items()}
    for q,c in ((a,linear[0]),(b,linear[1]),(variable(a,range(4)),M[0,0]),
                (variable(a,[4]),2*M[0,1]),(variable(b,[4]),M[1,1])):out=add(out,q,c)
    return {k:W*v for k,v in out.items()}

@lru_cache(maxsize=None)
def car_step(a,state):
    h,d=old.LINEAR_MASS[a]
    return old.car.quadratic({state:1.},h,d)

def mass(poly):
    out={}
    for (state,power),v in poly.items():
        for a in range(5):
            for target,z in car_step(a,state).items():
                pp=list(power);pp[a]+=1;key=(target,tuple(pp));out[key]=out.get(key,0)+v*z
    return out

def H(poly):return add(add(kinetic(poly),potential(poly)),mass(poly))

def jets(state):
    result=[{(state,ZERO):1.}]
    for _ in range(4):result.append(H(result[-1]))
    return result

def evaluate(poly,r,s):
    out={}
    for (state,p),v in poly.items():
        if p[1] or p[2] or p[3]:continue
        if state not in out:out[state]=np.zeros_like(r,dtype=complex)
        out[state]+=v*r**p[0]*s**p[4]
    return out

def inner(a,b):
    return sum(z.conj()*b[state] for state,z in a.items() if state in b)

def measure(jet,order):
    u,wu=np.polynomial.laguerre.laggauss(order);v,wv=np.polynomial.hermite.hermgauss(order)
    r,s=np.meshgrid(np.sqrt(u/ALPHA),v/np.sqrt(ALPHA),indexing='ij')
    weight=wu[:,None]*wv[None,:]*u[:,None]/np.sqrt(1+(r*r+s*s)/6);weight/=weight.sum()
    F=2/(1+(r*r+s*s)/6);A=np.sin(np.sqrt(F)*s)**2/4
    values=[evaluate(p,r,s) for p in jet]
    rows={}
    for degree in (0,1,2,3,4):
        integrand=(1j**degree)*sum((-1)**k*comb(degree,k)*inner(values[degree-k],values[k]) for k in range(degree+1))
        val=np.sum(weight*A*integrand)
        rows[str(degree)]=[float(val.real),float(val.imag)]
    identity4=sum((-1)**k*comb(4,k)*inner(values[4-k],values[k]) for k in range(5))
    return dict(derivatives=rows,identity_fourth=float(np.sum(weight*identity4).real))

if __name__=='__main__':
    first=jets(0);second=jets((1<<30)|(1<<31))
    results=dict(alpha=ALPHA,w=W,terms_empty=[len(x) for x in first],terms_pair=[len(x) for x in second],rows=[])
    for order in (32,48,64):
        a=measure(first,order);b=measure(second,order)
        results['rows'].append(dict(order=order,empty=a,pair=b,difference_fourth=b['derivatives']['4'][0]-a['derivatives']['4'][0]))
    target=Path(__file__).with_name('history_jet_probe_results.json')
    with target.open('x',encoding='utf8') as f:json.dump(results,f,ensure_ascii=False,indent=2)
    print(json.dumps(results,indent=2))
