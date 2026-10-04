"""Test whether lower-order native population response is preparation dependent.

Original Hamiltonian, explicitly different even anisotropic Gaussian readies.
One vertex only; this is diagnostic, not a connected-graph certificate.
"""
from pathlib import Path
import sys,json
from math import comb
import numpy as np
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'round744_drafts'))
import history_jet_probe as old
native_kinetic=old.kinetic

def weighted_radius(poly,alpha):
    out={}
    for i,a in enumerate(alpha):out=old.add(out,old.shift(poly,i,2),a)
    return out

def kinetic(poly,alpha):
    e={key:sum(a*n for a,n in zip(alpha,key[1]))*v for key,v in poly.items()}
    D=lambda p:old.add(old.euler(p),weighted_radius(p,alpha),-1)
    out=old.add(old.laplace(poly),e,-2)
    out=old.add(out,weighted_radius(poly,[a*a for a in alpha]))
    out=old.add(out,poly,-sum(alpha))
    out=old.add(out,D(D(poly)),1/6);out=old.add(out,D(poly),4/6)
    return {k:-v/2 for k,v in out.items()}

def measure(jet,beta,n=64):
    u,wu=np.polynomial.laguerre.laggauss(n);v,wv=np.polynomial.hermite.hermgauss(n)
    r,s=np.meshgrid(np.sqrt(u),v/np.sqrt(beta),indexing='ij')
    weight=wu[:,None]*wv[None,:]*u[:,None]/np.sqrt(1+(r*r+s*s)/6);weight/=weight.sum()
    A=np.sin(np.sqrt(2/(1+(r*r+s*s)/6))*s)**2/4
    values=[old.evaluate(p,r,s) for p in jet]
    rows=[]
    for degree in range(5):
        f=1j**degree*sum((-1)**k*comb(degree,k)*old.inner(values[degree-k],values[k]) for k in range(degree+1))
        rows.append(float(np.sum(weight*A*f).real))
    return rows

rows=[]
for beta in (.5,1.,2.):
    alpha=[1.]*4+[beta]
    old.kinetic=lambda p:kinetic(p,alpha)
    a,b=old.jets(0),old.jets((1<<30)|(1<<31))
    da,db=measure(a,beta),measure(b,beta)
    rows.append(dict(beta=beta,empty=da,pair=db,differences=[b-a for a,b in zip(da,db)]))
old.kinetic=native_kinetic
result=dict(scope='Same native one-vertex H, new normal symmetric ready parameters; not graph proof.',rows=rows)
with Path(__file__).with_name('local_population_jet_probe_results.json').open('x',encoding='utf8') as f:
    json.dump(result,f,ensure_ascii=False,indent=2)
print(json.dumps(rows))
