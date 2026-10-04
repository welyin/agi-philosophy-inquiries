"""Explore an analytic complete-sphere upper majorant of the original source.
No sampled moment is used as a rigorous sign or integration certificate.
"""
import sys
from pathlib import Path
import math
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import joint_full_holonomy_reduction as old


def coefficients(v):
    a=np.zeros((20,32,32),complex)
    for site in range(2):
        for k,t in enumerate(old.b.internal.T):
            z=np.zeros((32,32),complex);z[16*site:16*site+16,16*site:16*site+16]=t
            a[10*site+k]=(z+v.conj()@z@v.conj().T)/2
    q=np.einsum('aij,bij->ab',a.conj(),a).real/32
    return a,q


def moment(q,n=16):
    # Coefficients of det(I-2 diag(t I,s I)Q)^(-1/2), homogeneous order n.
    first=q.copy();first[10:]=0
    second=q.copy();second[:10]=0
    powers=[np.eye(20)];logs=[np.zeros(1)]
    for k in range(1,n+1):
        following=[]
        for a in range(k+1):
            x=np.zeros((20,20))
            if a<k:x+=second@powers[a]
            if a>0:x+=first@powers[a-1]
            following.append(x)
        powers=following
        logs.append(np.array([np.trace(x) for x in powers])*2**(k-1)/k)
    exp=[np.ones(1)]
    for k in range(1,n+1):
        row=np.zeros(k+1)
        for j in range(1,k+1):row+=j*np.convolve(logs[j],exp[k-j])
        exp.append(row/k)
    poch=[1]
    for k in range(n):poch.append(poch[-1]*(5+k))
    return float(math.factorial(n)/2**n*sum(exp[n][a]/(poch[a]*poch[n-a]) for a in range(n+1)))


def diagnose():
    rng=np.random.default_rng(69821)
    angles=rng.uniform(0,2*np.pi,(6,4));angles[0]=0
    hs,haar=old.probe.torus(angles)
    e=rng.normal(size=(16384,2,10));e/=np.linalg.norm(e,axis=2)[:,:,None]
    rows=[]
    for k,h in enumerate(hs):
        _,v=old.polar(h,True,True);coeff,q=coefficients(v)
        analytic=moment(q)
        qs=np.einsum('na,ab,nb->n',e.reshape(-1,20),q,e.reshape(-1,20))
        weights=old.probe.values(np.repeat(h[None,:],len(e),axis=0),e,True,True)
        physical=abs(np.linalg.det((np.eye(32)+v)/2)**2)
        actual=weights.real.mean()
        rows.append(dict(sample=k,angles=angles[k].tolist(),haar=float(haar[k]),
            positive_quadratic_min_eigenvalue=float(np.linalg.eigvalsh(q).min()),
            sphere_majorant=analytic,quadratic_empirical_mean=float(np.mean(qs**16)),
            scalar_upper=physical*analytic,scalar_empirical_mean=float(actual),
            ratio=float(physical*analytic/max(actual,1e-40))))
    return rows

if __name__=='__main__':
    import json
    print(json.dumps(diagnose(),indent=2),flush=True)
