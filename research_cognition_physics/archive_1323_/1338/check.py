"""Finite checks for the identical-reference equivariant class in round 1338.
Run with --write to create the result once. The default is a read-only replay.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def dag(a):
    return a.conj().T


def tensor(args):
    out=np.array([1.])
    for a in args:
        out=np.kron(out,a)
    return out


def rotation(n,rng):
    q,r=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
    d=np.diag(r)
    return q*(d/np.abs(d))


def move_slots(n, permutation):
    # Rows are output slots in the order given; columns retain the input order.
    channels=len(permutation)
    dim=n**channels
    indices=np.arange(dim).reshape((n,)*channels).transpose(permutation).ravel()
    return np.eye(dim)[indices]


def interaction(n,k):
    slots=k+1
    h=np.zeros((n**slots,n**slots),complex)
    for i in range(k):
        perm=list(range(slots))
        perm[i],perm[-1]=perm[-1],perm[i]
        h+=(.25+.17*i)*move_slots(n,perm)
    perm=move_slots(n,list(range(1,slots))+[0])
    weight=.31+.27j
    h+=weight*perm + weight.conjugate()*dag(perm)
    eigen,v=np.linalg.eigh(h)
    return (v*np.exp(-1j*eigen))@dag(v)


class Checks:
    def __init__(self):
        self.n={}
        self.err={}

    def close(self,key,value,ans,tol=2e-9):
        e=float(np.max(np.abs(np.asarray(value)-np.asarray(ans)))/max(1.,np.max(np.abs(ans))))
        if not np.isfinite(e) or e>tol:
            raise AssertionError((key,e))
        self.n[key]=self.n.get(key,0)+1
        self.err[key]=max(e,self.err.get(key,0.))


def partial_port(w,f,n):
    return np.kron(f.conj()[None,:],np.eye(n))@w


def masses(w,chi,t,m):
    n=w.shape[1]
    # Only WW* is physical for this check. Any completion implements the
    # same nonzero spectrum. A complete QR supplies the remaining columns.
    rest=np.linalg.qr(w,mode="complete")[0][:,n:]
    mat=np.vstack((t*np.kron(chi.conj()[None,:],np.eye(n)),m*dag(rest)))
    return np.linalg.svd(mat,compute_uv=False)[::-1]


def repeated_test(n,k,u,rng,log):
    q=rotation(n,rng)
    r,x=q[:,:2].T
    refs=tensor([r]*k)
    w=u@np.kron(refs[:,None],np.eye(n))
    incoming=tensor([r]*k+[x])
    final=u@incoming
    a=np.vdot(incoming,final)
    b=[]
    for j in range(k):
        orbit=[r]*(k+1)
        orbit[j]=x
        b.append(np.vdot(tensor(orbit),final))
    all_same=tensor([r]*(k+1))
    z=np.vdot(all_same,u@all_same)
    log.close("norm_of_one_excitation",abs(a)**2+np.linalg.norm(b)**2,1)
    log.close("symmetric_sector",z,a+sum(b))
    log.close("unit_phase",abs(z),1.)
    data_x=(final.reshape(-1,n)@x.conj())
    log.close("independent_calibration",np.vdot(data_x,data_x),abs(a)**2)
    if n<3:
        return {}
    joint=rotation(n**k,rng)
    summary={"n":n,"k":k,"calibration":float(abs(a)**2),"ports":[]}
    for chi,t in zip(joint[:,:2].T,(.27,.49)):
        gamma=np.vdot(chi,refs)
        l=partial_port(w,chi,n)
        log.close("rank_one_range",(np.eye(n)-np.outer(r,r.conj()))@(l-a*gamma*np.eye(n)),0.)
        log.close("response_on_reference",l@r,z*gamma*r)
        s=np.linalg.svd(l,compute_uv=False)[::-1]
        log.close("exact_middle",s[1:-1],np.repeat(abs(a*gamma),n-2))
        if s[0]<1e-8:
            raise AssertionError("ordinary test unexpectedly rank-deficient")
        log.close("common_curvature",s[0]*s[-1]/s[1]**2,1/abs(a))
        mass=1.8
        actual=masses(w,chi,t,mass)
        low=actual[:n]
        sigma=low*np.sqrt(mass**2+t**2-low**2)/(mass*t)
        log.close("finite_reconstruction",sigma,s)
        log.close("finite_flat_partners",actual[n:-n],np.repeat(mass,len(actual)-2*n))
        ratio=low[0]*low[-1]/low[1]**2
        log.close("finite_lower_bound",max(0.,1/abs(a)-ratio),0.)
        summary["ports"].append({"leading":float(s[0]*s[-1]/s[1]**2),"finite":float(ratio)})
    # Two fixed joint bras orthogonal to the all-r tensor exemplify the
    # universal-calibration obstruction. Existence for any pair is analytic.
    raw=joint[:,:2]-np.outer(refs,refs.conj()@joint[:,:2])
    for chi in np.linalg.qr(raw)[0].T:
        s=np.linalg.svd(partial_port(w,chi,n),compute_uv=False)
        log.close("zero_overlap_rank_at_most_one",s[1],0.)
    return summary


def span_test(n,k,u,d,rng,log):
    basis=rotation(n,rng)
    # Span exactly d even when some references repeat. d<n gives a probe.
    indices=np.arange(k)%d
    refs=tensor([basis[:,i] for i in indices])
    w=u@np.kron(refs[:,None],np.eye(n))
    chi=rotation(n**k,rng)[:,0]
    l=partial_port(w,chi,n)
    x=basis[:,d]
    scalar=np.vdot(x,l@x)
    proj=np.eye(n)-basis[:,:d]@dag(basis[:,:d])
    log.close("multi_direction_range",proj@(l-scalar*np.eye(n)),0.)
    if n>2*d:
        error=np.abs(np.linalg.svd(l,compute_uv=False)-abs(scalar))
        log.close("reference_span_multiplicity",np.sort(error)[:n-2*d],0.)
    return {"n":n,"k":k,"d":d,"formed":np.linalg.svd(l,compute_uv=False).tolist()}


def calculate():
    log=Checks()
    rng=np.random.default_rng(1338)
    equal=[]
    span=[]
    for n,k in ((2,2),(3,1),(3,2),(3,3),(4,2),(6,2)):
        u=interaction(n,k)
        r=rotation(n,rng)
        common=np.array([[1.]])
        for _ in range(k+1):
            common=np.kron(common,r)
        log.close("common_unitary_covariance",u@common,common@u)
        ret=repeated_test(n,k,u,rng,log)
        if ret:
            equal.append(ret)
        span.append(span_test(n,k,u,min(2,k,n-1),rng,log))
    return {"round":1338,"status":"PASS","equal_reference":equal,
            "reference_span":span,"checks":log.n,"max_errors":log.err}


def compare(a,b):
    if isinstance(b,dict):
        if a.keys()!=b.keys():
            raise AssertionError("key mismatch")
        for k in b:
            compare(a[k],b[k])
    elif isinstance(b,list):
        if len(a)!=len(b):
            raise AssertionError("list length mismatch")
        for x,y in zip(a,b):
            compare(x,y)
    elif isinstance(b,float):
        np.testing.assert_allclose(a,b,atol=2e-9,rtol=2e-8)
    elif a!=b:
        raise AssertionError((a,b))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    output=Path(__file__).with_name("results.json")
    ans=calculate()
    if args.write:
        output.write_text(json.dumps(ans,indent=2)+"\n",encoding="utf-8")
    else:
        compare(ans,json.loads(output.read_text(encoding="utf-8")))
    print(json.dumps(dict(status=ans["status"],count=sum(ans["checks"].values()),
                          sample=ans["equal_reference"][1])))
