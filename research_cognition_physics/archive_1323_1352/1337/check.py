"""Round 1337. Compute the formed code and tree spectrum directly from SWAP.
Default compares against results.json. Use --write only to generate a result.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def dag(a):
    return a.conj().T


def swap(n):
    eye = np.eye(n*n)
    return eye.reshape(n,n,n,n).transpose(1,0,2,3).reshape(n*n,n*n)


def pulse(n, theta):
    return np.cos(theta)*np.eye(n*n) - 1j*np.sin(theta)*swap(n)


def unitary(n,rng):
    q,r = np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
    d = np.diag(r)
    return q*(d/np.abs(d))


def spectral_prediction(n,theta,overlap):
    c=abs(np.cos(theta))
    s=abs(np.sin(theta))
    middle=c*overlap
    product=c*overlap**2
    trace=s*s+2*middle*middle
    high2=(trace+np.sqrt(max(0.,trace*trace-4*product*product)))/2
    if high2==0:
        return np.zeros(n)
    high=np.sqrt(high2)
    return np.r_[product/high, np.repeat(middle,n-2), high]


def root_spectrum(s, n, t, mass):
    total=mass**2+t**2
    prod=(mass*t*s)**2
    # (mass^2-t^2)^2+4 prod*(1/s^2-1) is well-defined even at equal mass.
    disc=(mass**2-t**2)**2+4*(mass*t)**2*np.maximum(1-s*s,0.)
    high2=(total+np.sqrt(disc))/2
    return np.sort(np.r_[np.sqrt(prod/high2),np.sqrt(high2),
                         np.repeat(mass,n*(n-2))])


class Audit:
    def __init__(self):
        self.counts={}
        self.worst={}

    def close(self,name,a,b,tol=3e-8):
        error=float(np.max(np.abs(np.asarray(a)-np.asarray(b)))/max(1.,np.max(np.abs(b))))
        if not np.isfinite(error) or error>tol:
            raise AssertionError((name,error,tol))
        self.counts[name]=1+self.counts.get(name,0)
        self.worst[name]=max(error,self.worst.get(name,0.))


def physical_case(log,n,basis,chi,theta,t,mass):
    # Tensor order is (reference, data), r is initially the first basis column.
    r=basis[:,0]
    formed=pulse(n,theta)@np.kron(basis,np.eye(n))
    w,perp=formed[:,:n],formed[:,n:]
    bra=np.kron(chi.conj()[None,:],np.eye(n))
    l=bra@w
    a=np.vdot(chi,r)
    calc=np.cos(theta)*a*np.eye(n)-1j*np.sin(theta)*np.outer(r,chi.conj())
    log.close("isometry",dag(w)@w,np.eye(n))
    log.close("formed_single_port",l,calc)
    s=np.linalg.svd(l,compute_uv=False)[::-1]
    expected=spectral_prediction(n,theta,abs(a))
    log.close("all_formed_singular_values",s,expected)
    full=np.vstack((t*bra,mass*dag(perp)))
    measured=np.linalg.svd(full,compute_uv=False)[::-1]
    log.close("complete_tree_spectrum",measured,root_spectrum(expected,n,t,mass))
    light=measured[:n]
    inferred=light*np.sqrt(np.maximum(t*t+mass*mass-light*light,0.))/(t*mass)
    log.close("reconstructed_formed_spectrum",inferred,s)
    c=abs(np.cos(theta))
    out=dict(formed=s.tolist(), tree_light=light.tolist())
    if n>=3 and abs(a)*c>1e-5:
        log.close("dark_multiplicity_formed",s[1:-1],np.repeat(abs(a)*c,n-2))
        log.close("formed_consecutive_ratio",s[-1]*s[0]/s[1]**2,1/c)
        log.close("deembedded_ratio",inferred[-1]*inferred[0]/inferred[1]**2,1/c)
        ratio=light[0]*light[-1]/light[1]**2
        log.close("convex_lower_bound_excess",max(0.,1/c-ratio),0.,3e-8)
        out.update(leading_ratio=float(s[-1]*s[0]/s[1]**2),finite_ratio=float(ratio))
    return out


def calibration(n,theta):
    # x orthogonal to r; ancillary reference is a second independently prepared r.
    r=np.eye(n)[:,0]
    x=np.eye(n)[:,1]
    state=(pulse(n,theta)@np.kron(r,x)).reshape(n,n)
    final=np.einsum("ij,k->ijk",state,r)
    antisym=(final-final.transpose(2,1,0))/2
    return float(np.linalg.norm(antisym)**2)


def check_partner(log):
    n=3
    rng=np.random.default_rng(20261118)
    basis=unitary(n,rng)
    r=basis[:,0]
    partner_dim=2
    ent=rng.normal(size=n*partner_dim)+1j*rng.normal(size=n*partner_dim)
    ent/=np.linalg.norm(ent)
    w=pulse(n,.8)@np.kron(r[:,None],np.eye(n))
    formed=np.kron(w,np.eye(partner_dim))@ent
    record=[]
    for chi in basis.T:
        l=np.kron(chi.conj()[None,:],np.eye(n))@w
        branch=np.kron(l,np.eye(partner_dim))@ent
        record.append(np.vdot(branch,branch).real)
    log.close("complete_record",sum(record),1.)
    before=ent.reshape(n,partner_dim)
    after=formed.reshape(n*n,partner_dim)
    log.close("partner_marginal",dag(before)@before,dag(after)@after)
    return record


def calculate():
    rng=np.random.default_rng(1337)
    log=Audit()
    cases=0
    for n in (2,3,4,6):
        basis=unitary(n,rng)
        for theta in (0.,.47,np.arccos(.6),np.pi/2,2.2):
            for x in (0.,.19,.73,1.):
                chi=x*basis[:,0]+np.sqrt(1-x*x)*np.exp(.73j)*basis[:,1]
                for t,mass in ((.3,1.7),(1.05,.8)):
                    physical_case(log,n,basis,chi,theta,t,mass)
                    cases+=1
            log.close("independent_contact_calibration",calibration(n,theta),
                      np.sin(theta)**2/2)
    n,theta=3,math.acos(3/5)
    basis=np.eye(n,dtype=complex)
    portal_u=np.array([3,4,0],complex)/5
    portal_d=np.array([5,0,12j],complex)/13
    fixed=dict(calibration_probability=calibration(n,theta))
    for name,chi,t in (("u",portal_u,.3),("d",portal_d,.45)):
        fixed[name]=physical_case(log,n,basis,chi,theta,t,1.7)
    benchmark=dict(up=.00127*171.7/.619**2, down=.00290*2.89/.055**2,
                   up_rectangular_upper=.00177*174.7/.535**2)
    if not(benchmark["up"]<1<benchmark["up_rectangular_upper"]):
        raise AssertionError("historical central vs rectangle")
    return dict(round=1337,status="PASS",cases=cases,finite_sample=fixed,
                historical_common_scale=benchmark,
                partner_example_branch_probabilities=check_partner(log),
                checks=log.counts,worst_normalized_errors=log.worst)


def compare(a,b):
    if isinstance(b,dict):
        if a.keys()!=b.keys():
            raise AssertionError("result keys differ")
        for key in a:
            compare(a[key],b[key])
    elif isinstance(b,(float,list)):
        np.testing.assert_allclose(a,b,rtol=5e-8,atol=3e-8)
    elif a!=b:
        raise AssertionError((a,b))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    path=Path(__file__).with_name("results.json")
    answer=calculate()
    if args.write:
        path.write_text(json.dumps(answer,indent=2)+"\n",encoding="utf-8")
    else:
        compare(answer,json.loads(path.read_text(encoding="utf-8")))
    print(json.dumps(dict(status=answer["status"],cases=answer["cases"],
                          checks=sum(answer["checks"].values()),
                          sample=answer["finite_sample"])))
