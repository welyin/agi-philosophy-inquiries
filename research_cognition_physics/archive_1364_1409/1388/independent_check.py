"""Orbit-intersection construction and exact cross-port endpoints, independent of author code."""
import argparse
import json
from math import cos, pi, sin, sqrt
from pathlib import Path

import numpy as np


def tracefree(x):
    return x-np.eye(3)*np.trace(x)/3


def basis(h,k):
    x=tracefree(h)
    x=x/np.linalg.norm(x)
    y=tracefree(k)
    y-=np.vdot(x,y).real*x
    return x,y/np.linalg.norm(y)


def zero_direction(pair):
    a,b=basis(*pair)
    def h(t):
        return a*cos(t)+b*sin(t)
    lo,hi=0.,pi
    p=np.linalg.det(h(lo)).real
    if abs(p)<2e-15:
        return a
    for _ in range(75):
        c=(hi+lo)/2
        v=np.linalg.det(h(c)).real
        if v*p<=0:
            hi=c
        else:
            lo=c
    z=h((lo+hi)/2)
    assert np.allclose(np.linalg.eigvalsh(z),[-1/sqrt(2),0.,1/sqrt(2)],rtol=0,atol=5e-14)
    return z


def dimension(*hs):
    gram=np.array([[np.vdot(tracefree(x),tracefree(y)).real for x in hs] for y in hs])
    return int(np.sum(np.linalg.eigvalsh(gram)>1e-10))


def plane_B(pair):
    a,b=basis(*pair)
    return np.linalg.norm(a@b-b@a)**2


def unitary(seed):
    rng=np.random.default_rng(seed)
    q,r=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
    return q@(np.diag(np.diag(r)/abs(np.diag(r))))


def align(p,q):
    zp,zq=zero_direction(p),zero_direction(q)
    _,ep=np.linalg.eigh(zp)
    _,eq=np.linalg.eigh(zq)
    u=ep@eq.conj().T
    assert np.linalg.norm(zp-u@zq@u.conj().T)<3e-13
    combined=dimension(*p,*(u@x@u.conj().T for x in q))
    assert combined<=3
    return u,combined


def edge(i,j,phase=1.):
    x=np.zeros((3,3),dtype=complex)
    x[i,j]=phase
    x[j,i]=np.conj(phase)
    return x+3*np.eye(3) # strictly positive, simple spectrum (2,3,4)


def compute():
    p=[edge(0,1),edge(0,1,1j)]
    q=[edge(0,2),edge(1,2)]
    pb,qb=plane_B(p),plane_B(q)
    assert np.allclose([pb,qb],[2,.5])
    assert dimension(*p,*q)==4
    u,combined=align(p,q)
    assert combined==3
    # A random strictly positive nondegenerate complex pair in each sector:
    # verify all four original spectra and internal transition probabilities survive.
    gen=[unitary(n) for n in range(4)]
    sample=[v@np.diag(vl)@v.conj().T for v,vl in zip(gen,[[1,3,8],[2,4,5],[1,4,11],[2,9,14]])]
    qq=sample[2:]
    matched,n=align(sample[:2],qq)
    a,b=(matched@x@matched.conj().T for x in qq)
    assert np.allclose(np.linalg.eigvalsh(a),[1,4,11]) and np.allclose(np.linalg.eigvalsh(b),[2,9,14])
    _,aa=np.linalg.eigh(a); _,bb=np.linalg.eigh(b)
    assert np.allclose(abs(aa.conj().T@bb)**2,abs(gen[2].conj().T@gen[3])**2)
    # Literal cross-sector record between internal mass eigenstates; vary only
    # the centralizer of their named common, nondegenerate anchor.
    anchor=np.diag([1.,2.,3.])
    fourier=np.exp(2j*pi*np.outer(np.arange(3),np.arange(3))/3)/sqrt(3)
    state=fourier[:,0]
    read1=fourier@np.diag([1.,2.,4.])@fourier.conj().T
    read2=fourier@np.diag([1.,3.,8.])@fourier.conj().T
    outcomes=[]
    for ph in ([0.,0.,0.],[0.,2*pi/3,4*pi/3]):
        op=np.diag(np.exp(1j*np.array(ph)))
        new=op@read2@op.conj().T
        assert np.allclose(op@anchor@op.conj().T,anchor)
        assert np.allclose(new@(op@state),op@state)
        assert np.allclose(read1@state,state)
        assert dimension(anchor,read1,anchor,new)<=3
        outcomes.append(float(abs(np.vdot(state,op@state))**2))
    assert np.allclose(outcomes,[1.,0.])
    return {"positive_pairs":{"fixed_frame_dimension":4,"optimal_free_frame_dimension":combined,
            "two_plane_B":[pb,qb],"traceless_anchor_eigenvalues":np.linalg.eigvalsh(zero_direction(p)).tolist()},
            "complex_full_spectra_example":{"aligned_span":n,"fixed_span":dimension(*sample)},
            "conditional_cross_eigenstate_probability_endpoints":outcomes}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list):
        return [rounded(v) for v in x]
    return float(format(x,".9g")) if isinstance(x,(float,np.floating)) else x


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--write",action="store_true")
    ap.add_argument("--check",action="store_true")
    out=Path(__file__).with_name("independent_results.json")
    result=rounded(compute())
    if ap.parse_args().write:
        with out.open("x",encoding="utf-8") as fp:
            json.dump(result,fp,indent=2,ensure_ascii=False)
            fp.write("\n")
    else:
        ref=json.loads(out.read_text(encoding="utf-8"))
        # The zero endpoint may differ at squared machine precision.
        assert np.allclose(result.pop("conditional_cross_eigenstate_probability_endpoints"),
                           ref.pop("conditional_cross_eigenstate_probability_endpoints"),atol=1e-24)
        assert result==ref
    print("1388 independent real-plane counts, crossed eigenstate endpoints and spectra checked.")
