"""Inspect independent basis-invariant one-loop degeneracies and spectral thresholds."""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent


def fourier(n):
    return np.exp(2j*np.pi*np.outer(np.arange(n),np.arange(n))/n)/np.sqrt(n)


def basis(n,seed):
    rng=np.random.default_rng(seed)
    q,r=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
    return q @ np.diag(np.diag(r)/np.abs(np.diag(r)))


def pmap(tensors,which,scale):
    # Build dephasing maps in Hermitian-matrix space, extended complex linearly.
    def linear(x):
        if which=="up":
            return sum(t@x.T@t.conj().T for t in tensors)/scale**2
        tmp=sum(t.conj().T@x@t for t in tensors)/scale**2
        return sum(t@tmp@t.conj().T for t in tensors)/scale**2
    return np.column_stack([linear(x.reshape(3,3)).reshape(-1) for x in np.eye(9)])


def witness(relative):
    p=basis(3,781)
    q=p@relative
    r=basis(3,1980)
    u,d=np.sqrt(8/5),np.sqrt(6/5)
    ups=np.array([u*np.outer(a,a) for a in p.T]+[np.zeros((3,3))])
    downs=np.array([d*np.outer(a,b.conj()) for a,b in zip(q.T,r.T)]+[np.zeros((3,3))])
    su,sd=np.array([1.e-5,.006,.79]),np.array([.000053,.0017,.42])
    def kernel(s,scale):
        v=s/scale
        return np.r_[v,np.sqrt(1-v@v)]
    hu,hd=kernel(su,u),kernel(sd,d)
    yu=np.einsum("a,aij->ij",hu,ups)
    yd=np.einsum("a,aij->ij",hd,downs)
    assert np.allclose(np.linalg.svd(yu,compute_uv=False)[::-1],su)
    assert np.allclose(np.linalg.svd(yd,compute_uv=False)[::-1],sd)
    assert np.allclose(p.conj().T@yu@p.conj(),np.diag(su))
    assert np.allclose(q.conj().T@yd@r,np.diag(sd))
    assert np.allclose(np.linalg.svd(yd.T,compute_uv=False),sd[::-1])
    residual=0.
    for arr,k in [(ups,u),(downs,d)]:
        for terms in [[a@a.conj().T for a in arr],[a.conj().T@a for a in arr]]:
            residual=max(residual,float(np.linalg.norm(sum(terms)-k*k*np.eye(3))))
        gram=np.einsum("aij,bij->ab",arr.conj(),arr)
        residual=max(residual,float(np.linalg.norm(gram-k*k*np.diag([1,1,1,0]))))
    assert residual<1.e-12
    pu=pmap(ups,"up",u)
    pd=pmap(downs,"down",d)
    for t in (pu,pd):
        assert np.allclose(t.conj().T,t) and np.allclose(t@t,t)
    invariant=float(np.trace(pu@pd).real)
    assert abs(invariant-np.sum(abs(relative)**4))<1.e-12
    for arr,fn,k,mp,seed in [(ups,"up",u,pu,888),(downs,"down",d,pd,340)]:
        rotated=np.einsum("ab,bij->aij",basis(4,seed),arr)
        assert np.allclose(pmap(rotated,fn,k),mp)
    j=float((relative[0,0]*relative[1,1]*(relative[0,1]*relative[1,0]).conjugate()).imag)
    return dict(higgs_inclusive_overlap=invariant,light_quartet=j,
                one_loop_family_max_residual=residual,projected_svd_u=list(su),projected_svd_d=list(sd))


def run():
    bounds=[]
    for k in (1,2):
        md,ms,le=map(F,("4.675","92.47","206.76838"))
        ed,es,el=map(F,(".056",".69",".00017"))
        low=(md-k*ed)*(le-k*el)/(ms+k*es)
        high=(md+k*ed)*(le+k*el)/(ms-k*es)
        ratio=float(low/9)
        frac=(np.sqrt(ratio)-1)/(np.sqrt(ratio)+1)
        bounds.append(dict(widths=k,D_low=float(low),D_high=float(high),
                           GJ_min_double_ratio_fraction=ratio-1,
                           GJ_min_equal_four_mass_fraction=frac))
    return dict(relative_orientation={name:witness(mat) for name,mat in
                [("aligned",np.eye(3)),("fourier",fourier(3)),("random",basis(3,411))]},direct_mass_checks=bounds)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    grp=parser.add_mutually_exclusive_group()
    grp.add_argument("--write",action="store_true")
    grp.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    output=ROOT/"independent_results.json"
    if args.write:
        with output.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    else:
        saved=json.loads(output.read_text(encoding="utf-8-sig"))
        def same(a,b):
            if isinstance(a,dict):
                return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
            if isinstance(a,list):
                return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
            return np.isclose(a,b,rtol=2e-10,atol=2e-13)
        assert same(result,saved)
        print("1480 independent inclusive tensor and mass checks passed")
