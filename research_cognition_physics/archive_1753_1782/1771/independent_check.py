"""Gauge-index projections and the least possible joint d=5 tensor on split flavour embeddings."""
import argparse
import itertools
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ALPHA, BETA = 2/3, 1/2


def physical_tensors(y,left,right):
    l = np.einsum('ij,kl->ijkl',y,left)
    r = np.einsum('ij,kl->ijkl',y,right)
    sym = (l+l.transpose(1,2,0,3)+l.transpose(2,0,1,3))/3
    return l-sym,(r-r.swapaxes(0,2))/2


def gauge_check():
    rng = np.random.default_rng(37)
    def draws(*shape):
        return rng.normal(size=shape)+1j*rng.normal(size=shape)
    y = draws(3,3)
    y = (y+y.T)/2
    a,b = draws(3,3),draws(3,3)
    q,l,u,d,e = draws(3,3,2),draws(3,2),draws(3,3),draws(3,3),draws(3)
    eps3 = np.zeros((3,3,3))
    for perm in itertools.permutations(range(3)):
        eps3[perm] = np.linalg.det(np.eye(3)[list(perm)])
    eps2 = np.array([[0,1],[-1,0]])
    op_l = np.einsum('abc,AB,CD,iaA,jbB,kcC,lD->ijkl',eps3,eps2,eps2,q,q,q,l,optimize=True)
    op_r = np.einsum('abc,ia,kb,lc,j->ijkl',eps3,u,u,d,e,optimize=True)
    k = np.linalg.norm(y)**2*np.eye(3)-y@y.conj().T
    ll,rr = physical_tensors(y,a,b)
    for projected,down,op,weight in [(ll,a,op_l,ALPHA),(rr,b,op_r,BETA)]:
        raw = np.einsum('ij,kl->ijkl',y,down)
        np.testing.assert_allclose(np.sum(op*(raw-projected)),0.,atol=5e-12)
        np.testing.assert_allclose(np.linalg.norm(projected)**2,
                                   weight*np.trace(down.conj().T@k@down).real)


def minimum(y,el,down,pl,pd):
    # A single Higgs-index insertion can shift the two triplet vertices by
    # the restriction of one unrestricted parent matrix to PL and PD.
    frame = np.column_stack((np.sqrt(ALPHA)*pl,np.sqrt(BETA)*pd))
    data = np.column_stack((np.sqrt(ALPHA)*el,np.sqrt(BETA)*down))
    shift = -data@np.linalg.pinv(frame)
    ll,rr = physical_tensors(y,el+shift@pl,down+shift@pd)
    squared = np.linalg.norm(ll)**2+np.linalg.norm(rr)**2
    x,sv,yh = np.linalg.svd(pl.conj().T@pd)
    shared = np.isclose(sv,1.,rtol=0,atol=2e-12)
    left, right = x[:,shared], yh.conj().T[:,shared]
    np.testing.assert_allclose(pl@left,pd@right,atol=1e-11)
    mismatch = el@left-down@right
    k = np.linalg.norm(y)**2*np.eye(3)-y@y.conj().T
    formula = (2/7)*np.trace(mismatch.conj().T@k@mismatch).real
    np.testing.assert_allclose(squared,formula,rtol=2e-8,atol=1e-24)
    return np.sqrt(squared),int(sum(shared))


def solve():
    gauge_check()
    previous = runpy.run_path(str(STAGE/'1768/check.py'))
    mat = previous['MAT']
    fits = json.loads((STAGE/'1766/independent_results.json').read_bytes())['refits']
    out = {}
    for name,row in fits.items():
        target = np.array([row['moduli'][k] for k in ('Vus','Vub','Vcb')])
        trial = target.copy()
        for _ in range(3):
            actual = mat['eval_point'](np.log(trial),name,.25)
            v = np.array([actual['moduli'][k] for k in ('Vus','Vub','Vcb')])
            if max(abs(v-target))<2e-10:
                break
            trial += target-v
        p = mat['load_mother'](*trial,.25)[0]['pt']
        d,e,y = p['D'], p['E'], np.diag(p['A'])
        q,g,f,w,_ = previous['tagged_frame'](row,name)
        ordinary_l = np.eye(3)+(w*(f-1))@w.conj().T
        pl = np.vstack((ordinary_l,np.sqrt(1-f*f)[:,None]*w.conj().T))
        ordinary_d = q*g
        j = np.flatnonzero(g<1)
        heavy = np.zeros((2,3))
        heavy[np.arange(2),j] = np.sqrt(1-g[j]**2)
        vals, dims = [], []
        for u in (np.eye(2),np.diag([1j,-1]),np.array([[0,1],[-1,0]]),
                  np.array([[1,1j],[1j,1]])/math.sqrt(2)):
            pd = np.vstack((ordinary_d,u@heavy))
            for embedding in [pl,pd]:
                np.testing.assert_allclose(embedding.conj().T@embedding,np.eye(3),atol=2e-14)
            val,dim = minimum(y,e@ordinary_l,d@ordinary_d,pl,pd)
            vals.append(val)
            dims.append(dim)
        k = np.linalg.norm(y)**2*np.eye(3)-y@y.conj().T
        kroot = np.diag(np.sqrt(k.diagonal().real))
        lo = math.sqrt(2/7)*np.linalg.svd(kroot@(e-d),compute_uv=False)[-1]*max(min(g),min(f))
        assert lo>0 and min(vals)>lo
        aligned,three = minimum(y,e,d,np.eye(3),np.eye(3))
        # With a third heavy direction the two three-planes may be disjoint.
        ext_l = np.vstack((np.eye(3),np.zeros((3,3))))
        ext_d = np.vstack((np.cos(.25)*np.eye(3),np.sin(.25)*np.eye(3)))
        escaped,zero = minimum(y,e,d,ext_l,ext_d)
        assert three==3 and zero==0 and escaped<1e-15
        out[name]=dict(common_plane_norm=aligned, four_relative_frames_norms=vals,
                       common_dimensions=dims, uniform_lower_bound=lo, extra_direction_norm=escaped)
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=json.loads(json.dumps(solve()),parse_float=lambda v:float(f'{float(v):.10g}'))
    path=HERE/'independent_results.json'
    if args.check:
        saved=json.loads(path.read_bytes())
        for name,row in result.items():
            assert set(saved[name])==set(row)
            for k,v in row.items():
                np.testing.assert_allclose(v,saved[name][k],rtol=1e-8,atol=1e-17)
        print("1771 independent: gauge-irreducible d=5 coefficients, split-plane theorem, six-family assignments")
    else:
        with path.open('x',encoding='utf8') as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write('\n')
