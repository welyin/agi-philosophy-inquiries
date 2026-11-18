"""Independent tensor-level audit for round 1338. Only NumPy and stdlib."""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


def kron_all(factors):
    product=np.array([1.])
    for factor in factors:
        product=np.kron(product,factor)
    return product


def perm(n,order):
    addresses=np.arange(n**len(order)).reshape((n,)*len(order))
    images=addresses.transpose(order).ravel()
    return np.eye(n**len(order))[images]


def ket(n,place):
    return np.eye(n)[:,place]


def unitary_parent(n,count):
    """Diagonalize a Hermitian combination of permutation intertwiners."""
    identity=list(range(count))
    swap=identity.copy()
    swap[0],swap[1]=swap[1],swap[0]
    p=perm(n,swap)
    cycle=perm(n,identity[1:]+[0])
    h=.43*p+.38*(cycle+cycle.T)+.19j*(cycle-cycle.T)
    eigen,v=np.linalg.eigh(h)
    return (v*np.exp(-1j*eigen))@v.conj().T


def cg_column(n,count,r,slot):
    """Coordinates of an arbitrary input at slot, background r elsewhere."""
    return np.column_stack([kron_all([ket(n,j) if t==slot else r for t in range(count)])
                            for j in range(n)])


def random_unit(rng,n):
    x=rng.normal(size=n)+1j*rng.normal(size=n)
    return x/np.linalg.norm(x)


def mutual_entanglement(phi,n,count):
    if count<3:
        return None
    return np.linalg.svd(phi.reshape(n,-1),compute_uv=False)[1]


def compute():
    rng=np.random.default_rng(133804)
    stats=dict(tensor_cases=0,total_readouts=0,entangled_readouts=0,
               max_covariance_error=0.,max_sector_amplitude_error=0.,
               max_rank_one_error=0.,max_middle_error=0.,max_leading_ratio_error=0.,
               max_complete_bound_excess=0.,max_strictly_middle_mass_spread=0.)
    rows=[]
    for n,count in [(3,2),(3,3),(3,4),(4,2),(4,3)]:
        u=unitary_parent(n,count)
        # A fixed generic rotation via QR, for direct common-rotation covariance.
        v,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
        rot=kron_all([v]*count)
        err=np.linalg.norm(u@rot-rot@u)
        stats['max_covariance_error']=max(stats['max_covariance_error'],err)
        for r,pick in [(ket(n,0),ket(n,1)),(v[:,0],v[:,1])]:
            w=u@cg_column(n,count,r,0)
            zero=kron_all([r]*count)
            phase=np.vdot(zero,u@zero)
            slots=np.column_stack([cg_column(n,count,r,j)@pick for j in range(count)])
            coeff=slots.conj().T@(w@pick)
            a=coeff[0]
            err=max(np.linalg.norm(w@pick-slots@coeff),abs(np.vdot(coeff,coeff)-1),
                    abs(sum(coeff)-phase),abs(abs(phase)-1))
            stats['max_sector_amplitude_error']=max(stats['max_sector_amplitude_error'],err)
            complement=np.linalg.qr(w,mode='complete')[0][:,n:]
            data=np.eye(n)
            phi_list=[kron_all([r]*(count-1)),random_unit(rng,n**(count-1))]
            if count>=3:
                # Incident all-r amplitude mixed coherently with an orthogonal joint branch.
                phi_list.append((kron_all([r]*(count-1))+1j*kron_all([pick]*(count-1)))/np.sqrt(2))
            for phi in phi_list:
                p=np.kron(data,phi[:,None])
                l=p.conj().T@w
                gamma=np.vdot(phi,kron_all([r]*(count-1)))
                b=abs(a*gamma)
                proj=np.outer(r,r.conj())
                error=max(np.linalg.norm((data-proj)@(l-a*gamma*data)),
                          np.linalg.norm(l@r-phase*gamma*r))
                sv=np.linalg.svd(l,compute_uv=False)[::-1]
                assert b>1e-6
                curv=sv[0]*sv[-1]/sv[1]**2
                stats['max_rank_one_error']=max(stats['max_rank_one_error'],error)
                stats['max_middle_error']=max(stats['max_middle_error'],max(abs(sv[1:-1]-b)))
                stats['max_leading_ratio_error']=max(stats['max_leading_ratio_error'],abs(curv-1/abs(a)))
                ent=mutual_entanglement(phi,n,count)
                if ent is not None and ent>1e-8:
                    stats['entangled_readouts']+=1
                for t,m in [(1.,3.),(.9,.8)]:
                    full=np.column_stack([t*p,m*complement])
                    low=np.linalg.svd(full,compute_uv=False)[::-1][:n]
                    curvature=low[0]*low[-1]/low[1]**2
                    excess=max(0.,1/abs(a)-curvature,curvature-(1+min(t/m,m/t)**2)/abs(a))
                    stats['max_complete_bound_excess']=max(stats['max_complete_bound_excess'],excess)
                    stats['max_strictly_middle_mass_spread']=max(stats['max_strictly_middle_mass_spread'],float(np.ptp(low[1:-1])))
                    stats['tensor_cases']+=1
                stats['total_readouts']+=1
            rows.append(dict(n=n,N=count,mod_a=abs(a),phase=[phase.real,phase.imag]))
    for key,value in stats.items():
        if key.startswith('max_'):
            assert value<2e-9,(key,value)
    assert stats['entangled_readouts']>0

    n,count=3,3
    r,s=ket(n,0),ket(n,1)
    p12=perm(n,[1,0,2]); p13=perm(n,[2,1,0])
    ident=np.eye(n**3)
    c1,c2=.75,.8
    u=(c2*ident-1j*math.sqrt(1-c2*c2)*p13)@(c1*ident-1j*math.sqrt(1-c1*c1)*p12)
    init=np.kron(np.eye(n),np.kron(r,s)[:,None])
    l=init.conj().T@u@init
    sv=np.linalg.svd(l,compute_uv=False)[::-1]
    curv=sv[0]*sv[-1]/sv[1]**2
    assert np.allclose(sv,[3/5,3/4,4/5],rtol=0,atol=1e-13)
    return dict(stats=stats,single_excitation_rows=rows,unequal_reference_boundary=dict(spectrum=sv.tolist(),curvature=curv))


def compare(a,b):
    if isinstance(a,dict):
        assert set(a)==set(b)
        for key in a:
            compare(a[key],b[key])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for first,second in zip(a,b):
            compare(first,second)
    else:
        assert np.isclose(a,b,rtol=2e-8,atol=5e-10),(a,b)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=compute()
    output=Path(__file__).with_name('joint_port_independent_results.json')
    if args.check:
        compare(result,json.loads(output.read_text(encoding='utf8')))
    else:
        with output.open('x',encoding='utf8') as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write('\n')
    print(json.dumps(result['stats'],allow_nan=False))
