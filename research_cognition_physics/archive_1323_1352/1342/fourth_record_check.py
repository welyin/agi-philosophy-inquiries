"""Round 1342: four tetrahedral complete records and an independent last-projector test."""
import argparse
import json
from pathlib import Path
import numpy as np


def unit(v):
    return v/np.linalg.norm(v)


def tetra(north,south):
    return [north/np.sqrt(3)+south*np.sqrt(2/3)*np.exp(2j*np.pi*j/3)
            for j in range(3)]+[north]


def distance(x,y):
    return np.sqrt(max(0.,1-abs(np.vdot(x,y))**2))


def projector(x):
    return np.outer(x,x.conj())


def compute():
    rng=np.random.default_rng(134204)
    def point(n):
        return unit(rng.normal(size=n)+1j*rng.normal(size=n))
    targets= tetra(*np.eye(2,dtype=complex))
    assert np.allclose(np.abs(np.array(targets).conj()@np.array(targets).T)**2,
                       np.eye(4)*2/3+np.ones((4,4))/3)
    exact=[];approx=[];max_excess=0.;max_l1_excess=0.
    for n in [2,5,10,19]:
        half=n//2
        north=np.r_[point(half),np.zeros(n-half)]
        south=np.r_[np.zeros(half),point(n-half)]
        desired=tetra(north,south)
        energies=np.r_[2*np.pi*np.arange(half),2*np.pi*np.arange(n-half)-2*np.pi/3]
        def state(t,energy=energies):
            return desired[0]*np.exp(-1j*t*energy)
        assert all(np.allclose(state(j),desired[j]) for j in range(3))
        times=[3.,3.1,3.5,4.427,20.1]
        accept=[float(abs(np.vdot(north,state(t)))**2) for t in times]
        assert max(accept)<1/3+2e-12 and abs(accept[0]-1/3)<2e-12
        exact.append(dict(n=n,fourth_accept=accept,minimum_fourth_distance=min(distance(state(t),north) for t in times)))
        perturbed=energies+2e-3*rng.normal(size=n)
        errors=[distance(state(j,perturbed),desired[j]) for j in range(4)]
        assert errors[3]+5*max(errors[:3])>=2/3-1e-10
        approx.append(dict(n=n,prefix_error=max(errors[:3]),fourth_error=errors[3]))
        for trial in range(18):
            basis=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))[0]
            raw=np.array([abs(basis.conj().T@d)**2 for d in desired])
            weights=abs(basis.conj().T@north)**2+abs(basis.conj().T@south)**2
            initials=point(n)
            record=np.abs(basis.conj().T@initials)**2
            max_l1_excess=max(max_l1_excess,float(np.max(2*weights/3-np.sum(abs(raw-record),axis=0))))
            tv=np.sum(abs(raw-record),axis=1)/2
            max_excess=max(max_excess,1/6-np.max(tv))
            assert min(tv)>=-1e-15
    assert max_excess<1e-12 and max_l1_excess<1e-12
    # Explicit hidden-clock model. A single time-independent H exists by spectral logarithm of this unitary shift.
    orientations=[np.column_stack([v,np.array([-v[1].conjugate(),v[0].conjugate()])]) for v in targets]
    shift=np.zeros((8,8),complex)
    for j in range(4):
        k=(j+1)%4
        shift[2*k:2*k+2,2*j:2*j+2]=orientations[k]@orientations[j].conj().T
    assert np.allclose(shift.conj().T@shift,np.eye(8)) and np.allclose(np.linalg.matrix_power(shift,4),np.eye(8))
    # Each spectral projector follows from the four-term discrete Fourier transform of shift^k.
    H=np.zeros((8,8),complex)
    for k in range(4):
        angle=2*np.pi*k/4
        P=sum((np.exp(-1j*angle*j)*np.linalg.matrix_power(shift,j) for j in range(4)))/4
        H-=angle*P
    assert np.allclose(H,H.conj().T)
    energy,v=np.linalg.eigh(H)
    full0=np.r_[targets[0],np.zeros(6)]
    local_accept=[]
    for j,target in enumerate(targets):
        full=(v*np.exp(-1j*j*energy))@(v.conj().T@full0)
        local_accept.append(float(np.sum(np.abs(full.reshape(4,2)@target.conj())**2)))
    assert np.allclose(local_accept,np.ones(4))
    return dict(exact_examples=exact,nearby_examples=approx,tetra_l1_examples=72,
                max_total_variation_bound_excess=max_excess,max_single_outcome_l1_excess=max_l1_excess,
                hidden_clock_local_acceptance=local_accept,
                exact_first_three_fourth_acceptance_upper=1/3,exact_pure_distance_lower=np.sqrt(2/3),all_four_trace_lower=1/6)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for v,w in zip(a,b):compare(v,w)
    else:assert np.isclose(a,b,rtol=2e-8,atol=3e-9),(a,b)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args(); res=compute()
    path=Path(__file__).with_name('fourth_record_results.json')
    if args.check:compare(res,json.loads(path.read_text(encoding='utf8')))
    else:
        with path.open('x',encoding='utf8') as stream:
            json.dump(res,stream,indent=2,allow_nan=False);stream.write('\n')
    print('PASS',res['all_four_trace_lower'],res['exact_first_three_fourth_acceptance_upper'])
