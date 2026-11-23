"""Round 1357: direct path x EM flux x cyclic flux dynamics with both outputs."""
import argparse
import json
from pathlib import Path
import numpy as np


def energy(n,m,alternative=False):
    if alternative:
        return n*n + int(n%3!=0) + int((n+m)%3!=0)
    return n*n + int(m!=0)


def battery_witness(alternate=False):
    # The compensating work ladder values are 0,1,2. These six physical basis
    # words contain all states accessed by the preparation and the controlled S.
    def field_energy(n,m):
        if not alternate:
            return energy(n,m)
        return n*n+int((n+m)%3==0 and n%3!=0)+int((n+m)%3!=0 and n==0)
    words=[(n,m,2-field_energy(n,m)) for n in [0,1] for m in range(3)]
    assert all(0<=b<=2 for n,m,b in words)
    n0=words.index((0,0,2))
    dual=np.zeros((6,6)); transfer=np.zeros((6,6))
    for m in range(3):
        j=next(i for i,(n,x,b) in enumerate(words) if (n,x)==(0,m))
        k=next(i for i,(n,x,b) in enumerate(words) if (n,x)==(0,(m+1)%3))
        dual[k,j]=1
        ell=next(i for i,(n,x,b) in enumerate(words) if (n,x)==(1,(m-1)%3))
        transfer[ell,j]=transfer[j,ell]=1
    assert np.allclose(transfer.T@transfer,np.eye(6))
    charge=np.diag([field_energy(n,m) for n,m,b in words]); work=np.diag([b for n,m,b in words])
    assert np.array_equal(charge+work,2*np.eye(6))
    prepared=evolve(dual+dual.T,np.pi/3)[:,n0]
    for_t={}
    for tag,wait in [("work_frozen",charge),("work_free",charge+work)]:
        ps=[]
        for t in [np.pi/2,np.pi]:
            clock=evolve(wait,t)
            amp=(clock@prepared-transfer.T@clock@transfer@prepared)/2
            ps.append(float(np.vdot(amp,amp).real))
        for_t[tag]=ps
    assert np.allclose(for_t["work_free"],[0,0])
    assert np.allclose(for_t["work_frozen"],[1/9,0] if alternate else [1/3,4/9])
    return for_t


def evolve(gener,phase):
    val,vec=np.linalg.eigh(gener)
    return (vec*np.exp(-1j*phase*val))@vec.conj().T


def run():
    # Exact finite-support computation. A 9-site cyclic grid supplies unitary
    # matrices outside the visited EM positions 0 and 1. Boundary is never used.
    flux=list(range(-4,5))
    labels=[(n,m) for n in flux for m in range(3)]
    dim=len(labels)
    def ket(n,m):
        x=np.zeros(dim,complex);x[labels.index((n,m))]=1
        return x
    S=np.zeros((dim,dim))
    D=S.copy()
    remap=S.copy()
    for i,(n,m) in enumerate(labels):
        S[labels.index((flux[(flux.index(n)+1)%len(flux)],(m-1)%3)),i]=1
        D[labels.index((n,(m+1)%3)),i]=1
        remap[labels.index((n,(n+m)%3)),i]=1
    assert np.allclose(remap@remap.T,np.eye(dim))
    assert np.allclose(remap@D@remap.T,D)
    S_mem=np.kron(np.roll(np.eye(len(flux)),1,axis=0),np.eye(3))
    assert np.allclose(remap@S@remap.T,S_mem)
    sector=np.diag([np.exp(2j*np.pi*(n+m)/3) for n,m in labels])
    assert np.allclose(sector@S,S@sector)
    assert np.allclose(sector@D,np.exp(2j*np.pi/3)*D@sector)
    anchor=ket(0,0)
    pulse=evolve(D+D.T,np.pi/3)
    pulsed=pulse@anchor
    assert np.allclose(abs(pulsed.reshape(-1,3)[flux.index(0)])**2,[1/9,4/9,4/9])
    had=np.kron(np.ones((2,2))-np.diag([0,2]),np.eye(dim))/np.sqrt(2)
    ctrl=np.block([[np.eye(dim),np.zeros((dim,dim))],[np.zeros((dim,dim)),S]])
    assert np.allclose(pulse@pulse,np.exp(2j*np.pi/3)*np.eye(dim))
    dephased=np.diag(abs(pulsed)**2)
    after_second=pulse@dephased@pulse.conj().T
    dephased_pop=np.diag(after_second).real.reshape(len(flux),3)[flux.index(0)]
    assert np.allclose(dephased_pop,np.array([11,8,8])/27)
    echo_checks=[]
    branches={}
    for alt in [False,True]:
        H=np.diag([energy(n,m,alt) for n,m in labels])
        # Under the coordinate change d=n+m (mod 3) the whole finite algebra is a rotor plus a neutral memory.
        aliasH=np.diag([energy(n,(d-n)%3,alt) for n,d in labels])
        assert np.allclose(remap@H@remap.T,aliasH)
        rows={}
        for name,v in [("no_pulse",anchor),("field_pulse",pulsed)]:
            assert np.isclose(np.vdot(v,H@v).real,0 if name=="no_pulse" else 8/9)
            start=had@np.r_[v,np.zeros(dim)]
            shifted=ctrl@start
            prob_n=np.sum(abs(shifted.reshape(2,len(flux),3))**2,axis=(0,2))
            assert np.allclose(prob_n[[n not in [0,1] for n in flux]],0.)
            balanced_work=np.vdot(shifted,np.kron(np.eye(2),H)@shifted).real-np.vdot(v,H@v).real
            readings={}
            for label,t in [("pi_over_2",np.pi/2),("pi",np.pi)]:
                proc=had@ctrl.T@np.kron(np.eye(2),evolve(H,t))@ctrl@had
                final=proc@np.r_[v,np.zeros(dim)]
                branches_at_t=final.reshape(2,len(flux),3)
                support=branches_at_t[:,flux.index(0),:]
                prob=abs(support)**2
                assert np.isclose(prob.sum(),1.)
                assert np.allclose(prob.sum(axis=0),abs(v.reshape(-1,3)[flux.index(0)])**2)
                dark=float(prob[1].sum())
                # Dephasing all magnetic branches before the path gives the same classical joint record.
                weights=abs(v.reshape(-1,3)[flux.index(0)])**2
                mixture=sum(weights[m]*abs((proc@np.r_[ket(0,m),np.zeros(dim)]).reshape(2,len(flux),3)[:,flux.index(0),:])**2 for m in range(3))
                assert np.allclose(mixture,prob)
                if not alt and name=="no_pulse":
                    echo_mixed=sum(dephased_pop[m]*np.linalg.norm((proc@np.r_[ket(0,m),np.zeros(dim)])[dim:])**2 for m in range(3))
                    echo_checks.append(float(echo_mixed))
                if alt or name=="no_pulse":
                    target=1. if label=="pi_over_2" else 0.
                else:
                    target=1/3 if label=="pi_over_2" else 4/9
                assert np.isclose(dark,target)
                conditional=None if dark<1e-13 else (prob[1]/dark).tolist()
                readings[label]=dict(dark=dark,bright=float(prob[0].sum()),dark_conditional_flux=conditional)
            rows[name]=dict(mean_controlled_shift_work_in_kappa=float(balanced_work),observations=readings)
        branches["alternate" if alt else "register"]=rows
    def shadow(n,m):
        return n*n+int(n%3!=0)
    assert all(energy(n,-n%3)==shadow(n,-n%3)==energy(n,-n%3,True) for n in flux)
    assert all(energy(0,m)==energy(0,m,True) for m in range(3))
    assert np.allclose(echo_checks,[5/9,8/27])
    return dict(branches=branches, double_pulse_with_intermediate_dephasing=echo_checks,
                same_three_level_battery={"registered":battery_witness(),"alternative":battery_witness(True)},
                after_pulse_probabilities=[1/9,4/9,4/9],
                injected_energy_in_kappa=8/9, single_frequency_robust_absolute_error=(np.sqrt(113)-7)/24)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif a is None:
        assert b is None
    else:
        assert np.allclose(a,b,atol=3e-12,rtol=3e-12),(a,b)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    r=run()
    out=Path(__file__).with_name("independent_results.json")
    if args.check:
        compare(r,json.loads(out.read_text(encoding="utf-8")))
    else:
        with out.open("x",encoding="utf-8") as f:
            json.dump(r,f,indent=2,allow_nan=False);f.write("\n")
    print("PASS",r["branches"]["register"]["field_pulse"])
