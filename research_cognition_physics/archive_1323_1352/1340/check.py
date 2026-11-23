"""Round 1340 finite parent. Uses only NumPy and the standard library.
Default replays the saved results without writing. --write generates results.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def dagger(a):
    return a.conj().T


def kron_many(factors):
    out=np.array([1.])
    for factor in factors:
        out=np.kron(out,factor)
    return out


def flavor_exchange(n, bank_count, index):
    axes=list(range(1+bank_count))
    axes[0],axes[index+1]=axes[index+1],axes[0]
    shape=(n,)*(1+bank_count)
    perm=np.arange(n**(1+bank_count)).reshape(shape).transpose(axes).ravel()
    return np.eye(len(perm))[perm]


def prepare(directions, copies):
    bank=kron_many([r for r,N in zip(directions,copies) for _ in range(N)])
    return np.kron(np.eye(len(directions[0])),bank[:,None])


def project(r):
    return np.outer(r,r.conj())


def parent_and_prediction(directions,copies,weights):
    n=len(directions[0])
    dim=n**(1+sum(copies))
    full=np.eye(dim)
    avg=np.eye(n,dtype=complex)
    variance=np.zeros((n,n),dtype=complex)
    start=0
    for r,N,w in zip(directions,copies,weights):
        if w:
            for j in range(start,start+N):
                full+=w/N*flavor_exchange(n,sum(copies),j)
            avg+=w*project(r)
            variance+=w*w/N*(np.eye(n)-project(r))
        start+=N
    return full,avg,variance


def trig(a,func,angle):
    eig,v=np.linalg.eigh(a)
    return (v*func(eig*angle))@dagger(v)


class Checks:
    def __init__(self):
        self.count={}
        self.error={}
    def close(self,key,value,target):
        diff=np.max(np.abs(np.asarray(value)-np.asarray(target)))
        e=float(diff/max(1.,np.max(np.abs(target))))
        if not np.isfinite(e) or e>4e-10:
            raise AssertionError((key,e))
        self.count[key]=self.count.get(key,0)+1
        self.error[key]=max(self.error.get(key,0.),e)
    def bound(self,key,value,target):
        self.close(key,max(value-target,0.),0.)


def rot(a,angle):
    vals,vec=np.linalg.eigh(a)
    return (vec*np.exp(-1j*angle*vals))@dagger(vec)


def study(log,directions,copies,rows,tag):
    n=len(directions[0])
    J=prepare(directions,copies)
    tu,td=.2,.25
    times=(.45,.65)
    matrix=[]
    for w,t,T in zip(rows,(tu,td),times):
        full,avg,var=parent_and_prediction(directions,copies,w)
        matrix.append((full,avg))
        diff=full@J-J@avg
        log.close("field_source_compression",dagger(J)@full@J,avg)
        log.close("retained_reference_variance",dagger(diff)@diff,var)
        eps=t*math.sqrt(sum(x*x/N for x,N in zip(w,copies)))
        used=next(i for i,v in enumerate(w) if v)
        offset=sum(copies[:used])
        term=sum(flavor_exchange(n,sum(copies),i) for i in range(offset,offset+copies[used]))/copies[used]
        for sign in (-1,1):
            # X-chirality sectors. The equality comparison includes the banks.
            actual=rot(full,sign*t*T)@J
            ideal=J@rot(avg,sign*t*T)
            log.bound("uniform_joint_state",np.linalg.norm(actual-ideal,2),eps*T)
            term_actual=dagger(actual)@term@actual
            term_want=dagger(ideal)@term@ideal
            log.bound("single_bank_source_response",np.linalg.norm(term_actual-term_want,2),2*eps*T)
            src=dagger(actual)@full@actual
            want=dagger(ideal)@full@ideal
            log.bound("whole_source_observable",np.linalg.norm(src-want,2),
                      2*(1+sum(abs(x) for x in w))*eps*T)
        ev,ef=np.linalg.eigh(t*full)
        ey,vy=np.linalg.eigh(t*avg)
        for i,e in enumerate(ey):
            psi=J@vy[:,i]
            prob=np.abs(dagger(ef)@psi)**2
            measured=np.dot(prob,(ev-e)**2)
            allowed=t*t*np.vdot(vy[:,i],var@vy[:,i]).real
            log.close("spectral_measure_variance",measured,allowed)
            window=.42
            log.bound("spectral_window",sum(prob[np.abs(ev-e)>=window]),allowed/window**2)
    up,down=matrix
    r=directions[0]
    data=J@r
    pre=trig(up[0],np.cos,tu*times[0])@data
    post=trig(down[0],np.sin,td*times[1])@pre
    ideal_pre=trig(up[1],np.cos,tu*times[0])@r
    ideal_post=trig(down[1],np.sin,td*times[1])@ideal_pre
    exact=float(np.vdot(post,post).real)
    predicted=float(np.vdot(ideal_post,ideal_post).real)
    other_chirality=trig(down[0],np.cos,td*times[1])@pre
    failed_source=trig(up[0],np.sin,tu*times[0])@data
    log.close("all_weak_records",exact+np.linalg.norm(other_chirality)**2+np.linalg.norm(failed_source)**2,1.)
    chain_error=sum(t*T*math.sqrt(sum(wj*wj/N for wj,N in zip(w,copies)))
                    for (t,T,w) in zip((tu,td),times,rows))
    # The u-L -> d-L successful source branch follows precisely these sectors.
    # Apply the norm-to-trace bound to the full source isometry with worst
    # sector errors; branch norm alone gives the weaker 2*chain_error.
    generic=(tu*times[0]*math.sqrt(sum(wj*wj/N for wj,N in zip(rows[0],copies)))+
             max(tu*math.sqrt(sum(wj*wj/N for wj,N in zip(rows[0],copies))),
                 td*math.sqrt(sum(wj*wj/N for wj,N in zip(rows[1],copies))))*times[1])
    log.bound("weak_joint_record_error",abs(exact-predicted),generic)
    return dict(tag=tag,N=copies,exact_joint=exact,effective_joint=predicted,
                sector_chain_norm_bound=chain_error,probability_bound=generic,
                up_effective_eigenvalues=np.linalg.eigvalsh(up[1]).tolist(),
                down_effective_eigenvalues=np.linalg.eigvalsh(down[1]).tolist())


def weak_source(log):
    # f=0,1 for up/down, ch=0,1 for L/R, s=0,1 for the source components.
    V=np.zeros((8,8))
    for f in range(2):
        for ch in range(2):
            for s in range(2):
                col=4*f+2*ch+s
                dest=(s,ch,f) if ch==0 else (f,ch,s)
                V[4*dest[0]+2*dest[1]+dest[2],col]=1.
    log.close("source_unitary",dagger(V)@V,np.eye(8))
    for sigma in (np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])):
        total=(np.kron(np.kron(sigma,np.diag([1,0])),np.eye(2))
               +np.kron(np.eye(4),sigma))/2
        log.close("weak_scalar_pulse",V@total,total@V)
    k=V.reshape(4,2,4,2)[:,0,:,1]
    expected=np.zeros((4,4))
    expected[2,0]=1
    log.close("selected_source_kraus",k,expected)


def preparation_bound(N,eps_fail,independent):
    n=3
    groups=(2 if independent else 1, 1)
    attempts=[math.ceil(math.log(sum(groups)*N/eps_fail)/(-math.log1p(-n**(-2*j))))
              for j in (1,2)]
    union=N*sum(k*(1-n**(-2*j))**L for j,L,k in zip((1,2),attempts,groups))
    l1,l2=attempts
    return dict(copies_each=N,live_bank_registers=N*(1+sum(groups)),attempt_caps=attempts,
                registered_seed_upper=N*(1+groups[0]*l1+l2),
                registered_Choi_upper=N*(groups[0]*l1+2*l2),formation_failure_upper=union)


def calculate():
    log=Checks()
    r=np.sqrt([.5,.3,.2])
    U=np.diag([1,1j,-1])
    r1,r2=U@r,U@U@r
    Z=np.diag([1,-1,1])
    shared=study(log,[r,r1,r2],[1,1,1],[[1,2,0],[0,1,3]],"common")
    rival=study(log,[r,r1,Z@r1,Z@r2],[1,1,1,1],[[1,2,0,0],[0,0,1,3]],"independent")
    log.close("equal_separate_spectra",shared["down_effective_eigenvalues"],rival["down_effective_eigenvalues"])
    # Additional size case with two copies per bank on a different flavor space.
    x=np.array([1.,0.])
    y=np.array([1.,1j])/np.sqrt(2)
    z=np.array([1.,-1.])/np.sqrt(2)
    paired=study(log,[x,y,z],[2,2,2],[[1,2,0],[0,1,3]],"two_per_bank")
    weak_source(log)
    worst=.25*math.sqrt(10)
    N=math.ceil(((.45+.65)*worst/.1)**2)
    common_budget=preparation_bound(N,.01,False)
    rival_budget=preparation_bound(N,.01,True)
    gap=abs(shared["effective_joint"]-rival["effective_joint"])
    if gap<=2*(.1+.01) or abs(shared["exact_joint"]-rival["exact_joint"])<.01:
        raise AssertionError("missing difference in finite predictions")
    return dict(round=1340,status="PASS",shared=shared,rival=rival,two_per_bank=paired,
                finite_budget=dict(shared=common_budget,rival=rival_budget,uniform_state_bound=(.45+.65)*worst/math.sqrt(N),
                                   certified_joint_probability_gap=gap-2*(.1+.01)),
                checks=log.count,worst_errors=log.error)


def compare(x,y):
    if isinstance(y,dict):
        if x.keys()!=y.keys():
            raise AssertionError("keys differ")
        for k in y:
            compare(x[k],y[k])
    elif isinstance(y,(float,list)):
        np.testing.assert_allclose(x,y,rtol=2e-8,atol=3e-10)
    elif x!=y:
        raise AssertionError((x,y))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--write",action="store_true")
    args=p.parse_args()
    f=Path(__file__).with_name("results.json")
    result=calculate()
    if args.write:
        f.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    else:
        compare(result,json.loads(f.read_text(encoding="utf-8")))
    print(json.dumps(dict(status=result["status"],checks=sum(result["checks"].values()),
                          exact=[result[k]["exact_joint"] for k in ("shared","rival")],
                          effective=[result[k]["effective_joint"] for k in ("shared","rival")],
                          budget=result["finite_budget"])))
