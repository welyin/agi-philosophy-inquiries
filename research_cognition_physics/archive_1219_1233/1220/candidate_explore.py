"""Exploratory finite-source Hamiltonian; not acceptance evidence."""
from pathlib import Path
import math
import numpy as np

def basis(K,N,parity):
    return [(h,p,a,b) for h in (0,1) for p in range(-K,K+1)
            for a in range(N+1) for b in range(N+1-a) if (h+a+b)%2==parity]

def model(K=12,N=4,parity=0,m=25.,mass=20.,gap=.5,k=1.,kappa=2.,g=.15,rel=True):
    bs=basis(K,N,parity); ix={x:i for i,x in enumerate(bs)}; d=len(bs)
    H=np.zeros((d,d),complex); Hdiag=np.zeros(d); pobs=np.zeros(d); qobs=np.zeros(d)
    for i,(h,p,a,b) in enumerate(bs):
        pp=k*p; qq=-k*(p+a-b); mm=m+gap*h
        kin=math.sqrt(mm*mm+pp*pp)-m+math.sqrt(mass*mass+qq*qq)-mass if rel else gap*h+pp*pp/(2*mm)+qq*qq/(2*mass)
        Hdiag[i]=kin+k*(a+b)+kappa; H[i,i]=Hdiag[i]; pobs[i]=pp;qobs[i]=qq
        for dp in (-1,1):
            j=ix.get((h,p+dp,a,b))
            if j is not None:H[j,i]+=-kappa/2
        for aa,bb,dp,fac in [(a+1,b,-1,math.sqrt(a+1)),(a-1,b,1,math.sqrt(a)),(a,b+1,1,math.sqrt(b+1)),(a,b-1,-1,math.sqrt(b))]:
            j=ix.get((1-h,p+dp,aa,bb))
            if j is not None:H[j,i]+=g*fac
    assert np.linalg.norm(H-H.conj().T)<1e-10
    return bs,ix,H,pobs,qobs

def run(K=12,N=4,hin=0,kappa=2.,g=.15,T=3*math.pi/4,rel=True):
    bs,ix,H,p,q=model(K,N,hin,kappa=kappa,g=g,rel=rel)
    w,v=np.linalg.eigh(H)
    initial=np.zeros((len(bs),3),complex)
    alpha=np.array([.5,1/math.sqrt(2),.5])*np.exp(-1j*np.array([-1,0,1])*math.pi/2)
    for j,pp in enumerate((-1,0,1)):initial[ix[(hin,pp,0,0)],j]=alpha[j]
    final=v@(np.exp(-1j*T*w)[:,None]*(v.conj().T@initial))
    rho0=np.sum(initial,axis=1); psi=np.sum(final,axis=1)
    prob=abs(psi)**2; incoh=(abs(final)**2).sum(axis=1)
    event=np.array([h!=hin for h,*_ in bs])
    photonplus=np.array([a>b for h,p0,a,b in bs]); photonminus=np.array([b>a for h,p0,a,b in bs])
    qp=np.array([q0<0 for q0 in q]);qm=np.array([q0>0 for q0 in q])
    initE=float(np.vdot(rho0,H@rho0).real)
    return dict(dim=len(bs),hin=hin,kappa=kappa,rel=rel,up=float(prob@event),inc=float(incoh@event),delta=float((prob-incoh)@event),p=float(prob@p),q=float(prob@q),E=initE,
       joint_qnegative=float(prob@(event&qp)),inc_qnegative=float(incoh@(event&qp)),joint_qpositive=float(prob@(event&qm)),
       plus=float(prob@(event&photonplus)),minus=float(prob@(event&photonminus)),
       momedge=float(sum(prob[i] for i,x in enumerate(bs) if abs(x[1])==K)),fockedge=float(sum(prob[i] for i,x in enumerate(bs) if x[2]+x[3]==N)))

if __name__=='__main__':
    import json
    for hin in (0,1):
        for kap in (2.,0.):
            print(json.dumps(run(hin=hin,kappa=kap),sort_keys=True),flush=True)
        print(json.dumps(run(hin=hin,rel=False),sort_keys=True),flush=True)
