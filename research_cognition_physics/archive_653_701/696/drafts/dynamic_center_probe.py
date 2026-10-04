"""Exploratory full S9 integral on a dynamic weak-centre spatial-loop background.
Fixed original local objects; one spatial site, two AP time sites, all16 channels.
"""
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import joint_gauss_boundary_functional as b

def background(first=True,second=False):
    r=b.prior.rep(np.eye(3),-np.eye(2),1)
    n=128;gs=np.kron(np.kron(np.eye(2),b.internal.spin.G5),np.eye(16))
    g0=np.kron(np.kron(np.eye(2),b.internal.spin.GAMMA[3]),np.eye(16))
    ts=np.zeros((n,n),complex)
    ts[:64,:64]=np.kron(np.eye(4),r if first else np.eye(16))
    ts[64:,64:]=np.kron(np.eye(4),r if second else np.eye(16))
    tt=np.kron(np.array([[0,1],[-1,0]]),np.eye(64))
    x=np.eye(n)-(ts+ts.conj().T)/2+g0@tt
    h=gs@x;ev,vec=np.linalg.eigh(h);d=(np.eye(n)+gs@(vec*np.sign(ev))@vec.conj().T)/2
    u,v=vec[:,ev<0],vec[:,ev>0]
    f=np.kron(np.kron(np.eye(2),np.array([[1,0],[1,0],[0,1],[0,1]])/np.sqrt(2)),np.eye(16))
    k=np.kron(np.kron(np.eye(2),b.internal.spin.G5@b.internal.spin.GAMMA[1]),np.eye(16))
    lam,w=np.linalg.eigh(f.T@h@f);un,vn=w[:,lam<0],w[:,lam>0]
    coeff=np.zeros((2,10,32,32),complex)
    for site in range(2):
        for a in range(10):
            m=np.zeros((128,128),complex);m[64*site:64*(site+1),64*site:64*(site+1)]=np.kron(b.internal.B,b.internal.T[a])
            coeff[site,a]=un.T@f.T@m@k@f@vn
    o=b.prior.vector_rotation(r);pc=np.where(np.diag(o)>.5)[0];pw=np.where(np.diag(o)<-.5)[0]
    assert len(pc)==6 and len(pw)==4 and np.max(abs(o-np.diag(np.diag(o))))<1e-12
    rng=np.random.default_rng(69601);factor=None;ratios=[]
    for j in range(3):
        e=rng.normal(size=(2,10));e/=np.linalg.norm(e,axis=1)[:,None]
        c=np.einsum('xa,xaij->ij',e,coeff);det=np.linalg.det(c)
        full=b.regular(u,v,d,b.fixed_matrices(e,b.mass.car.PHI[:2]),lam=0)['weight']
        if factor is None:factor=full/det
        ratios.append(float(abs(full/(factor*det)-1)))
    return coeff,pc,pw,factor,dict(Wilson_gap=float(min(abs(ev))),factor=[float(factor.real),float(factor.imag)],crosscheck=ratios)

def jacobi(n,alpha,beta):
    k=np.arange(n,dtype=float);s=2*k+alpha+beta
    diag=(beta*beta-alpha*alpha)/(s*(s+2))
    k=np.arange(1,n,dtype=float);s=2*k+alpha+beta
    off=2/s*np.sqrt(k*(k+alpha)*(k+beta)*(k+alpha+beta)/((s-1)*(s+1)))
    x,v=np.linalg.eigh(np.diag(diag)+np.diag(off,1)+np.diag(off,-1))
    return x,v[0]**2

def integral(coeff,pc,pw,nr=9,na=17):
    r,wr=jacobi(nr,1,2);r=(r+1)/2
    u,wu=jacobi(na,1.5,1.5);v,wv=jacobi(na,.5,.5)
    uu,vv=np.meshgrid(u,v,indexing='ij');ww=(wu[:,None]*wv[None,:]).ravel()
    total=0j;minimum=np.inf;maximum=-np.inf
    for i,a in enumerate(r):
        for j,c in enumerate(r):
            e=np.zeros((2,na*na,10));e[0,:,pc[0]]=np.sqrt(a);e[0,:,pw[0]]=np.sqrt(1-a)
            e[1,:,pc[0]]=np.sqrt(c)*uu.ravel();e[1,:,pc[1]]=np.sqrt(c)*np.sqrt(1-uu.ravel()**2)
            e[1,:,pw[0]]=np.sqrt(1-c)*vv.ravel();e[1,:,pw[1]]=np.sqrt(1-c)*np.sqrt(1-vv.ravel()**2)
            mats=np.einsum('xna,xaij->nij',e,coeff);values=np.linalg.det(mats)
            total+=wr[i]*wr[j]*np.dot(ww,values)
            minimum=min(minimum,float(min(values.real)));maximum=max(maximum,float(max(values.real)))
    return total,minimum,maximum

if __name__=='__main__':
    c,pc,pw,f,d=background();print('background',d,flush=True)
    for nr,na in ((9,17),(10,18)):
        val,lo,hi=integral(c,pc,pw,nr,na)
        print('full',nr,na,f*val,'det-range',lo,hi,flush=True)
