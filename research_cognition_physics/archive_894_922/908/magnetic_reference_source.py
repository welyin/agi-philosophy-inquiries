"""908: off-shell first-jet derivative of the actual h/s/p/Mh material reference.
The background may solve902; variations here are independent spacetime jets.
Coordinate-density sources, not765 fiber-normalized sources, are returned.
"""
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent/'902'))
import common_background_evolution as ev

def bracket(a,b):
    z=np.zeros(np.broadcast_shapes(a.shape,b.shape),dtype=np.result_type(a,b))
    for c,i,j,f in ev.TERMS:z[...,c]+=f*a[...,i]*b[...,j]
    return z

def geometry(g,phi,dphi,A,dA):
    # dphi[mu,A], dA[mu,nu,a], A[mu,a]; independent off-shell jets.
    G=np.linalg.inv(g);h=np.sqrt(np.sum(phi[...,:4]**2,axis=-1));radial=phi[...,:4]/h[...,None]
    n=np.einsum('...A,...mA->...m',radial,dphi[...,:,:4]);v=np.einsum('...mn,...n->...m',G,n)
    q=np.sum(n*v,axis=-1);P=G-v[..., :,None]*v[...,None,:]/q[...,None,None]
    F=dA-np.swapaxes(dA,-3,-2)+bracket(A[..., :,None,:],A[...,None,:,:])
    Fc=F[...,:8];B=np.einsum('...mr,...ns,...rsa->...mna',P,P,Fc)
    M=.5*np.sum(B*Fc,axis=(-3,-2,-1))
    H=np.einsum('...ns,...mna,...rsa->...mr',P,Fc,Fc)
    metric=-np.einsum('...am,...mn,...nb->...ab',P,H,P)
    clock=-2*np.einsum('...am,...mn,...n->...a',P,H,v)/q[...,None]
    X=np.stack((h,phi[...,4],phi[...,5],M),axis=-1)
    return dict(X=X,h=h,radial=radial,n=n,q=q,P=P,F=F,B=B,metric=metric,clock=clock)

def variation(data,phi,dphi,A,dg,dphi0,ddphi,dA0,ddA):
    dh=np.sum(data['radial']*dphi0[...,:4],axis=-1)
    dr=(dphi0[...,:4]-data['radial']*dh[...,None])/data['h'][...,None]
    dn=np.einsum('...A,...mA->...m',dr,dphi[...,:,:4])+np.einsum('...A,...mA->...m',data['radial'],ddphi[...,:,:4])
    dF=ddA-np.swapaxes(ddA,-3,-2)+bracket(dA0[..., :,None,:],A[...,None,:,:])+bracket(A[..., :,None,:],dA0[...,None,:,:])
    dM=np.sum(data['metric']*dg,axis=(-2,-1))+np.sum(data['clock']*dn,axis=-1)+np.sum(data['B']*dF[...,:8],axis=(-3,-2,-1))
    return np.stack((dh,dphi0[...,4],dphi0[...,5],dM),axis=-1)

def first_jet_source(data,phi,dphi,A,current,k):
    # Pairing: current^mu_a delta A_mu^a + k_b delta X^b.
    # This forms all derivative potentials before adjoint spacetime divergence.
    km=k[...,3];metric=km[...,None,None]*data['metric'];clock=km[...,None]*data['clock']
    Pphi=np.zeros(phi.shape,dtype=np.result_type(phi,current,k));Dphi=np.zeros(dphi.shape,dtype=Pphi.dtype)
    tangential=dphi[...,:,:4]-data['n'][..., :,None]*data['radial'][...,None,:]
    Pphi[...,:4]=k[...,0,None]*data['radial']+np.einsum('...m,...mA->...A',clock,tangential)/data['h'][...,None]
    Pphi[...,4]=k[...,1];Pphi[...,5]=k[...,2]
    Dphi[...,:,:4]=clock[..., :,None]*data['radial'][...,None,:]
    PA=np.array(current,copy=True);DA=np.zeros(A.shape[:-2]+(4,4,12),dtype=PA.dtype)
    DA[...,:8]=2*km[...,None,None,None]*data['B']
    # 2k B^{mu nu} i[A_mu,delta A_nu] as a coordinate partial coefficient.
    for c,i,j,f in ev.TERMS:
        PA[..., :,j]+=f*np.einsum('...mn,...m->...n',DA[..., :, :,c],A[..., :,i])
    return dict(g=metric,phi=Pphi,dphi=Dphi,A=PA,dA=DA)

def pair(source,v):
    return (
        np.sum(source['g']*v['g'],axis=(-2,-1))+np.sum(source['phi']*v['phi'],axis=-1)
        +np.sum(source['dphi']*v['dphi'],axis=(-2,-1))+np.sum(source['A']*v['A'],axis=(-2,-1))
        +np.sum(source['dA']*v['dA'],axis=(-3,-2,-1)))
