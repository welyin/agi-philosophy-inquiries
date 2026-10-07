"""905: Euler source densities -> harmonic/temporal canonical forcing and constraints.
Verification uses a declared covariant gauge-kinetic variation f(p)F^2;
it is not asserted to be the quantum source872 or a new model principle.
"""
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'902'))
import common_background_evolution as ev

def convert(y,z,jg,jphi,jA):
    mu=z['alpha']*z['vol'];s=2*np.einsum('...ma,...ab,...bn->...mn',y['g'],jg,y['g'])/mu[...,None,None]
    trace=np.einsum('...mn,...mn->...',z['ig'],s)
    out={k:np.zeros_like(v) for k,v in y.items()}
    out['gd']=(-2*s+y['g']*trace[...,None,None])/z['ig'][...,0,0,None,None]
    out['pi']=jphi;out['E']=jA[...,1:,:]
    normal=np.concatenate((np.ones(mu.shape+(1,)),-z['beta']),axis=-1)/z['alpha'][...,None]
    nn=np.einsum('...m,...mn,...n->...',normal,s,normal);ni=np.einsum('...m,...mi->...i',normal,s[..., :,1:])
    return out,dict(H=-2*nn,M=ni,Gauss=jA[...,0,:]),s

def ym_parts(model,y,z,m):
    K=model.K;iv=z['invgamma'];gam=z['gamma'];e=m['e'];curv=m['curv'];up=m['up']
    e2=np.einsum('...ij,...ia,...ja,a->...',iv,e,e,K);f2=np.einsum('...ija,...ija,a->...',curv,up,K)
    rho=.5*e2+.25*f2
    stress=-np.einsum('...ia,...ja,a->...ij',e,e,K)+np.einsum('...kl,...ika,...jla,a->...ij',iv,curv,curv,K)+gam*(.5*e2-.25*f2)[...,None,None]
    mom=np.einsum('...ja,...ija->...i',y['E'],curv)
    alpha=z['alpha'];beta=z['beta'];vol=z['vol']
    T=np.zeros_like(y['g']);T[...,1:,1:]=stress
    T[...,0,1:]=T[...,1:,0]=alpha[...,None]*mom/vol[...,None]+np.einsum('...j,...ji->...i',beta,stress)
    T[...,0,0]=alpha**2*rho+2*alpha*np.einsum('...i,...i->...',beta,mom)/vol+np.einsum('...i,...j,...ij->...',beta,beta,stress)
    flux=(alpha*vol)[...,None,None,None]*K*up
    flux+=beta[..., :,None,None]*y['E'][...,None,:,:]-beta[...,None,:,None]*y['E'][..., :,None,:]
    return dict(T=T,rho=rho,mom=mom,curvature_square=f2-2*e2,flux=flux)

def div_flux(model,A,flux):
    return np.einsum('...jjia->...ia',model.grad(flux))+sum(ev.bracket(A[...,j,None,:],flux[...,j,:,:]) for j in range(3))

def gauge_kinetic_source(model,y):
    rhs,z,m=model.rhs(y,aux=True);ym=ym_parts(model,y,z,m);profile=y['phi'][...,5];pdot=rhs['phi'][...,5]
    mu=z['alpha']*z['vol'];s=profile[...,None,None]*ym['T']
    jg=.5*mu[...,None,None]*np.einsum('...am,...mn,...nb->...ab',z['ig'],s,z['ig'])
    jphi=np.zeros_like(y['phi']);jphi[...,5]=-.25*mu*ym['curvature_square']
    jA=np.zeros(y['A'].shape[:3]+(4,12))
    fE=profile[...,None,None]*y['E']
    jA[...,0,:]=np.einsum('...iia->...a',model.grad(fE))+sum(ev.bracket(y['A'][...,i,:],fE[...,i,:]) for i in range(3))
    jA[...,1:,:]=div_flux(model,y['A'],profile[...,None,None,None]*ym['flux'])-profile[...,None,None]*rhs['E']-pdot[...,None,None]*y['E']
    return jg,jphi,jA,ym

def deformed_velocity_in_base_variables(model,y,lam):
    """Independent equations of S_base-lambda/4 int mu p K F^2.
    Canonical E_total=(1+lambda*p) E_base. Returned derivative is of E_base.
    Full deformed physical source and new scalar force are retained.
    """
    rhs,z,m=model.rhs(y,aux=True);ym=ym_parts(model,y,z,m)
    profile=y['phi'][...,5];Z=1+lam*profile;assert Z.min()>0
    # Metric force from the actual changed stress, not from source conversion.
    T=m['T']+lam*profile[...,None,None]*ym['T'];trace=np.einsum('...mn,...mn->...',z['ig'],T)
    dgd=model.grad(y['gd']);wave=2*z['Q']-2*T+y['g']*trace[...,None,None]
    wave-=2*np.einsum('...i,...imn->...mn',z['ig'][...,0,1:],dgd)
    wave-=model.lapweighted(y['g'],z['ig'][...,1:,1:])
    out={k:v.copy() for k,v in rhs.items()};out['gd']=wave/z['ig'][...,0,0,None,None]
    out['pi'][...,5]-=.25*lam*z['alpha']*z['vol']*ym['curvature_square']
    # YM flux uses E_total and K_total; the same scalar current is subtracted once.
    oldYM=div_flux(model,y['A'],ym['flux']);matter_current=oldYM-rhs['E']
    total_E_dot=div_flux(model,y['A'],Z[...,None,None,None]*ym['flux'])-matter_current
    out['E']=(total_E_dot-lam*rhs['phi'][...,5,None,None]*y['E'])/Z[...,None,None]
    return out

def deformed_constraints(model,y,lam):
    z=model.geometry(y);m=model.matter(y,z);ym=ym_parts(model,y,z,m);p=y['phi'][...,5]
    Etotal=(1+lam*p)[...,None,None]*y['E']
    gauss=np.einsum('...iia->...a',model.grad(Etotal))+sum(ev.bracket(y['A'][...,i,:],Etotal[...,i,:]) for i in range(3))+np.einsum('...A,...aA->...a',y['pi'],m['rp'])
    return dict(H_extra=-2*lam*p*ym['rho'],M_extra=lam*p[...,None]*ym['mom']/z['vol'][...,None],Gauss=gauss)
