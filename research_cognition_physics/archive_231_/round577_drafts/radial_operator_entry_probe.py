"""577 entry diagnostic only: adjacent-node transfer without an added pointer.
Radial four-coordinate operator discretization tests the algebra; not full Gauss evolution.
"""
import json
from pathlib import Path
import numpy as np
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_pointer_matter_compatibility as spectrum
old=spectrum.old
TARGET=HERE/'radial_operator_entry_probe_results.json'


def ratio(z):
    small=abs(z)<1e-5
    safe=np.where(small,1e-4,z)
    r=np.arccosh(1+safe)/np.sqrt(safe*(2+safe))
    dr=(1-(1+safe)*r)/(safe*(2+safe))
    r=np.where(small,1-z/3+2*z*z/15-2*z**3/35,r)
    dr=np.where(small,-1/3+4*z/15-6*z*z/35,dr)
    return r,dr


def arrays(n):
    width=.32;dx=2*width/n
    axes=[c+(np.arange(n)-n//2)*dx for c in np.tile(spectrum.VAC,2)]
    coords=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1)
    x,y=coords[...,:2],coords[...,2:]
    fx=2-np.sum(x*x,axis=-1)/6;fy=2-np.sum(y*y,axis=-1)/6
    def inverse(q,f):return f[...,None,None]*(np.eye(2)-q[..., :,None]*q[...,None,:]/12)
    kx,ky=inverse(x,fx),inverse(y,fy)
    mux=x[...,0]**3*fx**-3;muy=y[...,0]**3*fy**-3
    mu=mux*muy
    ax,ay=np.sqrt(fx),np.sqrt(fy)
    z=(np.sum((x-y)**2,axis=-1)/12+(ax-ay)**2/2)/(ax*ay)
    edge=12*np.arcsinh(np.sqrt(np.maximum(z,0)/2))**2
    L=spectrum.MAT;u=spectrum.U0
    d=x*x-u;e=y*y-u
    onsite_x=np.einsum('...i,ij,...j->...',d,L,d)/(4*fx*fx)
    onsite_y=np.einsum('...i,ij,...j->...',e,L,e)/(4*fy*fy)
    V=onsite_x+onsite_y+edge
    r,dr=ratio(z);C=1+z;den=ax*ay
    CA=-y/(6*den[...,None])+C[...,None]*x/(6*fx[...,None])
    CB=-x/(6*den[...,None])+C[...,None]*y/(6*fy[...,None])
    CAB=-np.eye(2)/(6*den[...,None,None])-y[..., :,None]*y[...,None,:]/(36*(fy*den)[...,None,None])
    CAB+=x[..., :,None]*CB[...,None,:]/(6*fx[...,None,None])
    hess=6*(dr[...,None,None]*CA[..., :,None]*CB[...,None,:]+r[...,None,None]*CAB)
    dT=np.stack((x[...,0],np.zeros_like(fx)),axis=-1)
    db=np.stack((np.zeros_like(fy),np.cos(y[...,1])),axis=-1)
    va=np.einsum('...ij,...j->...i',kx,dT);vb=np.einsum('...ij,...j->...i',ky,db)
    coefficient=-np.einsum('...i,...ij,...j->...',va,hess,vb)
    offsets=coords-np.tile(spectrum.VAC,2)
    sigma=.075;radius=.24
    t=np.sum((offsets/radius)**2,axis=-1);cut=np.maximum(1-t,1e-300)
    logamp=-np.sum(offsets*offsets,axis=-1)/(4*sigma*sigma)-.3/cut
    amp=np.exp(np.maximum(logamp,-745));amp[t>=1]=0
    amp/=np.sqrt(np.sum(mu*amp*amp)*dx**4)
    return dict(dx=dx,coords=coords,kx=kx,ky=ky,mux=mux,muy=muy,mu=mu,V=V,
                amp=amp,T=x[...,0]**2/2,b=np.sin(y[...,1]),coefficient=coefficient,
                preparation_energy_density=np.einsum('...i,...ij,...j->...',dT,kx,dT))


def operator_probe(n,hbar=.07,theta=.2):
    q=arrays(n);dx=q['dx'];mu=q['mu']
    def D(f,axis):return (np.roll(f,-1,axis=axis)-np.roll(f,1,axis=axis))/(2*dx)
    def H(psi):
        out=q['V']*psi
        for start,mat,measure in ((0,q['kx'],q['mux']),(2,q['ky'],q['muy'])):
            gradients=[D(psi,start+a) for a in range(2)]
            div=sum(D(measure*sum(mat[...,i,j]*gradients[j] for j in range(2)),start+i) for i in range(2))
            out=out-hbar*hbar*div/(2*measure)
        return out
    def inner(a,b):return np.sum(mu*a.conj()*b)*dx**4
    psi=q['amp']*np.exp(1j*theta*q['T']/hbar)
    h1=H(psi);h2=H(h1);h3=H(h2);b=q['b']
    d1=-2*inner(h1,b*psi).imag/hbar
    d2=-2*(inner(h2,b*psi).real-inner(h1,b*h1).real)/(hbar*hbar)
    d3=(2*inner(h3,b*psi).imag-6*inner(h2,b*h1).imag)/(hbar**3)
    coefficient=float(np.sum(mu*q['amp']**2*q['coefficient'])*dx**4)
    h0=H(q['amp']);energy=float(inner(psi,h1).real-inner(q['amp'],h0).real)
    continuum_cost=theta*theta*np.sum(mu*q['amp']**2*q['preparation_energy_density'])*dx**4/2
    symmetry=abs(inner(psi,h1)-inner(h1,psi))
    return dict(n=n,hbar=hbar,theta=theta,first_derivative=float(d1),second_derivative=float(d2),
        third_derivative=float(d3),predicted_third_derivative=theta*coefficient,
        relative_third_derivative_error=abs(d3-theta*coefficient)/abs(theta*coefficient),
        phase_preparation_energy=energy,continuum_preparation_energy=float(continuum_cost),
        Hermiticity_error=float(symmetry))


def run():
    rows=[operator_probe(n) for n in (13,19,27)]
    assert rows[-1]['relative_third_derivative_error']<rows[0]['relative_third_derivative_error']
    assert all(r['third_derivative']<0 and abs(r['first_derivative'])<1e-10 for r in rows)
    return dict(entry_only=True,radial_diagnostic_not_full_Gauss_simulation=True,rows=rows)


if __name__=='__main__':
    result=run()
    if '--write-results' in sys.argv:
        with TARGET.open('x',encoding='utf8') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
