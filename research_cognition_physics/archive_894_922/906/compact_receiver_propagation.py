"""906: compact original neutral receiver mode and its bilinear source on902.
Continuum prescription is fixed across grids; numerical results are diagnostics.
No absolute reference stress, read vertex, or full872 feedback is substituted.
"""
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.dont_write_bytecode=True
sys.path.insert(0,str(STAGE/'903'));sys.path.insert(0,str(STAGE/'905'))
import curved_fermion_propagation as shared
import canonical_source_bridge as bridge
ev=shared.ev
SIG=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
BETA=np.diag([1.,1.,-1.,-1.]);ALPHA=np.zeros((4,4,4),complex);ALPHA[0]=np.eye(4)
for i in range(3):ALPHA[i+1,:2,2:]=SIG[i];ALPHA[i+1,2:,:2]=SIG[i]
GAMMA=np.array([1j*BETA@a for a in ALPHA])
MASS=1.;CENTER=np.array([0.,np.pi/2,np.pi/4]);INNER=.75;OUTER=2.25

def mat(a,u):return np.einsum('ab,...b->...a',a,u)
def profile(r):
    r=np.asarray(r);out=np.ones_like(r,dtype=float);out[r>=OUTER]=0
    mask=(r>INNER)&(r<OUTER);t=(r[mask]-INNER)/(OUTER-INNER)
    left=np.exp(-1/t);right=np.exp(-1/(1-t));out[mask]=right/(left+right)
    return out

def normalization(order=256):
    # Radial quadrature only for the selected compact preparation, not a grid renormalization.
    x,w=np.polynomial.legendre.leggauss(order);r=INNER+(OUTER-INNER)*(x+1)/2
    integral=4*np.pi*(INNER**3/3+(OUTER-INNER)/2*np.sum(w*r*r*profile(r)**2))
    return float(integral**-.5)
NORM=normalization()

def coordinates(N):
    x=2*np.pi*np.arange(N)/N;q=np.stack(np.meshgrid(x,x,x,indexing='ij'),axis=-1)
    d=(q-CENTER+np.pi)%(2*np.pi)-np.pi
    return q,np.linalg.norm(d,axis=-1)

def initial(N):
    _,r=coordinates(N);u=np.zeros(r.shape+(4,),complex);u[...,0]=NORM*profile(r)
    return u

def connection(fr):
    out=np.zeros(fr['low'].shape[:-3]+(4,4,4),complex)
    for a in range(4):
        for b in range(a+1,4):out+=fr['low'][..., :,a,b,None,None]*(GAMMA[a]@GAMMA[b]/2)
    return out

def potential(y,z,fr,u,mass=MASS):
    out=mass*z['alpha'][...,None]*mat(BETA,u)
    for a in range(4):
        for b in range(a+1,4):
            B=GAMMA[a]@GAMMA[b]/2
            for A in range(4):
                M=-.5j*(ALPHA[A]@B-B.conj().T@ALPHA[A])
                if np.max(abs(M))>1e-14:
                    c=np.einsum('...m,...m->...',fr['W'][..., :,A],fr['low'][..., :,a,b])
                    out+=c[...,None]*mat(M,u)
    return out

def apply(model,y,z,fr,u,mass=MASS):
    out=potential(y,z,fr,u,mass)
    for i in range(3):
        du=shared.dcomplex(model,u,i);flux=np.zeros_like(u)
        for A in range(4):
            c=fr['W'][...,i+1,A,None]
            out-=.5j*c*mat(ALPHA[A],du);flux+=c*mat(ALPHA[A],u)
        out-=.5j*shared.dcomplex(model,flux,i)
    return out

def rhs(model,y,u):
    dy,z,_=model.rhs(y,aux=True);fr=shared.frame_data(y,z)
    return dy,-1j*apply(model,y,z,fr,u)

def step(model,y,u,dt):
    a,ua=rhs(model,y,u);b,ub=rhs(model,ev.add(y,a,dt/2),u+dt/2*ua)
    c,uc=rhs(model,ev.add(y,b,dt/2),u+dt/2*ub);d,ud=rhs(model,ev.add(y,c,dt),u+dt*uc)
    return {k:y[k]+dt/6*(a[k]+2*b[k]+2*c[k]+d[k]) for k in y},u+dt/6*(ua+2*ub+2*uc+ud)

def bilinear_jet(y,z,fr,u,du):
    # du has coordinate time and all spatial derivatives of half-density u.
    root=np.sqrt(z['vol']);psi=u/root[...,None]
    dlvol=.5*np.einsum('...ij,...mji->...m',z['invgamma'],z['d'][..., :,1:,1:])
    deriv=du/root[...,None,None]-.5*dlvol[..., :,None]*psi[...,None,:]
    cov=deriv+np.einsum('...mab,...b->...ma',connection(fr),psi)
    inv=fr['W']/z['alpha'][...,None,None]
    Cup=np.einsum('...mA,Aab->...mab',inv,ALPHA)
    Clow=np.einsum('...mn,...nab->...mab',y['g'],Cup)
    Z=np.einsum('...a,...mab,...nb->...mn',psi.conj(),Clow,cov)
    onshell_stress=.5*np.imag(Z+np.swapaxes(Z,-1,-2))
    scalar=np.einsum('...a,ab,...b->...',psi.conj(),BETA,psi).real
    current=np.einsum('...a,...mA,Aab,...b->...m',u.conj(),fr['W'],ALPHA,u).real
    dirac=1j*np.einsum('...mab,...mb->...a',Cup,cov)-MASS*mat(BETA,psi)
    lagrangian=-np.imag(np.einsum('...a,...mab,...mb->...',psi.conj(),Cup,cov))-MASS*scalar
    stress=onshell_stress+y['g']*lagrangian[...,None,None]
    return dict(scalar=scalar,stress=stress,onshell_stress=onshell_stress,lagrangian=lagrangian,current_density=current,dirac_residual=dirac)

def bilinears(model,y,u):
    z=model.geometry(y);fr=shared.frame_data(y,z)
    ut=-1j*apply(model,y,z,fr,u)
    du=np.stack([ut,*[shared.dcomplex(model,u,i) for i in range(3)]],axis=-2)
    return bilinear_jet(y,z,fr,u,du),z,fr,ut

def norm(u):return float((2*np.pi)**3*np.mean(np.sum(abs(u)**2,axis=-1)))

def source(model,y,u):
    data,z,fr,ut=bilinears(model,y,u);mu=z['alpha']*z['vol'];T=data['stress']
    jg=.5*mu[...,None,None]*np.einsum('...am,...mn,...nb->...ab',z['ig'],T,z['ig'])
    force,load,recovered=bridge.convert(y,z,jg,np.zeros_like(y['phi']),np.zeros(y['g'].shape[:3]+(4,12)))
    assert np.max(abs(recovered-T))<1e-12
    return data,z,force,load

def diagnostics(model,y,u,time,dt=2e-5):
    data,z,force,load=source(model,y,u);dy,ut=rhs(model,y,u)
    sides=[]
    for sign in (1,-1):
        yy=ev.add(y,dy,sign*dt);uu=u+sign*dt*ut
        dd,zz,_,_=bilinears(model,yy,uu);mu=zz['alpha']*zz['vol']
        sides.append((dd,mu[...,None,None]*np.einsum('...ma,...an->...mn',zz['ig'],dd['stress'])))
    mu=z['alpha']*z['vol'];mixed=np.einsum('...ma,...an->...mn',z['ig'],data['stress'])
    density=mu[...,None,None]*mixed
    div=(sides[0][1][...,0,:]-sides[1][1][...,0,:])/(2*dt)
    for i in range(3):div+=model.grad(density[...,i+1,:])[...,i,:]
    div-=mu[...,None]*np.einsum('...rmn,...mr->...n',z['Gamma'],mixed)
    cur=(sides[0][0]['current_density'][...,0]-sides[1][0]['current_density'][...,0])/(2*dt)
    for i in range(3):cur+=model.grad(data['current_density'][...,i+1])[...,i]
    trace=np.einsum('...mn,...mn->...',z['ig'],data['stress'])+MASS*data['scalar']
    _,r=coordinates(model.N);mask=r<=.45
    densitynorm=np.sum(abs(u)**2,axis=-1);outside=r>OUTER+2*time
    return dict(time=time,norm=norm(u),sampled_inner_scalar_min=float(data['scalar'][mask].min()),
        sampled_inner_scalar_max=float(data['scalar'][mask].max()),
        sampled_spatial_speed_bound=float(np.max(np.linalg.norm(z['beta'],axis=-1)+z['alpha']/np.sqrt(np.linalg.eigvalsh(z['gamma'])[...,0]))),
        outside_radius_test_probability=float((2*np.pi)**3*np.mean(densitynorm*outside)),
        covariant_Dirac_residual_max=float(np.max(abs(data['dirac_residual']))),
        current_density_divergence_max=float(np.max(abs(cur))),
        stress_density_divergence_max=float(np.max(abs(div))),
        on_shell_trace_error_max=float(np.max(abs(trace))),
        source_stress_max=float(np.max(abs(data['stress']))),
        offshell_completion_max=float(np.max(abs(data['stress']-data['onshell_stress']))),
        gravity_acceleration_force_max=float(np.max(abs(force['gd']))),
        normal_constraint_load_max={k:float(np.max(abs(v))) for k,v in load.items()},
        background_constraints=model.constraints(y))

def flow(N,steps=4,Tend=.01):
    assert N%2==1
    with ev.ResearchRuntime(ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=ev.Model(N,old);y,_=model.initial();u=initial(N);norm0=norm(u)
        for _ in range(steps):y,u=step(model,y,u,Tend/steps)
        row=diagnostics(model,y,u,Tend)
        row.update(N=N,steps=steps,initial_quadrature_norm=norm0,norm_drift=abs(norm(u)-norm0))
        return row,u
