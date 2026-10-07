"""903: original 64-component Dirac/Majorana propagation on the902 joint background.
Full spatial PDE; finite collocation diagnostics, not a quantum-state/source certificate.
"""
from pathlib import Path
import sys,json,argparse,time
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.dont_write_bytecode=True
sys.path.insert(0,str(STAGE/'902'));sys.path.insert(0,str(STAGE/'802'))
import common_background_evolution as ev
import original_bff_vertex as vertex
G=np.array(vertex.original.GAMMA);ALPHA=np.array([np.eye(64),*G]);ETA=vertex.ETA
MASS=np.array(vertex.NUMERATORS)
CHARGE=np.r_[np.arange(32,64),np.arange(32)]

def internal_generators():
    q=np.zeros((12,32,32),complex);sig=vertex.original.old.chiral.SIG
    for name,sl in vertex.original.old.matter.SLICES.items():
        n=sl.stop-sl.start
        for a in range(8):
            if name in ('Q','u','d'):q[a,sl,sl]=np.kron(ev.T[a],np.eye(n//3))
        for a in range(3):
            if name in ('Q','L'):q[8+a,sl,sl]=np.kron(np.eye(n//4),np.kron(sig[a]/2,np.eye(2)))
        q[11,sl,sl]=vertex.original.old.chiral.CHARGES[name]*np.eye(n)
    out=np.zeros((12,64,64),complex);out[:,:32,:32]=q;out[:,32:,32:]=-q.conj()
    return out
Q=internal_generators()
SPIN=[]
for a in range(4):
    for b in range(a+1,4):
        # sum over both a,b of omega_ab gamma^a gamma^b/4
        B=vertex.GAMMA[a]@vertex.GAMMA[b]/2
        for A in range(4):
            mat=-.5j*(ALPHA[A]@B-B.conj().T@ALPHA[A])
            if np.max(abs(mat))>1e-13:SPIN.append((A,a,b,mat))

def mm(mat,u):
    # Constant matrices have at most two nonzeros in each row.
    out=np.zeros_like(u)
    for row in range(64):
        cols=np.flatnonzero(abs(mat[row])>1e-14)
        if len(cols):out[...,row,:]=np.einsum('j,...jr->...r',mat[row,cols],u[...,cols,:])
    return out

def dcomplex(model,u,axis):
    h=np.fft.fftn(u,axes=(0,1,2));k=model.k[...,axis].reshape(model.k.shape[:3]+(1,)*(u.ndim-3))
    return np.fft.ifftn(1j*k*h,axes=(0,1,2))

def frame_data(y,z,rotation=None):
    """Time-gauge tetrad with symmetric positive spatial square root.
    Analytic matrix-root derivative handles repeated eigenvalues by sums, not differences.
    rotation is optional (Lambda, dLambda), only for independent gauge-covariance checks.
    """
    gamma=z['gamma'];vals,O=np.linalg.eigh(gamma);roots=np.sqrt(vals)
    B=np.einsum('...ai,...i,...bi->...ab',O,roots,O)
    dg=z['d'];dd=dg[...,1:,1:]
    local=np.einsum('...ai,...mab,...bj->...mij',O,dd,O)
    local/=roots[...,None,:,None]+roots[...,None,None,:]
    dB=np.einsum('...ai,...mij,...bj->...mab',O,local,O)
    ig=z['ig'];di=-np.einsum('...ab,...mbc,...cd->...mad',ig,dg,ig)
    alpha=z['alpha'];beta=z['beta']
    da=.5*alpha[...,None]**3*di[...,0,0]
    db=-di[...,0,1:]/ig[...,0,0,None,None]+ig[...,0,1:][...,None,:]*di[...,0,0,None]/ig[...,0,0,None,None]**2
    E=np.zeros(y['g'].shape);E[...,0,0]=alpha;E[...,1:,1:]=B;E[...,1:,0]=np.einsum('...ai,...i->...a',B,beta)
    dE=np.zeros(dg.shape);dE[...,0,0]=da;dE[...,1:,1:]=dB
    dE[...,1:,0]=np.einsum('...mai,...i->...ma',dB,beta)+np.einsum('...ai,...mi->...ma',B,db)
    if rotation is not None:
        L,dL=rotation
        dE=np.einsum('...mab,...bi->...mai',dL,E)+np.einsum('...ab,...mbi->...mai',L,dE)
        E=np.einsum('...ab,...bi->...ai',L,E)
    inv=np.linalg.inv(E)
    omega=np.einsum('...ar,...rms,...sb->...mab',E,z['Gamma'],inv)-np.einsum('...mar,...rb->...mab',dE,inv)
    low=np.einsum('ab,...mbc->...mac',ETA,omega)
    W=alpha[...,None,None]*inv
    dinv=-np.einsum('...ia,...mab,...bj->...mij',inv,dE,inv)
    dW=da[..., :,None,None]*inv[...,None,:,:]+alpha[...,None,None,None]*dinv
    return dict(E=E,dE=dE,W=W,dW=dW,low=low)

def terms(y,z,f):
    alpha=z['alpha'];phi=y['phi'][...,:5];F=2-np.sum(phi*phi,axis=-1)/6
    for a in range(5):yield alpha*phi[...,a]/np.sqrt(F),MASS[a]
    for a in range(12):
        for A in range(4):
            coeff=np.einsum('...i,...i->...',f['W'][...,1:,A],y['A'][..., :,a])
            yield coeff,ALPHA[A]@Q[a]
    for A,a,b,mat in SPIN:yield np.einsum('...m,...m->...',f['W'][..., :,A],f['low'][..., :,a,b]),mat

def potential(y,z,f,u,omit_spin=False,omit_color=False):
    out=np.zeros_like(u)
    for num,(c,mat) in enumerate(terms(y,z,f)):
        if omit_spin and num>=53:continue
        if omit_color and 5<=num<37:continue
        if np.max(abs(c))>0:out+=c[...,None,None]*mm(mat,u)
    return out

def apply(model,y,z,f,u,omit_spin=False,omit_color=False):
    # Split derivative is discretely skew-adjoint even when collocation lacks Leibniz.
    out=potential(y,z,f,u,omit_spin,omit_color)
    for i in range(3):
        du=dcomplex(model,u,i);flux=np.zeros_like(u)
        for a in range(4):
            c=f['W'][...,i+1,a,None,None]
            out+=-.5j*c*mm(ALPHA[a],du)
            flux+=c*mm(ALPHA[a],u)
        out+=-.5j*dcomplex(model,flux,i)
    return out

def point_operator(y,z,f,u,du):
    result=potential(y,z,f,u)
    for i in range(3):
        for a in range(4):
            result-=1j*f['W'][...,i+1,a,None,None]*mm(ALPHA[a],du[i])
            result-=.5j*f['dW'][...,i+1,i+1,a,None,None]*mm(ALPHA[a],u)
    return result

def initial_modes(model):
    # Smooth propagator test vectors, not a replacement for730's prepared covariance.
    n=model.N;x=2*np.pi*np.arange(n)/n;X,Y,Z=np.meshgrid(x,x,x,indexing='ij')
    u=np.zeros((n,n,n,64,2),complex)
    u[...,30,0]=(1+.2*np.cos(X)+.1*np.sin(Y+Z))*np.exp(1j*X)
    u[...,0,1]=(1+.15*np.cos(Z)+.08*np.sin(X-Y))*np.exp(1j*(Y-Z))
    for r in range(2):u[...,r]/=np.sqrt(np.mean(np.sum(abs(u[...,r])**2,axis=-1)))
    return u

def rhs(model,y,u):
    dy,z,_=model.rhs(y,aux=True);f=frame_data(y,z)
    return dy,-1j*apply(model,y,z,f,u)

def step(model,y,u,dt):
    a,ua=rhs(model,y,u)
    b,ub=rhs(model,ev.add(y,a,dt/2),u+dt/2*ua)
    c,uc=rhs(model,ev.add(y,b,dt/2),u+dt/2*ub)
    d,ud=rhs(model,ev.add(y,c,dt),u+dt*uc)
    return {k:y[k]+dt/6*(a[k]+2*b[k]+2*c[k]+d[k]) for k in y},u+dt/6*(ua+2*ub+2*uc+ud)

def gram(u,v=None):
    if v is None:v=u
    return u.reshape(-1,u.shape[-1]).conj().T@v.reshape(-1,v.shape[-1])/np.prod(u.shape[:3])

def flow(N,steps,Tend=.01):
    assert N%2==1,'Odd Fourier grids preserve conjugation without a Nyquist convention.'
    with ev.ResearchRuntime(ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=ev.Model(N,old);y,_=model.initial();u=initial_modes(model);first=u.copy();initial_gram=gram(u)
        for _ in range(steps):y,u=step(model,y,u,Tend/steps)
        coeff=np.fft.fftn(u,axes=(0,1,2))/N**3
        targets={','.join(map(str,k)):[[float(z.real),float(z.imag)] for z in coeff[tuple(j%N for j in k)][[30,31,0,12,62,32],[0,0,1,1,0,1]]] for k in ((1,0,0),(0,1,-1),(1,1,0))}
        row=dict(N=N,steps=steps,Tend=Tend,gram_error=float(np.max(abs(gram(u)-initial_gram))),
            mode_change_norms=np.sqrt(np.mean(np.sum(abs(u-first)**2,axis=-2),axis=(0,1,2))).tolist(),
            selected_Fourier_coefficients=targets,background_constraints=model.constraints(y))
        return row,model,y,u
