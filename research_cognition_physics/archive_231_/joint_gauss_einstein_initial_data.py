"""572: same gauge source, a compact CMC obstruction, and a changed-source completion.

Four-dimensional Einstein-frame action and canonical matching are inputs.
Numerical elliptic solves illustrate an analytic existence argument, not GR emergence.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_gauss_continuum_sampling as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauss_einstein_initial_data_results.json'
M02=2.
VOL=old.LBOX**3


def fft(x):return np.fft.fftn(x,axes=(0,1,2))
def ifft(x):return np.fft.ifftn(x,axes=(0,1,2)).real


def waves(N):
    k=np.fft.fftfreq(N,d=old.LBOX/N)*2*np.pi
    return np.stack(np.meshgrid(k,k,k,indexing='ij'),axis=-1)


def derivative(x,mu):
    N=x.shape[0];k=waves(N)[...,mu]
    while k.ndim<x.ndim:k=k[...,None]
    return ifft(1j*k*fft(x))


def laplace(x):return ifft(-np.sum(waves(x.shape[0])**2,axis=-1)*fft(x))


def make_source(N,repair=True):
    dx=old.LBOX/N
    grid=np.stack(np.meshgrid(*([np.arange(N)*dx]*3),indexing='ij'),axis=-1)
    f=old.fields(grid);x,y,z=np.moveaxis(grid,-1,0)
    hstar=np.sqrt(old.PAR['h2']);mat,u,_=old.scalar.parameters();sstar=np.sqrt(u[1])
    dh=np.zeros(grid.shape);ds=np.zeros_like(dh)
    dh[...,0]=dh[...,1]=.05*hstar*np.cos(x+y);ds[...,1]=-.06*sstar*np.sin(y)
    da=np.zeros(grid.shape[:-1]+(3,3,3));da0=np.zeros(grid.shape[:-1]+(3,3))
    da[...,0,0,0]=-old.A0*old.U0*np.sin(x);da[...,1,1,1]=.2*old.B0*np.cos(y)
    da[...,2,2,2]=-.12*np.sin(z)
    da0[...,1,0]=-.03*np.sin(y);da0[...,0,1]=.08*np.cos(x)
    PX=f['PX'].copy();Ps=f['Ps'].copy()
    H=dh.reshape(-1,3);S=ds.reshape(-1,3);gram=dx**3*(H.T@H+S.T@S)
    expected=VOL*np.array([.001*hstar+.0028,.001*hstar,0.])
    coeff=np.linalg.pinv(gram,rcond=1e-12)@expected
    ph=-np.einsum('...i,i->...',dh,coeff);ps=-np.einsum('...i,i->...',ds,coeff)
    if repair:PX[...,1]+=ph;Ps+=ps
    DX=np.zeros(grid.shape[:-1]+(3,2),complex);Fw=np.zeros(grid.shape[:-1]+(3,3,3))
    F0=np.zeros(grid.shape[:-1]+(3,3));mom=np.zeros(grid.shape)
    canonical=np.zeros_like(mom)
    for i in range(3):
        DX[...,i,1]=dh[...,i]
        DX[...,i,:]+=1j*np.einsum('...a,aij,...j->...i',f['a'][...,i,:],old.T,f['X'])+3j*f['a0'][...,i,None]*f['X']
        mom[...,i]=np.sum(PX.conj()*DX[...,i,:],axis=-1).real+Ps*ds[...,i]
        canonical[...,i]=PX[...,1].real*dh[...,i]+Ps*ds[...,i]
        for j in range(3):
            Fw[...,i,j,:]=da[...,i,j,:]-da[...,j,i,:]-np.cross(f['a'][...,i,:],f['a'][...,j,:])
            F0[...,i,j]=da0[...,i,j]-da0[...,j,i]
            mom[...,i]+=np.sum(f['E'][...,j,:]*Fw[...,i,j,:],axis=-1)+f['E0'][...,j]*F0[...,i,j]
            canonical[...,i]+=np.sum(f['E'][...,j,:]*da[...,i,j,:],axis=-1)+f['E0'][...,j]*da0[...,i,j]
    phi=np.concatenate((f['X'].real,f['X'].imag,f['s'][...,None]),axis=-1)
    p=np.concatenate((PX.real,PX.imag,Ps[...,None]),axis=-1)
    Dphi=np.concatenate((DX.real,DX.imag,ds[...,None]),axis=-1)
    F=M02-np.sum(phi*phi,axis=-1)/6
    # K^{-1}=F(I-phi phi^T/(6 M0²)); this acts on all five real components.
    pKp=F*(np.sum(p*p,axis=-1)-np.sum(p*phi,axis=-1)**2/(6*M02))
    B=np.sum(Dphi*Dphi,axis=(-1,-2))/F+np.sum(np.sum(Dphi*phi[...,None,:],axis=-1)**2,axis=-1)/(6*F**2)
    delta=np.stack((f['h']**2-u[0],f['s']**2-u[1]),axis=-1)
    V=np.einsum('...i,ij,...j->...',delta,mat,delta)/4
    U=V/F**2
    bw,b0=old.PAR['b'][1:];kw,k0=old.PAR['K'][1:]
    Y=bw*np.sum(f['E']**2,axis=(-1,-2))+b0*np.sum(f['E0']**2,axis=-1)
    for i in range(3):
        for j in range(i+1,3):Y+=kw*np.sum(Fw[...,i,j,:]**2,axis=-1)/2+k0*F0[...,i,j]**2/2
    fmin=M02-((1.05*hstar)**2+(1.06*sstar)**2)/6
    delta_bound=np.array([u[0]*(1.05**2-1),u[1]*(1.06**2-1)])
    ubound=float(np.linalg.eigvalsh(mat)[-1]*(delta_bound@delta_bound)/(4*fmin**2))
    tau2=3*ubound+3
    ymin=kw*(old.A0*old.B0*(1-old.U0)*.8)**2/2
    return dict(N=N,dx=dx,grid=grid,f=f,PX=PX,Ps=Ps,phi=phi,p=p,Dphi=Dphi,F=F,
                mom=mom,canonical=canonical,expected=expected,gram=gram,counterflow=coeff,
                counterflow_norm=float(np.sqrt(dx**3*np.sum(ph*ph+ps*ps))),
                B=B,pKp=pKp,U=U,Y=Y,tau2=tau2,C=2*tau2/3-2*U,
                F_lower=fmin,U_upper=ubound,Y_lower=ymin)


def solve_momentum(q):
    M=q['mom'];k=waves(q['N']);k2=np.sum(k*k,axis=-1);safe=k2.copy();safe[0,0,0]=1
    mh=fft(M);dot=np.sum(k*mh,axis=-1)
    wh=(mh-k*dot[...,None]/(4*safe[...,None]))/safe[...,None]
    wh[0,0,0,:]=0
    divw=1j*np.sum(k*wh,axis=-1)
    ah=1j*(k[..., :,None]*wh[...,None,:]+k[...,None,:]*wh[..., :,None])
    ah-=2/3*divw[...,None,None]*np.eye(3)
    A=ifft(ah)
    divA=ifft(1j*np.einsum('...j,...ij->...i',k,ah))
    return A,dict(momentum_residual=float(np.max(abs(divA+M))),
                  trace_residual=float(np.max(abs(np.trace(A,axis1=-2,axis2=-1)))))


def cg(operator,b,precondition,tol=1e-12):
    x=np.zeros_like(b);r=b.copy();z=precondition(r);p=z.copy();rz=float(np.sum(r*z))
    bnorm=float(np.linalg.norm(b))
    if bnorm==0:return x,0
    for it in range(180):
        ap=operator(p);den=float(np.sum(p*ap));assert den>0
        alpha=rz/den;x+=alpha*p;r-=alpha*ap
        if np.linalg.norm(r)<tol*bnorm:return x,it+1
        z=precondition(r);new=float(np.sum(r*z));p=z+(new/rz)*p;rz=new
    raise AssertionError('CG did not converge')


def solve_hamiltonian(q,initial=1.):
    tensor,info=solve_momentum(q);A=np.sum(tensor*tensor,axis=(-1,-2))+q['pKp']
    B,Y,C=q['B'],q['Y'],q['C'];psi=np.full(B.shape,initial)
    def reaction(z):return C*z**5-B*z-A*z**-7-2*Y*z**-3
    def residual(z):return -8*laplace(z)+reaction(z)
    history=[];k2=np.sum(waves(q['N'])**2,axis=-1)
    for iteration in range(24):
        rr=residual(psi);error=float(np.max(abs(rr)));history.append(error)
        if error<2e-10:break
        J=5*C*psi**4-B+7*A*psi**-8+6*Y*psi**-4
        assert np.min(J)>0
        op=lambda v:-8*laplace(v)+J*v
        pre=lambda v:ifft(fft(v)/(8*k2+np.mean(J)))
        step,_=cg(op,-rr,pre)
        fraction=1.
        while fraction>2**-18:
            trial=psi+fraction*step
            if np.min(trial)>0 and np.max(abs(residual(trial)))<error:break
            fraction*=.5
        assert fraction>2**-18
        psi=trial
    else:raise AssertionError('Newton did not converge')
    # Direct original Hamiltonian reconstruction, not just the rearranged equation.
    curvature=-8*psi**-5*laplace(psi)
    shear_norm=psi**-12*np.sum(tensor*tensor,axis=(-1,-2))
    rho=.5*psi**-12*q['pKp']+.5*psi**-4*B+q['U']+psi**-8*Y
    einstein=curvature-shear_norm+2*q['tau2']/3-2*rho
    lower=min(.5,(q['Y_lower']/(2*(2*q['tau2']/3)))**(1/8))
    # This upper barrier uses collocation extrema only, unlike the analytic general theorem.
    upper=(1+(float(np.max(A))+float(np.max(B))+2*float(np.max(Y)))/float(np.min(C)))**.25
    assert np.max(reaction(np.full_like(psi,lower)))<0
    assert np.min(reaction(np.full_like(psi,upper)))>0
    assert lower<np.min(psi)<=np.max(psi)<upper
    return psi,tensor,dict(**info,Newton_iterations=len(history)-1,
        elliptic_residual=history[-1],original_Hamiltonian_residual=float(np.max(abs(einstein))),
        psi_min=float(np.min(psi)),psi_max=float(np.max(psi)),
        lower_barrier=lower,grid_upper_barrier=upper,min_curvature=float(np.min(curvature)),
        max_curvature=float(np.max(curvature)),tau=float(np.sqrt(q['tau2'])),
        density_min=float(np.min(rho)),density_max=float(np.max(rho)))


def inherited_probe_check():
    path=HERE/'round572_drafts/momentum_entry_probe.py'
    spec=importlib.util.spec_from_file_location('probe572',path);probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
    got=probe.run();saved=json.loads((path.with_name('momentum_entry_probe_results.json')).read_text('utf8'))
    assert got==saved
    return dict(saved_probe_reproduced=True,covariant_vs_canonical_error=got['Gauss_integration_by_parts_error'])


def original_source_obstruction_check():
    q=make_source(12,False);p=q['dx']**3*np.sum(q['mom'],axis=(0,1,2))
    pc=q['dx']**3*np.sum(q['canonical'],axis=(0,1,2))
    assert np.max(abs(p-q['expected']))<1e-12 and np.max(abs(p-pc))<1e-12
    return dict(total_momentum=p.tolist(),expected=q['expected'].tolist(),
                impossible_for_same_densities_flat_conformal_periodic_CMC=True)


def counterflow_and_Gauss_check():
    q=make_source(12);raw=make_source(12,False);qw,q0=old.charges(q['f']['X'],q['PX'])
    gw,g0=old.charges(q['f']['X'],raw['PX'])
    change=max(np.max(abs(qw-gw)),np.max(abs(q0-g0)))
    mean=q['dx']**3*np.sum(q['mom'],axis=(0,1,2))
    assert change<1e-14 and np.max(abs(mean))<1e-12
    residual=np.array([0.,0.,1.])-q['gram']@np.linalg.pinv(q['gram'])@np.array([0.,0.,1.])
    assert np.linalg.norm(residual)==1
    return dict(counterflow=q['counterflow'].tolist(),flat_norm=q['counterflow_norm'],
        total_momentum_after=mean.tolist(),Gauss_change=float(change),Gram_rank=int(np.linalg.matrix_rank(q['gram'])),
        z_target_outside_range=True)


def frame_and_positive_coefficients_check():
    q=make_source(12);phi=q['phi'][2,3,4];F=q['F'][2,3,4]
    df=-phi/3;K=np.eye(5)/F+1.5*np.outer(df,df)/F**2
    inverse=F*(np.eye(5)-np.outer(phi,phi)/(6*M02))
    err=float(np.max(abs(K@inverse-np.eye(5))))
    X=q['f']['X'][2,3,4];tangent=1j*old.T[1]@X;k=np.r_[tangent.real,tangent.imag,0.]
    identity=float(np.max(abs(K@k-k/F)))
    assert err<1e-13 and identity<1e-13
    assert q['F_lower']>0 and np.min(q['F'])>=q['F_lower']
    assert np.max(q['U'])<q['U_upper'] and np.min(q['C'])>2-1e-14
    assert np.min(q['Y'])>q['Y_lower']>0 and np.min(q['B'])>=0 and np.min(q['pKp'])>=0
    return dict(M0_squared=M02,inverse_error=err,gauge_tangent_error=identity,
        F_analytic_lower=q['F_lower'],U_analytic_upper=q['U_upper'],Y_analytic_lower=q['Y_lower'],
        tau_squared=q['tau2'],C_min=float(np.min(q['C'])))


def vector_constraint_check():
    q=make_source(16);A,info=solve_momentum(q)
    assert info['momentum_residual']<1e-12 and info['trace_residual']<1e-12
    # Independent real-space finite-difference consistency for the smooth tensor.
    dx=q['dx'];fd=sum((np.roll(A[...,i],-1,axis=i)-np.roll(A[...,i],1,axis=i))/(2*dx) for i in range(3))
    return dict(**info,centered_difference_residual=float(np.max(abs(fd+q['mom']))))


def coupled_constraint_check():
    rows=[]
    for N in (12,16,24):
        q=make_source(N);_,_,info=solve_hamiltonian(q)
        assert info['momentum_residual']<1e-12 and info['original_Hamiltonian_residual']<3e-8
        rows.append(dict(N=N,**info))
    return dict(rows=rows,analytic_existence_is_separate_from_collocation=True)


def uniqueness_and_bounds_check():
    q=make_source(12);a,_,_=solve_hamiltonian(q,.7);b,_,_=solve_hamiltonian(q,1.3)
    err=float(np.max(abs(a-b)));assert err<1e-10
    # Algebraic sign in the maximum-ratio uniqueness argument.
    ratio=1.2;positive=q['C']*(ratio**5-ratio)*a**5+2*q['Y']*(ratio-ratio**-3)*a**-3
    assert np.min(positive)>0
    return dict(two_initial_guesses_max_difference=err,positive_ratio_term_min=float(np.min(positive)))


def maximal_slice_obstruction_check():
    q=make_source(12)
    assert np.min(q['U'])>=-1e-15 and np.min(q['Y'])>0
    # With tau=0, C=-2U<=0: every integrated reaction term is nonpositive
    # and the Yang-Mills term is strictly negative for every positive psi.
    rows=[]
    tensor,_=solve_momentum(q);A=np.sum(tensor*tensor,axis=(-1,-2))+q['pKp']
    for value in (.4,1.,2.):
        reaction=-2*q['U']*value**5-q['B']*value-A*value**-7-2*q['Y']*value**-3
        integral=float(q['dx']**3*np.sum(reaction));assert integral<0
        rows.append(dict(constant_probe=value,integrated_reaction=integral))
    return dict(rows=rows,proof_covers_every_positive_nonconstant_psi_by_periodic_Laplacian_integral=True,
        restriction='same flat conformal CMC branch and nonnegative matched potential; not all GR')


def run():
    checks=('inherited_probe_check','original_source_obstruction_check','counterflow_and_Gauss_check',
            'frame_and_positive_coefficients_check','vector_constraint_check','coupled_constraint_check',
            'uniqueness_and_bounds_check','maximal_slice_obstruction_check')
    evidence={k:globals()[k]() for k in checks}
    deps=('joint_gauss_continuum_sampling.py','joint_gauss_continuum_sampling_results.json',
          'joint_quotient_gauge_completion.py','round572_drafts/momentum_entry_probe.py',
          'round572_drafts/momentum_entry_probe_results.json')
    return dict(round=572,tests_run=len(checks),failures=0,errors=0,checks=list(checks),evidence=evidence,
                dependency_hashes={k:hashlib.sha256((HERE/k).read_bytes()).hexdigest() for k in deps},
                scope='Einstein-frame canonical matching, prescribed four-dimensional two-derivative gravity, flat conformal periodic seed, CMC; original source obstructed, finite changed source admits unique positive conformal factor; no time/quantum/GR-generation theorem')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
