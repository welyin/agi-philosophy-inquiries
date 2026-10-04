"""571: exact finite-graph Gauss data from a declared smooth continuum source.

Fixed periodic geometry and nowhere-zero Higgs in one global unitary gauge.
This is an initial-data map, not a flow/quantum/geometry-generation theorem.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_gauss_nonlinear_compatibility as prev
import joint_quotient_gauge_completion as old
import joint_scalar_propagation_matching as scalar

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauss_continuum_sampling_results.json'
LBOX=2*np.pi
T=old.generators(2)
PAR=old.parameters()
A0,B0,U0,AMP=.41,.29,.38,.17


def rotate(a,v,factor):
    """Ordinary SO(3) rotation with generator +a cross; Ad(exp(i a.t)) is minus."""
    length=np.linalg.norm(a,axis=-1,keepdims=True)
    axis=np.divide(a,length,out=np.zeros_like(a),where=length>0)
    theta=factor*length
    return v*np.cos(theta)+np.cross(axis,v)*np.sin(theta)+axis*np.sum(axis*v,axis=-1,keepdims=True)*(1-np.cos(theta))


def su2(a,factor):
    length=np.linalg.norm(a,axis=-1)
    ratio=np.divide(np.sin(factor*length/2),length,
                    out=np.full_like(length,factor/2),where=length>0)
    return np.cos(factor*length/2)[...,None,None]*np.eye(2)+2j*ratio[...,None,None]*np.einsum('...a,aij->...ij',a,T)


def fields(points,mode='charged'):
    d=points.shape[-1];shape=points.shape[:-1];x,y=points[...,0],points[...,1]
    f=1+U0*np.cos(x);cc=1+1.5*U0**2
    a=np.zeros(shape+(d,3));E=np.zeros_like(a);aa=np.zeros(shape+(d,));EE=np.zeros_like(aa)
    r=np.zeros(shape+(3,))
    if mode=='charged':
        a[...,0,0]=A0*f;a[...,1,1]=B0*(1+.2*np.sin(y))
        E[...,0,1]=AMP*(f*f-cc);E[...,1,1]=.11*np.cos(x)
        r[...,1]=-2*AMP*U0*f*np.sin(x);r[...,2]=-A0*f*E[...,0,1]
        EE[...,0]=6*A0*AMP*((2*U0-.75*U0**3)*np.sin(x)+.75*U0**2*np.sin(2*x)+U0**3*np.sin(3*x)/12)+.13
        if d==3:
            a[...,2,2]=.12*np.cos(points[...,2]);E[...,2,2]=.09*np.sin(y)
    else:
        a[...,0,2]=.17*f;a[...,1,2]=.13*(1+.2*np.sin(y))
        E[...,0,2]=.04*np.sin(x);E[...,1,2]=.11*np.cos(x)
        r[...,2]=.04*np.cos(x);EE[...,0]=-.24*np.sin(x)+.13
    EE[...,1]=.07*np.cos(x)
    aa[...,0]=.03*np.cos(y);aa[...,1]=.08*np.sin(x)
    h=np.sqrt(PAR['h2'])*(1+.05*np.sin(x+y))
    _,u,_=scalar.parameters();s=np.sqrt(u[1])*(1+.06*np.cos(y))
    X=np.zeros(shape+(2,),complex);X[...,1]=h
    TX=1j*np.einsum('aij,...j->...ai',T,X)
    radial=.04*np.cos(x+y)
    PX=np.zeros_like(X);PX[...,1]=radial
    PX-=4/h[...,None]**2*np.einsum('...a,...ai->...i',r,TX)
    return dict(a=a,E=E,a0=aa,E0=EE,h=h,X=X,s=s,PX=PX,Ps=.025*np.cos(x),
                r=r,radial=radial)


def charges(X,P):
    TX=1j*np.einsum('aij,...j->...ai',T,X)
    qw=np.einsum('...i,...ai->...a',P.conj(),TX).real
    q0=np.sum(P.conj()*(3j*X),axis=-1).real
    return qw,q0


def incidence_transpose(edge):
    d=edge.shape[-1]
    return sum(edge[...,mu]-np.roll(edge[...,mu],1,axis=mu) for mu in range(d))


def gradient(node):
    return np.stack([node-np.roll(node,-1,axis=mu) for mu in range(node.ndim)],axis=-1)


def covdiv(a,p,eps):
    d=a.shape[-2]
    return sum(p[...,mu,:]-np.roll(rotate(a[...,mu,:],p[...,mu,:],eps),1,axis=mu)
               for mu in range(d))


def poisson_gradient(rho):
    N=rho.shape[0];d=rho.ndim
    lam=np.zeros(rho.shape)
    for mu in range(d):
        shape=[1]*d;shape[mu]=N
        lam+=4*np.sin(np.pi*np.fft.fftfreq(N)).reshape(shape)**2
    spectrum=np.fft.fftn(rho);safe=lam.copy();safe[(0,)*d]=1.
    potential=spectrum/safe;potential[(0,)*d]=0
    return gradient(np.fft.ifftn(potential).real)


def sample(N,d=2,mode='charged'):
    eps=LBOX/N;grid=np.stack(np.meshgrid(*([np.arange(N)*eps]*d),indexing='ij'),axis=-1)
    f=fields(grid,mode);a=np.empty(grid.shape[:-1]+(d,3));E=np.empty_like(a)
    a0=np.empty(grid.shape[:-1]+(d,));E0=np.empty_like(a0)
    for mu in range(d):
        mid=grid.copy();mid[...,mu]+=eps/2;fm=fields(mid,mode)
        a[...,mu,:]=fm['a'][...,mu,:];E[...,mu,:]=fm['E'][...,mu,:]
        a0[...,mu]=fm['a0'][...,mu];E0[...,mu]=fm['E0'][...,mu]
    p=eps**(d-1)*rotate(a,E,-eps/2);p0=eps**(d-1)*E0;PX=eps**d*f['PX']
    c=np.array([0.,0.,1.])-rotate(a,np.broadcast_to([0.,0.,1.],a.shape),-eps)
    M=float(np.sum(c*c));S=float(np.sum(c*p))
    dp=np.zeros_like(p) if mode=='neutral' else -(S/M)*c
    assert mode!='neutral' or (M==0 and abs(S)<1e-14)
    phat=p+dp;rhat=covdiv(a,phat,eps)
    qw,_=charges(f['X'],PX);weak_residual=qw+rhat
    TX=1j*np.einsum('aij,...j->...ai',T,f['X'])
    dPX=-4/f['h'][...,None]**2*np.einsum('...a,...ai->...i',weak_residual,TX)
    rho=-6*rhat[...,2]-incidence_transpose(p0)
    dp0=poisson_gradient(rho)
    return dict(N=N,d=d,eps=eps,grid=grid,f=f,a=a,a0=a0,p=p,p0=p0,PX=PX,
                dp=dp,dp0=dp0,dPX=dPX,c=c,M=M,S=S,rho=rho,
                phat=phat,p0hat=p0+dp0,PXhat=PX+dPX)


def norm_squared(q,P,p,p0):
    eps,d=q['eps'],q['d'];bw,b0=PAR['b'][1:]
    return float(np.sum(abs(P)**2)/eps**d+2*eps**(2-d)*(bw*np.sum(p*p)+b0*np.sum(p0*p0)))


def residual(q,corrected=True):
    P=q['PXhat'] if corrected else q['PX'];p=q['phat'] if corrected else q['p']
    p0=q['p0hat'] if corrected else q['p0']
    qw,q0=charges(q['f']['X'],P)
    return qw+covdiv(q['a'],p,q['eps']),q0+incidence_transpose(p0)


def energy(q,corrected=True):
    eps,d=q['eps'],q['d'];f=q['f'];X=f['X'];s=f['s']
    P,p,p0=(q['PXhat'],q['phat'],q['p0hat']) if corrected else (q['PX'],q['p'],q['p0'])
    kinetic=norm_squared(q,P,p,p0)/2+eps**d*np.sum(f['Ps']**2)/2
    mat,u,_=scalar.parameters();delta=np.stack((f['h']**2-u[0],s*s-u[1]),axis=-1)
    potential=eps**d*np.sum(np.einsum('...i,ij,...j->...',delta,mat,delta))/4
    W=su2(q['a'],eps);z=np.exp(1j*eps*q['a0'])
    for mu in range(d):
        neighbour=np.roll(X,-1,axis=mu)
        difference=X-z[...,mu,None]**3*np.einsum('...ij,...j->...i',W[...,mu,:,:],neighbour)
        potential+=eps**(d-2)*np.sum(abs(difference)**2)/2
        potential+=eps**(d-2)*np.sum((s-np.roll(s,-1,axis=mu))**2)/2
    dg,coef=old.coefficients(PAR)
    for mu in range(d):
        for nu in range(mu+1,d):
            wm,wn=W[...,mu,:,:],W[...,nu,:,:]
            wc=wm@np.roll(wn,-1,axis=mu)@np.swapaxes(np.roll(wm,-1,axis=nu).conj(),-1,-2)@np.swapaxes(wn.conj(),-1,-2)
            zc=z[...,mu]*np.roll(z[...,nu],-1,axis=mu)/np.roll(z[...,mu],-1,axis=nu)/z[...,nu]
            tw=np.trace(wc,axis1=-2,axis2=-1)
            chi=PAR['wq']*(3*tw*zc+3*(zc**-4+zc**2))+PAR['wl']*(tw*zc**-3+zc**6+1)
            vf=dg*(12*PAR['wq']+4*PAR['wl']-chi.real)+coef[1]*(4-abs(tw)**2)+coef[2]*(1-(zc**6).real)
            potential+=eps**(d-4)*np.sum(vf)
    return float(kinetic+potential)


def continuum_source_check():
    rng=np.random.default_rng(571);points=rng.uniform(0,LBOX,size=(31,3));step=2e-5;rows=[]
    for mode in ('charged','neutral'):
        f=fields(points,mode);div=np.zeros((31,3));div0=np.zeros(31)
        for mu in range(3):
            shift=np.zeros(3);shift[mu]=step
            plus,minus=fields(points+shift,mode),fields(points-shift,mode)
            div+=(plus['E'][:,mu]-minus['E'][:,mu])/(2*step)
            div0+=(plus['E0'][:,mu]-minus['E0'][:,mu])/(2*step)
        rr=div-np.sum(np.cross(f['a'],f['E']),axis=1)
        qw,q0=charges(f['X'],f['PX'])
        err=max(np.max(abs(qw+rr)),np.max(abs(q0+div0)),np.max(abs(rr-f['r'])))
        assert err<1e-9
        rows.append(dict(mode=mode,independent_finite_difference_Gauss_error=float(err)))
    return dict(rows=rows)


def half_transport_check():
    rng=np.random.default_rng(5711);worst=0.
    for _ in range(10):
        a,v=rng.normal(size=(2,3));U=old.group_exp(.23*a,2)
        got=prev.adjoint(U)@v;expected=rotate(a,v,-.23)
        worst=max(worst,float(np.max(abs(got-expected))))
    rows=[]
    for N in (8,12,18,26):
        q=sample(N);v=q['eps']**2
        diff=covdiv(q['a'],q['p'],q['eps'])/v-q['f']['r']
        err=float(np.sqrt(np.mean(np.sum(diff**2,axis=-1))))
        rows.append(dict(N=N,error=err,error_over_eps_squared=err/q['eps']**2))
    assert worst<1e-13 and rows[-1]['error']<rows[0]['error']/8
    return dict(adjoint_matrix_vs_rotation_error=worst,rows=rows)


def exact_sampling_check():
    rows=[]
    for d,Ns in ((2,(8,12,18,26)),(3,(8,12,16))):
        for N in Ns:
            q=sample(N,d);gw,g0=residual(q)
            after=max(np.max(abs(gw)),np.max(abs(g0)))/q['eps']**d
            before=max(np.max(abs(x)) for x in residual(q,False))/q['eps']**d
            correction=np.sqrt(norm_squared(q,q['dPX'],q['dp'],q['dp0']))
            assert after<2e-11 and abs(q['S'])>1e-5
            rows.append(dict(d=d,N=N,initial_Gauss_density_error=float(before),
                final_Gauss_density_error=float(after),sample_total_weak3=q['S'],
                correction_kinetic_norm=float(correction),
                correction_norm_over_eps_squared=float(correction/q['eps']**2)))
    return dict(rows=rows)


def source_degrees_check():
    q=sample(18);eps=q['eps'];h=q['f']['h']
    radial=np.sum(q['f']['X'].conj()*q['dPX'],axis=-1).real
    # Any periodic curl plus a harmonic edge field is orthogonal to a node gradient.
    psi=np.sin(q['grid'][...,0]+2*q['grid'][...,1]);probe=np.zeros_like(q['p0'])
    probe[...,0]=psi-np.roll(psi,1,axis=1)
    probe[...,1]=np.roll(psi,1,axis=0)-psi
    probe+=np.array([.13,-.21])
    inner=float(np.sum(probe*q['dp0']))
    harm=float(np.max(abs(np.mean(q['dp0'],axis=(0,1)))))
    assert np.max(abs(radial))<1e-14 and abs(inner)<1e-12 and harm<1e-14
    assert np.max(abs(incidence_transpose(probe)))<1e-14
    gamma=eps**0*q['M'];expected=LBOX**2*(A0**2*(1+U0**2/2)+B0**2*1.02)
    return dict(radial_momentum_change=float(np.max(abs(radial))),transverse_probe_pairing=inner,
        harmonic_momentum_change=harm,normalized_gap=gamma,continuum_gap=expected,
        min_Higgs_amplitude=float(np.min(h)))


def neutral_branch_check():
    rows=[]
    for N in (8,16,24):
        q=sample(N,mode='neutral');gw,g0=residual(q);v=q['eps']**2
        err=max(np.max(abs(gw)),np.max(abs(g0)))/v
        assert q['M']==0 and np.max(abs(q['dp']))==0 and err<1e-11
        corr=np.sqrt(norm_squared(q,q['dPX'],q['dp'],q['dp0']))
        rows.append(dict(N=N,gap=q['M'],Gauss_density_error=float(err),
                         correction_norm_over_eps_squared=float(corr/q['eps']**2)))
    return dict(rows=rows,zero_gap_branch_is_supported_without_division=True)


def global_balance_check():
    rows=[]
    for N in (8,12,18,26):
        q=sample(N);eps=q['eps'];c=q['c'];r=covdiv(q['a'],q['p'],eps)
        identity=abs(float(np.sum(r[...,2]))-q['S'])
        after=float(np.sum(c*q['phat']))
        gap=eps**0*q['M']
        predicted=LBOX**2*(A0**2*(1+U0**2/2)+B0**2*1.02)
        rho_mean=float(np.mean(q['rho']))
        dp0=q['dp0'];bound=np.linalg.norm(q['rho'])/(2*np.sin(np.pi/N))
        assert identity<2e-14 and abs(after)<2e-14 and abs(rho_mean)<1e-14
        assert np.linalg.norm(dp0)<=bound+1e-12
        rows.append(dict(N=N,total_charge_identity_error=identity,total_charge_after=after,
            normalized_gap=gap,continuum_gap=predicted,S_over_eps_squared=q['S']/eps**2,
            circle_correction_norm=float(np.linalg.norm(dp0)),Poisson_bound=float(bound)))
    return dict(rows=rows)


def full_energy_check():
    rows=[]
    for N in (8,12,18,26):
        q=sample(N);raw=energy(q,False);fixed=energy(q,True)
        norm=np.sqrt(norm_squared(q,q['PX'],q['p'],q['p0']))
        delta=np.sqrt(norm_squared(q,q['dPX'],q['dp'],q['dp0']))
        bound=norm*delta+delta**2/2
        assert abs(fixed-raw)<=bound+1e-11 and min(raw,fixed)>0
        rows.append(dict(N=N,full_trial_H=raw,full_physical_H=fixed,
                         energy_change=fixed-raw,proven_change_bound=float(bound)))
    return dict(rows=rows,full_569_positive_magnetic_function_used=True)


def discarded_source_counterexample():
    rows=[]
    for N in (8,16,24):
        q=sample(N);eps=q['eps'];r=covdiv(q['a'],q['phat'],eps)
        rebuilt=poisson_gradient(-6*r[...,2])
        dropped=q['p0hat']-rebuilt
        norm=np.sqrt(2*PAR['b'][2]*np.sum(dropped*dropped))
        exact=np.sqrt(2*PAR['b'][2]*LBOX**2*(.13**2+.07**2/2))
        radial=eps**2*q['f']['radial']
        radial_norm=np.sqrt(np.sum(radial**2)/eps**2)
        expected_radial=.04*LBOX/np.sqrt(2)
        assert abs(norm-exact)<1e-12 and abs(radial_norm-expected_radial)<1e-12
        rows.append(dict(N=N,photon_information_discarded_norm=float(norm),
                         radial_information_discarded_norm=float(radial_norm)))
    return dict(rows=rows,
        Gauss_alone_does_not_authorize_discarding_transverse_harmonic_or_radial_source=True)


def run():
    names=('continuum_source_check','half_transport_check','exact_sampling_check','source_degrees_check',
           'neutral_branch_check','global_balance_check','full_energy_check','discarded_source_counterexample')
    evidence={name:globals()[name]() for name in names}
    deps=('joint_gauss_nonlinear_compatibility.py','joint_quotient_gauge_completion.py',
          'joint_scalar_propagation_matching.py','joint_singlet_common_mass_rg_results.json')
    return dict(round=571,tests_run=len(names),failures=0,errors=0,checks=list(names),evidence=evidence,
                dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
                scope='fixed smooth periodic classical initial data; nonzero Higgs in a global unitary gauge; fixed positive charged-background gap or exactly neutral branch; no time/quantum/gravity convergence')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
