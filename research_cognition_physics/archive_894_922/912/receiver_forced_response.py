"""912: same906 semidiscrete receiver source, compatible leading initial A,
full904 forced tangent, and explicit product-rule defects. Diagnostics are
not rigorous continuum/physical-observable error certificates.
"""
from pathlib import Path
import argparse,hashlib,json,sys,time
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(HERE))
import receiver_initial_completion as previous
r=previous.r;tangent=previous.tangent
TARGET=HERE/'receiver_forced_response_results.json'

def norm(v):return float(np.max(abs(np.asarray(v))))
def source_data(model,y,u,z=None):
    if z is None:z=model.geometry(y)
    fr=r.shared.frame_data(y,z);ut=-1j*r.apply(model,y,z,fr,u)
    dx=[r.shared.dcomplex(model,u,i) for i in range(3)]
    data=r.bilinear_jet(y,z,fr,u,np.stack([ut,*dx],axis=-2))
    mu=z['alpha']*z['vol'];T=data['stress']
    jg=.5*mu[...,None,None]*np.einsum('...am,...mn,...nb->...ab',z['ig'],T,z['ig'])
    force,load,recovered=r.bridge.convert(y,z,jg,np.zeros_like(y['phi']),np.zeros(y['g'].shape[:3]+(4,12)))
    assert norm(recovered-T)<1e-12
    return data,force,load,ut,z,fr,dx

def product_defect(model,y,u):
    data,_,_,ut,z,fr,dx=source_data(model,y,u)
    localH=r.potential(y,z,fr,u);leibniz=np.zeros_like(u)
    for i in range(3):
        flux=np.zeros_like(u);cdu=np.zeros_like(u);dcu=np.zeros_like(u)
        for a in range(4):
            c=fr['W'][...,i+1,a,None];dc=fr['dW'][...,i+1,i+1,a,None]
            flux+=c*r.mat(r.ALPHA[a],u)
            cdu+=c*r.mat(r.ALPHA[a],dx[i]);dcu+=dc*r.mat(r.ALPHA[a],u)
        leibniz+=r.shared.dcomplex(model,flux,i)-cdu-dcu
        localH-=1j*cdu+.5j*dcu
    localut=-1j*localH;delta=-.5*leibniz
    point=r.bilinear_jet(y,z,fr,u,np.stack([localut,*dx],axis=-2))
    root=np.sqrt(z['vol']);psi=u/root[...,None];deltaR=delta/root[...,None]
    Cup=np.einsum('...mA,Aab->...mab',fr['W']/z['alpha'][...,None,None],r.ALPHA)
    Clow=np.einsum('...mn,...nab->...mab',y['g'],Cup)
    b=np.imag(np.einsum('...a,...mab,...b->...m',psi.conj(),Clow,deltaR))
    dl=-np.imag(np.einsum('...a,...ab,...b->...',psi.conj(),Cup[...,0,:,:],deltaR))
    predicted=y['g']*dl[...,None,None]
    predicted[..., :,0]+=.5*b;predicted[...,0,:]+=.5*b
    dDirac=1j*np.einsum('...ab,...b->...a',Cup[...,0,:,:],deltaR)
    # Cauchy-Schwarz/Frobenius upper bound, evaluated on the finite grid.
    amp=np.linalg.norm(psi,axis=-1)*np.linalg.norm(deltaR,axis=-1)
    cm=np.linalg.norm(Clow,axis=(-2,-1));ct=np.linalg.norm(Cup[...,0,:,:],axis=(-2,-1))
    upper=abs(y['g'])*(amp*ct)[...,None,None]
    upper[..., :,0]+=.5*amp[...,None]*cm;upper[...,0,:]+=.5*amp[...,None]*cm
    assert norm(ut-localut-delta)<1e-12
    assert norm(data['stress']-point['stress']-predicted)<1e-12
    assert norm(data['dirac_residual']-point['dirac_residual']-dDirac)<1e-12
    assert np.max(abs(predicted)-upper)<1e-12
    result=dict(time_jet_identity_error=norm(ut-localut-delta),
        source_identity_error=norm(data['stress']-point['stress']-predicted),
        Dirac_identity_error=norm(data['dirac_residual']-point['dirac_residual']-dDirac),
        product_defect_max=norm(leibniz),point_Dirac_defect=norm(point['dirac_residual']),
        split_Dirac_defect=norm(data['dirac_residual']),source_correction_max=norm(predicted),
        finite_grid_source_correction_upper_max=norm(upper),
        continuum_source_error_certified=False)
    if norm(u[...,1:])<1e-14 and norm(u[...,0].imag)<1e-14:
        f=u[...,0].real;s=np.array([0.,0.,1.]);df=model.grad(f)
        curl=-.25*np.cross(model.grad(f*f),s)
        pointJ=z['vol'][...,None]*point['stress'][...,0,1:]
        splitJ=z['vol'][...,None]*data['stress'][...,0,1:]
        correction=.25*np.cross(model.grad(f*f)-2*f[...,None]*df,s)+z['vol'][...,None]*predicted[...,0,1:]
        result.update(initial_spin_curl_decomposition_error=norm(splitJ-curl-correction),
            point_momentum_identity_error=norm(pointJ+.5*f[...,None]*np.cross(df,s)),
            split_to_curl_max=norm(splitJ-curl),curl_zero_mode=norm(curl.mean(axis=(0,1,2))),
            extra_physical_counterflow_required=False)
        assert result['initial_spin_curl_decomposition_error']<1e-12
    return result

def initial(model,analytic,old):
    geo=old.geo;N=model.N;q,p0,_,_,_=old.completed(N,1.);x=q['grid'][...,0]
    base=dict(q);base['B']=q['B']+.0004*np.cos(x)**2;base['C']=q['C']-.0004*np.sin(x)**2;base['U']=q['U']+.0002*np.sin(x)**2
    psi,At,_=geo.solve_hamiltonian(base,initial=p0);y,_=model.initial();u=r.initial(N)
    data,_,loads,_,z,_,_=source_data(model,y,u)
    rho=-loads['H']/2;J=z['vol'][...,None]*loads['M'];mean=J.mean(axis=(0,1,2))
    # Original continuous source has zero total momentum (912 spin-curl).
    # Solve nonzero numerical modes, expose the zero-mode defect; do not
    # turn a quadrature/Leibniz artifact into an additional physical flow.
    dAt,_=geo.solve_momentum(dict(base,mom=J))
    AA=np.sum(At*At,axis=(-1,-2))+q['pKp']
    dAA=2*np.sum(At*dAt,axis=(-1,-2))
    V=5*base['C']*psi**4-base['B']+7*AA*psi**-8+6*base['Y']*psi**-4
    rhs=dAA*psi**-7+2*rho*psi**5;k2=np.sum(geo.waves(N)**2,axis=-1)
    op=lambda v:-8*geo.laplace(v)+V*v
    dpsi,it=geo.cg(op,rhs,lambda v:geo.ifft(geo.fft(v)/(8*k2+V.mean())),tol=2e-13)
    A={k:np.zeros_like(v) for k,v in y.items()}
    dgam=4*psi[...,None,None]**3*dpsi[...,None,None]*np.eye(3)
    dK=-2*psi[...,None,None]**-3*dpsi[...,None,None]*At+psi[...,None,None]**-2*dAt+np.sqrt(q['tau2'])/3*dgam
    A['g'][...,1:,1:]=dgam;A['gd'][...,1:,1:]=-2*dK
    A['gd'][...,0,1:]=A['gd'][...,1:,0]=-2*(model.grad(dpsi)/psi[...,None]-model.grad(psi)*dpsi[...,None]/psi[...,None]**2)
    h=1e-24;hc=tangent.constraints(analytic,{k:v.astype(complex)+1j*h*A[k] for k,v in y.items()})['harmonic'].imag/h
    A['gd'][...,0,0]+=2*hc[...,0];A['gd'][...,0,1:]+=hc[...,1:];A['gd'][...,1:,0]+=hc[...,1:]
    momentum_residual=sum(geo.derivative(dAt[..., :,i],i) for i in range(3))+J
    result=dict(N=N,source_provider='same906_Hermitian_semidiscrete_at_all_times',
        additional_matter_counterflow_max=0.,numerical_momentum_zero_mode=mean.tolist(),
        total_coordinate_momentum=((2*np.pi)**3*mean).tolist(),
        conformal_residual=norm(op(dpsi)-rhs),conformal_momentum_residual=norm(momentum_residual),
        momentum_residual_minus_zero_mode=norm(momentum_residual-mean),cg_iterations=it,
        physical_constraints_certified=False)
    assert result['conformal_residual']<2e-10 and result['momentum_residual_minus_zero_mode']<2e-11
    return y,u,A,result

def total_constraints(model,analytic,y,u,A):
    h=1e-24;c=tangent.constraints(analytic,{k:v.astype(complex)+1j*h*A[k] for k,v in y.items()})
    _,_,load,_,_,_,_=source_data(model,y,u)
    out={k:v.imag/h for k,v in c.items()}
    out['H']+=load['H'];out['M']+=load['M'];out['Gauss']+=load['Gauss']
    return {k:norm(v) for k,v in out.items()}

def rhs(model,analytic,y,u,A):
    dy,z,_=model.rhs(y,aux=True);_,force,_,ut,_,_,_=source_data(model,y,u,z)
    LA=analytic.jvp(y,A)
    return dy,ut,{k:LA[k]+force[k] for k in A}

def step(model,analytic,y,u,A,dt):
    stages=[]
    for c,index in ((0,None),(.5,0),(.5,1),(1,2)):
        yy=y if index is None else r.ev.add(y,stages[index][0],c*dt)
        uu=u if index is None else u+c*dt*stages[index][1]
        AA=A if index is None else r.ev.add(A,stages[index][2],c*dt)
        stages.append(rhs(model,analytic,yy,uu,AA))
    weights=(1,2,2,1)
    return ({k:y[k]+dt/6*sum(w*s[0][k] for w,s in zip(weights,stages)) for k in y},
        u+dt/6*sum(w*s[1] for w,s in zip(weights,stages)),
        {k:A[k]+dt/6*sum(w*s[2][k] for w,s in zip(weights,stages)) for k in A})

def finite_amplitude_rhs(model,y,u,eps):
    # Independent Einstein equation with T_total=T_boson+eps*T_receiver.
    # This validates the finite ODE tangent, not nonlinear constraint completion.
    dy,z,m=model.rhs(y,aux=True);data,_,_,ut,_,_,_=source_data(model,y,u,z)
    T=m['T']+eps*data['stress'];trace=np.einsum('...mn,...mn->...',z['ig'],T)
    wave=2*z['Q']-2*T+y['g']*trace[...,None,None]
    wave-=2*np.einsum('...i,...imn->...mn',z['ig'][...,0,1:],model.grad(y['gd']))
    wave-=model.lapweighted(y['g'],z['ig'][...,1:,1:])
    dy['gd']=wave/z['ig'][...,0,0,None,None]
    return dy,ut

def finite_step(model,y,u,eps,dt):
    s=[]
    for c,j in ((0,None),(.5,0),(.5,1),(1,2)):
        yy=y if j is None else r.ev.add(y,s[j][0],c*dt)
        uu=u if j is None else u+c*dt*s[j][1]
        s.append(finite_amplitude_rhs(model,yy,uu,eps))
    weights=(1,2,2,1)
    return ({k:y[k]+dt/6*sum(w*t[0][k] for w,t in zip(weights,s)) for k in y},
        u+dt/6*sum(w*t[1] for w,t in zip(weights,s)))

def flow(N,steps=2,Tend=.00025,finite_difference=False):
    start=time.time()
    with r.ev.ResearchRuntime(r.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=r.ev.Model(N,old);analytic=tangent.AnalyticModel(N,old)
        y,u,A,init=initial(model,analytic,old);y0=y;u0=u;A0=A
        defect0=product_defect(model,y,u);c0=total_constraints(model,analytic,y,u,A);n0=r.norm(u)
        for _ in range(steps):y,u,A=step(model,analytic,y,u,A,Tend/steps)
        cf=total_constraints(model,analytic,y,u,A);df=product_defect(model,y,u)
        fd=[]
        if finite_difference:
            for eps in (.1,.05):
                sides=[]
                for sign in (1,-1):
                    yy=r.ev.add(y0,A0,sign*eps);uu=u0.copy()
                    for _ in range(steps):yy,uu=finite_step(model,yy,uu,sign*eps,Tend/steps)
                    sides.append(yy)
                errors={k:norm((sides[0][k]-sides[1][k])/(2*eps)-A[k]) for k in A}
                fd.append(dict(source_amplitude=eps,all_components_error=errors,maximum_error=max(errors.values())))
            assert fd[-1]['maximum_error']<1e-6, fd
        row=dict(N=N,steps=steps,Tend=Tend,initial=init,initial_defect=defect0,final_defect=df,
            initial_total_constraints=c0,final_total_constraints=cf,
            initial_response_max={k:norm(v) for k,v in A0.items()},final_response_max={k:norm(v) for k,v in A.items()},
            response_change_max={k:norm(A[k]-A0[k]) for k in A},receiver_norm_drift=abs(r.norm(u)-n0),
            independent_finite_amplitude_checks=fd,elapsed_seconds=round(time.time()-start,3),
            all104_mixed_boson_variables_retained=True,same_receiver_source_initial_and_evolution=True,
            physical_Pi_projection_implemented=False,continuum_error_certified=False,
            total872_response_computed=False)
        print(json.dumps(row,ensure_ascii=False),flush=True)
        return row,y,u,A

def run():
    rows=[];fields={}
    for N,steps,Tend,fd in ((17,2,.00025,True),(17,4,.00025,False),(17,2,-.00025,False),(25,2,.00025,False)):
        row,y,u,A=flow(N,steps,Tend,fd);rows.append(row);fields[(N,steps,Tend)]=(y,u,A)
    a=fields[(17,2,.00025)];b=fields[(17,4,.00025)]
    diff=dict(background={k:norm(a[0][k]-b[0][k]) for k in a[0]},receiver=norm(a[1]-b[1]),
        forced_A={k:norm(a[2][k]-b[2][k]) for k in a[2]})
    hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        (Path(__file__),HERE/'receiver_initial_completion.py',STAGE/'906/compact_receiver_propagation.py',STAGE/'904/coupled_boson_tangent.py',STAGE/'905/canonical_source_bridge.py')}
    return dict(round=912,date='2026-10-06',rows=rows,timestep_diagnostic=diff,
        same_semidiscrete_source_used_throughout=True,full104_forced_tangent_computed=True,
        product_rule_source_defect_identity_checked=True,zero_mode_defect_exposed_not_compensated=True,
        actual_finite_observable_error_certified=False,full872_numerical_feedback_completed=False,
        full_goal_completed=False,source_hashes=hashes)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not TARGET.exists()
    result=run()
    if args.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'timestep_diagnostic':result['timestep_diagnostic'],'saved':args.write},indent=2))
