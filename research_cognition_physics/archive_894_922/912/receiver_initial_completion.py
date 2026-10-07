"""912 working: original906 receiver stress ->731/754 compatible initial A.
Actual source coefficient, not a prescribed substitute profile. Numerical
constraint/finite-difference checks, not a rigorous continuum error bound.
"""
from pathlib import Path
import argparse,hashlib,json,sys,time
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(STAGE/'906'));sys.path.insert(0,str(STAGE/'904'))
import compact_receiver_propagation as r
import coupled_boson_tangent as tangent
TARGET=HERE/'receiver_initial_completion_results.json'
def norm(a):return float(np.max(abs(np.asarray(a))))

def initial_stress(model,y):
    z=model.geometry(y);fr=r.shared.frame_data(y,z);N=model.N
    points,rr=r.coordinates(N);d=(points-r.CENTER+np.pi)%(2*np.pi)-np.pi
    psi=r.initial(N);df=np.zeros_like(rr);sel=(rr>r.INNER)&(rr<r.OUTER)
    t=(rr[sel]-r.INNER)/(r.OUTER-r.INNER);f=r.profile(rr[sel])
    df[sel]=-f*(1-f)*(t**-2+(1-t)**-2)/(r.OUTER-r.INNER)
    dx=np.zeros(rr.shape+(3,4),complex);safe=np.where(rr>0,rr,1.)
    dx[..., :,0]=r.NORM*df[...,None]*d/safe[...,None]
    H=r.potential(y,z,fr,psi)
    for i in range(3):
        for a in range(4):
            H-=1j*fr['W'][...,i+1,a,None]*r.mat(r.ALPHA[a],dx[...,i,:])
            H-=.5j*fr['dW'][...,i+1,i+1,a,None]*r.mat(r.ALPHA[a],psi)
    data=r.bilinear_jet(y,z,fr,psi,np.concatenate(((-1j*H)[...,None,:],dx),axis=-2))
    assert norm(data['dirac_residual'])<2e-13
    # alpha=1,beta=0 on the original initial slice.
    assert norm(z['alpha']-1)<1e-14 and norm(z['beta'])==0
    return data,z,data['stress'][...,0,0],z['vol'][...,None]*data['stress'][...,0,1:]

def actual(N,finite_difference=False):
    start=time.time()
    with r.ev.ResearchRuntime(r.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        import joint_source_constraint_response as completion
        model=r.ev.Model(N,old);analytic=tangent.AnalyticModel(N,old);geo=old.geo
        q,p0,_,_,_=old.completed(N,1.);x=q['grid'][...,0]
        base=dict(q);base['B']=q['B']+.0004*np.cos(x)**2;base['C']=q['C']-.0004*np.sin(x)**2;base['U']=q['U']+.0002*np.sin(x)**2
        psi,At,_=geo.solve_hamiltonian(base,initial=p0);y,_=model.initial()
        source,z,rho,J=initial_stress(model,y)
        dp=np.zeros_like(q['p']);dE=np.zeros_like(q['f']['E']);dE0=np.zeros_like(q['f']['E0'])
        dh=np.stack([geo.derivative(q['f']['h'],i) for i in range(3)],axis=-1)
        ds=np.stack([geo.derivative(q['f']['s'],i) for i in range(3)],axis=-1)
        menus=[]
        for i in (0,1):
            pm=np.zeros_like(dp);pm[...,1]=dh[...,i];pm[...,4]=ds[...,i]
            menus.append((pm,np.zeros_like(dE),np.zeros_like(dE0)))
        fw,_=completion.curvature(q);zz=q['grid'][...,2];weighted=np.sin(2*zz)[...,None]*fw[...,0,2,:]
        em=np.zeros_like(dE);em[...,0,:]=completion.covariant(weighted,2,q);em[...,2,:]=-completion.covariant(weighted,0,q)
        menus.append((np.zeros_like(dp),em,np.zeros_like(dE0)))
        columns=np.column_stack([q['dx']**3*completion.momentum(q,*m).sum(axis=(0,1,2)) for m in menus])
        before=q['dx']**3*J.sum(axis=(0,1,2));coeff=np.linalg.solve(columns,-before)
        for c,(pm,em,e0m) in zip(coeff,menus):dp+=c*pm;dE+=c*em;dE0+=c*e0m
        dmom=completion.momentum(q,dp,dE,dE0)+J
        dAt,_=geo.solve_momentum(dict(base,mom=dmom))
        AA=np.sum(At*At,axis=(-1,-2))+q['pKp']
        dAA=2*np.sum(At*dAt,axis=(-1,-2))+2*np.sum(completion.kinverse(q,q['p'])*dp,axis=-1)
        _,bw,b0=geo.old.PAR['b']
        dY=2*bw*np.sum(q['f']['E']*dE,axis=(-1,-2))+2*b0*np.sum(q['f']['E0']*dE0,axis=-1)
        V=5*base['C']*psi**4-base['B']+7*AA*psi**-8+6*base['Y']*psi**-4
        rhs=dAA*psi**-7+2*dY*psi**-3+2*rho*psi**5
        k2=np.sum(geo.waves(N)**2,axis=-1)
        op=lambda v:-8*geo.laplace(v)+V*v
        dpsi,it=geo.cg(op,rhs,lambda v:geo.ifft(geo.fft(v)/(8*k2+V.mean())),tol=2e-13)
        A={k:np.zeros_like(v) for k,v in y.items()}
        dgamma=4*psi[...,None,None]**3*dpsi[...,None,None]*np.eye(3)
        dK=-2*psi[...,None,None]**-3*dpsi[...,None,None]*At+psi[...,None,None]**-2*dAt+np.sqrt(q['tau2'])/3*dgamma
        A['g'][...,1:,1:]=dgamma;A['gd'][...,1:,1:]=-2*dK
        A['gd'][...,0,1:]=A['gd'][...,1:,0]=-2*(model.grad(dpsi)/psi[...,None]-model.grad(psi)*dpsi[...,None]/psi[...,None]**2)
        A['pi'][...,:5]=dp;A['E'][...,8:11]=dE;A['E'][...,11]=dE0
        h=1e-24;raw=tangent.constraints(analytic,{k:v.astype(complex)+1j*h*A[k] for k,v in y.items()})
        linear={k:v.imag/h for k,v in raw.items()}
        harmonic_before=linear['harmonic'].copy()
        #902 stores C_mu (lower index). At alpha=1,beta=0,
        #d C_0 / d gd_00=-1/2, d C_i / d gd_0j=-delta_ij.
        #Only gauge velocities change; intrinsic metric,K,pi,E stay fixed.
        A['gd'][...,0,0]+=2*harmonic_before[...,0]
        A['gd'][...,0,1:]+=harmonic_before[...,1:]
        A['gd'][...,1:,0]+=harmonic_before[...,1:]
        raw=tangent.constraints(analytic,{k:v.astype(complex)+1j*h*A[k] for k,v in y.items()})
        linear={k:v.imag/h for k,v in raw.items()}
        total=dict(H=linear['H']-2*rho,M=linear['M']+J/z['vol'][...,None],Gauss=linear['Gauss'],harmonic=linear['harmonic'])
        finite=[]
        if finite_difference:
            for eps in (.004,.002):
                solutions=[]
                for sign in (1,-1):
                    e=sign*eps;qnew=dict(base);pp=q['p']+e*dp;EE=q['f']['E']+e*dE;EE0=q['f']['E0']+e*dE0
                    qnew['pKp']=np.sum(pp*completion.kinverse(q,pp),axis=-1)
                    qnew['Y']=base['Y']+bw*np.sum(EE*EE-q['f']['E']**2,axis=(-1,-2))+b0*np.sum(EE0*EE0-q['f']['E0']**2,axis=-1)
                    qnew['C']=base['C']-2*e*rho;qnew['mom']=q['mom']+e*dmom
                    solutions.append(geo.solve_hamiltonian(qnew,initial=psi)[0])
                finite.append(dict(source_amplitude=eps,centered_conformal_response_error=norm((solutions[0]-solutions[1])/(2*eps)-dpsi)))
        row=dict(N=N,source_is_original906_mode=True,analytic_initial_mode_spatial_jet_used=True,
            initial_mode_quadrature_norm=r.norm(r.initial(N)),source_Dirac_jet_residual=norm(source['dirac_residual']),
            source_lagrangian_residual=norm(source['lagrangian']),
            original_mode_energy_integral=float(q['dx']**3*np.sum(z['vol']*rho)),
            maximum_normal_energy=norm(rho),maximum_coordinate_momentum_density=norm(J),
            total_coordinate_momentum_before=before.tolist(),counterflow_coefficients=coeff.tolist(),
            total_coordinate_momentum_after=(q['dx']**3*dmom.sum(axis=(0,1,2))).tolist(),
            conformal_linear_equation_residual=norm(op(dpsi)-rhs),
            conformal_momentum_equation_residual=norm(sum(geo.derivative(dAt[..., :,i],i) for i in range(3))+dmom),
            harmonic_before_discrete_gauge_completion=norm(harmonic_before),
            independent_canonical_linear_constraints_with_source={k:norm(v) for k,v in total.items()},
            zero_A_constraint_defect=dict(H=norm(2*rho),M=norm(J/z['vol'][...,None])),
            actual_initial_response_max={k:norm(v) for k,v in A.items()},
            finite_difference_checks=finite,cg_iterations=it,elapsed_seconds=round(time.time()-start,3),
            full_continuum_error_certified=False,full_initial_A_physical_projection_implemented=False,
            source_preparation_is_absolute_renormalized_stress=False)
        assert row['conformal_linear_equation_residual']<2e-10
        assert row['conformal_momentum_equation_residual']<2e-10
        assert norm(row['total_coordinate_momentum_after'])<2e-10
        assert row['independent_canonical_linear_constraints_with_source']['Gauss']<2e-11
        assert row['independent_canonical_linear_constraints_with_source']['harmonic']<2e-10
        return row,y,A

def run():
    rows=[]
    for N in (17,25):
        row,_,_=actual(N,finite_difference=N==17);rows.append(row)
        print(json.dumps(row,ensure_ascii=False),flush=True)
    return dict(round=912,status='working',formal_rounds=911,cumulative_numbered_groups=3696,date='2026-10-06',rows=rows,
        actual_receiver_to_existing_constraint_completion_implemented=True,
        full_A_forced_time_evolution_computed=False,full872_numerical_feedback_completed=False,
        actual_finite_observable_error_certified=False,full_goal_completed=False,
        source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                      (Path(__file__),STAGE/'906/compact_receiver_propagation.py',STAGE/'904/coupled_boson_tangent.py',STAGE/'902/common_background_evolution.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    result=run()
    if a.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
