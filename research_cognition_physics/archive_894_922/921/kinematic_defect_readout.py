"""921: original Legendre tangent residual -> full104 response -> original readout.
This computes one explicitly isolated residual contribution at the same fixed
approximate background; it is NOT a full background/tangent error certificate.
"""
from pathlib import Path
from types import SimpleNamespace
import sys,argparse,json,time,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(STAGE/'920'))
import common_legendre_defect as defect
prior=defect.prior;weak=prior.weak;c=defect.c;response=defect.response;r=defect.r;mb=defect.mb
TARGET=HERE/'kinematic_defect_readout_results.json'

def coefficients(minus,plus,T):
    ym,dm=minus;yp,dp=plus
    return np.stack(((yp+ym)/2-T*(dp-dm)/4,3*(yp-ym)/4-T*(dp+dm)/4,T*(dp-dm)/4,(T*(dp+dm)-(yp-ym))/4),axis=0)

class Family:
    def __init__(self,T,old):
        self.T=T;self.model=r.ev.Model(17,old);self.analytic=response.tangent.AnalyticModel(17,old)
        y0,u0,A0,_=response.initial(self.model,self.analytic,old);ends=[]
        for sign in (-1,1):
            y={k:v.copy() for k,v in y0.items()};u=u0.copy();A={k:v.copy() for k,v in A0.items()}
            for _ in range(2):y,u,A=response.step(self.model,self.analytic,y,u,A,sign*T/2)
            dy,_,dA=response.rhs(self.model,self.analytic,y,u,A);ends.append((y,dy,A,dA))
        self.co={k:coefficients((ends[0][0][k],ends[0][1][k]),(ends[1][0][k],ends[1][1][k]),T) for k in y0}
        pack=lambda y:np.concatenate((y['pi'],y['E'].reshape(y['pi'].shape[:3]+(36,))),axis=-1)
        mc=defect.MomentumField((pack(ends[0][0]),pack(ends[0][1])),(pack(ends[1][0]),pack(ends[1][1])),self.model.k,T)
        mv=defect.MomentumField((pack(ends[0][2]),pack(ends[0][3])),(pack(ends[1][2]),pack(ends[1][3])),self.model.k,T)
        self.mc=SimpleNamespace(k=mc.k,c=mc.co,T=T);self.mv=SimpleNamespace(k=mv.k,c=mv.co,T=T)
        self.par=SimpleNamespace(K=self.model.K,L=self.model.L,u=self.model.u)
    def y(self,t):
        s=t/self.T;w=np.array([1.,s,s*s,s*s*s]);return {k:np.einsum('p,p...->...',w,v) for k,v in self.co.items()}

def first_jets(grid):
    g,phi,A=mb.unpack(grid.get());_,dphi,dA=mb.unpack(np.stack([grid.get((mu,)) for mu in range(4)],axis=1))
    return dict(g=g,phi=phi,A=A,dphi=dphi,dA=dA)

def kinetic_residual(p,packed,par):
    g=p['g'];ig=np.linalg.inv(g);gamma=g[:,1:,1:];vol=np.sqrt(np.linalg.det(gamma));alpha=1/np.sqrt(-ig[:,0,0]);beta=-ig[:,0,1:]/ig[:,0,0,None]
    G,*_=c.target(p['phi'],par);Gi=np.linalg.inv(G);pi=packed[:,:6];E=packed[:,6:].reshape((-1,3,12))
    A=p['A'][:,1:,:];D=p['dphi'][:,1:,:]+np.einsum('bia,aAB,bB->biA',A,c.REP,p['phi'])
    dA=p['dA'][:,1:,1:,:];F=dA-dA.swapaxes(1,2)+c.bracket(A[:,:,None,:],A[:,None,:,:])
    phidot=alpha[:,None]/vol[:,None]*np.einsum('bAB,bB->bA',Gi,pi)+np.einsum('bi,biA->bA',beta,D)
    Adot=alpha[:,None,None]/vol[:,None,None]*np.einsum('bij,bja->bia',gamma,E)/par.K+np.einsum('bj,bjia->bia',beta,F)
    return dict(phi=p['dphi'][:,0,:]-phidot,A=p['dA'][:,0,1:,:]-Adot)

def projected_residual(fam,bg,field,t,M):
    p=first_jets(weak.GridField(bg,t,M));v=first_jets(weak.GridField(field,t,M))
    mc=weak.GridField(fam.mc,t,M).get();mv=weak.GridField(fam.mv,t,M).get();h=1e-24
    base=kinetic_residual(p,mc,fam.par)
    shifted=kinetic_residual({k:p[k].astype(complex)+1j*h*v[k] for k in p},mc.astype(complex)+1j*h*mv,fam.par)
    tangent={k:a.imag/h for k,a in shifted.items()};N=fam.model.N
    ik=np.rint(fam.model.k).astype(int);ids=tuple(ik[...,i]%M for i in range(3))
    def project(a):
        co=np.fft.fftn(a.reshape((M,M,M)+a.shape[1:]),axes=(0,1,2))/M**3
        return np.fft.ifftn(co[ids]*N**3,axes=(0,1,2)).real
    force={k:np.zeros_like(v) for k,v in fam.y(0.).items()}
    for k,v in tangent.items():force[k]=-project(v)
    stats=dict(t=t,M=M,base_residual_sample_max={k:c.maximum(a) for k,a in base.items()},tangent_residual_sample_max={k:c.maximum(a) for k,a in tangent.items()},
       projected_tangent_force_max={k:c.maximum(force[k]) for k in ('phi','A')})
    return force,stats

def interpolate_force(values,t,T):
    s=t/T;weights=np.array([s*(s-1)/2,1-s*s,s*(s+1)/2])
    return {k:sum(w*a[k] for w,a in zip(weights,values)) for k in values[0]}

def integrate(fam,forcing):
    T=fam.T;ends=[];stats=[]
    def rhs(t,u):
        L=fam.analytic.jvp(fam.y(t),u);f=interpolate_force(forcing,t,T)
        return {k:L[k]+f[k] for k in u}
    for sign in (-1,1):
        u={k:np.zeros_like(a) for k,a in fam.y(0.).items()};dt=sign*T/2;t=0.
        for _ in range(2):
            stages=[]
            for d,j in ((0,None),(.5,0),(.5,1),(1.,2)):
                trial=u if j is None else {k:u[k]+d*dt*stages[j][k] for k in u}
                stages.append(rhs(t+d*dt,trial))
            u={k:u[k]+dt/6*sum(w*s[k] for w,s in zip((1,2,2,1),stages)) for k in u};t+=dt
        du=rhs(t,u);ends.append((mb.pack(u),mb.pack(du)))
        stats.append(dict(time=t,all104_response_max={k:c.maximum(a) for k,a in u.items()},
          linear_canonical_constraints=response.tangent.tangent_constraints(fam.analytic,fam.y(t),u)))
    field=c.ex.previous.rebuilt.work.interpolant(17,T,ends[0],ends[1],fam.model.k,0.)
    return field,stats

def run():
    start=time.time();bg,original,_=c.ex.previous.rebuilt.build_pair();T=bg.T
    with r.ev.ResearchRuntime(r.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        fam=Family(T,old);forcing=[];force_stats=[]
        for t in (-T,0.,T):
            f,s=projected_residual(fam,bg,original,t,25);forcing.append(f);force_stats.append(s)
            print(json.dumps(dict(stage='force',stats=s,elapsed=round(time.time()-start,2))),flush=True)
        nodal,nodalstats=projected_residual(fam,bg,original,0.,17)
        collocation_gap={k:c.maximum(nodal[k]-forcing[1][k]) for k in ('phi','A')}
        checks=[]
        for t in (-T/2,T/2):
            actual,stats=projected_residual(fam,bg,original,t,25);f=interpolate_force(forcing,t,T)
            checks.append(dict(t=t,interpolation_difference_not_bound={k:c.maximum(actual[k]-f[k]) for k in ('phi','A')}))
        correction,evolution=integrate(fam,forcing)
    src=c.ex.loop.LoopSource(bg,2,8);src.kernels()
    pair=src.response(lambda x,p:correction.jets(x))
    direct,details=c.aj.previous.exact_record_derivative(bg,correction,src)
    source_gap=direct-np.array(pair['total'])
    print(json.dumps(dict(stage='actual_readout',source=pair,direct=direct.tolist(),elapsed=round(time.time()-start,2))),flush=True)
    # Same nonconstant homogeneous receiver test and same physical projector
    # algebra as914/915; the propagator itself is not claimed to be G_Pi.
    weighted=c.aj.prior.weighted;u,z,_=weighted.receiver_modes();w,dw,wstats=weighted.weight_jets(u,z,src.points)
    b,da,jstats=c.aj.extractor(bg,correction,src.points)
    xi,dxi,J,p,v,alpha,rho,residual=b;base=(xi,dxi,p,v,alpha,rho,residual)
    ward,gauge,corr=c.aj.prior.ward_pair(src,base,da,w,dw)
    projected=weighted.weighted_pair(src,weighted.projected_inputs(base,da),w,dw)
    ward_error=c.maximum(ward-projected-gauge)
    assert ward_error<1e-13
    zero={k:np.zeros_like(a) for k,a in v.items()}
    zero_pair=src.response(lambda x,p:zero)
    assert max(abs(a) for a in zero_pair['total'])==0.
    # Exact pairing for finite source arrays: retain original reference terms.
    assert pair['jet_pairing_identity_residual']<1e-13
    return dict(round=921,date='2026-10-06',N=17,projection_quadrature_M=25,time_half_window=T,
      preparation_correction_at_t0_zero=True,original_background_source_and_receiver_retained=True,
      force_is_negative_tangent_kinematic_residual=True,force_slots=['phi','A'],all104_dynamical_variables_propagated=True,
      force_samples=force_stats,original_grid_force=nodalstats,collocation_vs_projection_difference=collocation_gap,
      time_force_interpolation_checks=checks,finite_dynamical_contribution=evolution,
      original_closed_source_pairing=pair,independent_finite_Wilson_derivative=direct.tolist(),
      finite_path_curvature_source_difference_not_error_bound=source_gap.tolist(),independent_record_jet_checks=details,
      original_weighted_Ward_source=weighted.complex_list(ward),direct_projected_weighted_pairing=weighted.complex_list(projected),
      pure_gauge_quadrature_defect=weighted.complex_list(gauge),weighted_Ward_identity_error=ward_error,
      receiver_weight=wstats,extractor_jets=jstats,
      integrated_field_first_jets_used=True,pointwise_jet_replacement_used=False,
      same_source_reference_terms_retained=True,
      contribution_scope='fixed approximate background; tangent kinematic residual only; original source and original weighted source',
      background_error_Hessian_and_readout_variation_included=False,
      initial_geometric_constraint_error_included=False,remaining_momentum_and_receiver_residuals_included=False,
      continuum_projection_error_certified=False,physical_GPi_inverse_implemented=False,
      actual_full_observable_error_certified=False,full872_response_computed=False,full_goal_completed=False,
      elapsed_seconds=round(time.time()-start,3),source_hashes={str(q.relative_to(STAGE)):hashlib.sha256(q.read_bytes()).hexdigest() for q in
       (Path(__file__),STAGE/'920/common_legendre_defect.py',STAGE/'919/weak_gauss_transport.py',STAGE/'908/relational_loop_source.py',STAGE/'904/coupled_boson_tangent.py',STAGE/'915/analytic_material_jets.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    result=run()
    if a.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
