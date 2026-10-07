"""920: one same-action Legendre defect for all joint matter sources.
Quantities are evaluated at the original material-loop points. Completing the
square is an identity, not a proof that the defect's observable impact is small.
"""
from pathlib import Path
from types import SimpleNamespace
import sys,json,argparse,hashlib,time
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(STAGE/'919'))
import canonical_gauss_commutator as prior
weak=prior.weak;c=prior.c;response=prior.response;r=prior.r;mb=prior.mb
TARGET=HERE/'common_legendre_defect_results.json'

class MomentumField:
    def __init__(self,minus,plus,k,T):
        ym,dm=minus;yp,dp=plus;self.T=T;self.k=k.reshape((-1,3));N=k.shape[0]
        coeff=np.stack(((yp+ym)/2-T*(dp-dm)/4,3*(yp-ym)/4-T*(dp+dm)/4,T*(dp-dm)/4,(T*(dp+dm)-(yp-ym))/4),axis=-2)
        self.co=np.fft.fftn(coeff,axes=(0,1,2)).reshape((-1,4,42))/N**3
    def evaluate(self,x):
        s=x[:,0]/self.T;powers=np.stack((s*0+1,s,s*s,s*s*s),axis=1);ans=[]
        for lo in range(0,len(x),64):
            hi=lo+64;phase=x[lo:hi,1:]@self.k.T
            z=np.cos(phase)@self.co.real.reshape((-1,168))-np.sin(phase)@self.co.imag.reshape((-1,168))
            ans.append(np.einsum('bp,bpf->bf',powers[lo:hi],z.reshape((-1,4,42))))
        packed=np.concatenate(ans)
        return dict(pi=packed[:,:6],E=packed[:,6:].reshape((-1,3,12)))

def momentum_pair(T):
    pack=lambda y:np.concatenate((y['pi'],y['E'].reshape(y['pi'].shape[:3]+(36,))),axis=-1)
    with r.ev.ResearchRuntime(r.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=r.ev.Model(17,old);analytic=response.tangent.AnalyticModel(17,old)
        y0,u0,A0,_=response.initial(model,analytic,old);ends=[]
        for sign in (-1,1):
            y={k:v.copy() for k,v in y0.items()};u=u0.copy();A={k:v.copy() for k,v in A0.items()}
            for _ in range(2):y,u,A=response.step(model,analytic,y,u,A,sign*T/2)
            dy,du,dA=response.rhs(model,analytic,y,u,A)
            ends.append(((pack(y),pack(dy)),(pack(A),pack(dA))))
        return MomentumField(ends[0][0],ends[1][0],model.k,T),MomentumField(ends[0][1],ends[1][1],model.k,T)

def quantities(p,Pi,par):
    # Only the first-jet/ADM dictionary is needed; the second-jet fields in p
    # are used by c.geometry to provide the same inherited geometry object.
    z=c.geometry(p);G,*_=c.target(p['phi'],par);Gi=np.linalg.inv(G)
    rp=np.einsum('aij,bj->bai',c.REP,p['phi'])
    D=p['dphi']+np.einsum('bma,bai->bmi',p['A'],rp)
    F=p['dA']-p['dA'].swapaxes(1,2)+c.bracket(p['A'][:,:,None,:],p['A'][:,None,:,:])
    vphi=D[:,0,:]-np.einsum('bi,biA->bA',z['beta'],D[:,1:,:])
    vA=F[:,0,1:,:]-np.einsum('bj,bjia->bia',z['beta'],F[:,1:,1:,:])
    pi_r=z['vol'][:,None]/z['alpha'][:,None]*np.einsum('bAB,bB->bA',G,vphi)
    E_r=z['vol'][:,None,None]/z['alpha'][:,None,None]*par.K*np.einsum('bij,bja->bia',z['invgamma'],vA)
    dp=pi_r-Pi['pi'];dE=E_r-Pi['E'];weight=z['alpha']/z['vol']
    qform=lambda P,E:np.einsum('bA,bAB,bB->b',P,Gi,P)+np.einsum('bia,bij,bja,a->b',E,z['gamma'],E,1/par.K)
    Hkin=.5*weight*qform(Pi['pi'],Pi['E'])
    L1=np.sum(Pi['pi']*vphi,axis=1)+np.sum(Pi['E']*vA,axis=(1,2))-Hkin
    L2=.5/weight*(np.einsum('bA,bAB,bB->b',vphi,G,vphi)+np.einsum('bia,bij,bja,a->b',vA,z['invgamma'],vA,par.K))
    square=.5*weight*qform(dp,dE)
    rvphi=weight[:,None]*np.einsum('bAB,bB->bA',Gi,dp)
    rvA=weight[:,None,None]*np.einsum('bij,bja->bia',z['gamma'],dE)/par.K
    return dict(L1=L1,L2=L2,square=square,dp=dp,dE=dE,velocity_phi_defect=rvphi,velocity_A_defect=rvA,
      momentum_pi_equation=vphi-weight[:,None]*np.einsum('bAB,bB->bA',Gi,Pi['pi']),
      momentum_E_equation=vA-weight[:,None,None]*np.einsum('bij,bja->bia',z['gamma'],Pi['E'])/par.K)

def run():
    start=time.time();par=c.parameters();bg,field,_=c.ex.previous.rebuilt.build_pair();mc,mv=momentum_pair(bg.T)
    src=c.ex.loop.LoopSource(bg,2,8);x=src.points;p=c.aj.Cache(bg,x).jets();v=c.aj.Cache(field,x).jets()
    Pi=mc.evaluate(x);dPi=mv.evaluate(x);q=quantities(p,Pi,par);h=1e-24
    shifted=quantities({k:p[k].astype(complex)+1j*h*v[k] for k in p},{k:Pi[k].astype(complex)+1j*h*dPi[k] for k in Pi},par)
    dq={k:z.imag/h for k,z in shifted.items()}
    square_identity=c.maximum(q['L2']-q['L1']-q['square'])
    first_variation_identity=c.maximum(dq['L2']-dq['L1']-dq['square'])
    mom_eq_error=max(c.maximum(q['velocity_phi_defect']-q['momentum_pi_equation']),c.maximum(q['velocity_A_defect']-q['momentum_E_equation']))
    assert square_identity<1e-12 and first_variation_identity<1e-12 and mom_eq_error<1e-12
    # Actual source dictionary: off-shell A0 test zeta with an independent
    # first spatial jet. This is a variation direction, not a new preparation.
    rng=np.random.default_rng(920);zeta=rng.normal(size=(len(x),12));dzeta=rng.normal(size=(len(x),3,12))
    ap={k:vv.astype(complex).copy() for k,vv in p.items()}
    ap['A'][:,0,:]+=1j*h*zeta;ap['dA'][:,1:,0,:]+=1j*h*dzeta
    actual=quantities(ap,Pi,par)['square'].imag/h
    Dz=dzeta+c.bracket(p['A'][:,1:,:],zeta[:,None,:])
    expected=np.einsum('bA,aAB,bB,ba->b',q['dp'],c.REP,p['phi'],zeta)-np.sum(q['dE']*Dz,axis=(1,2))
    a0_error=c.maximum(actual-expected)
    assert a0_error<1e-13
    # Independently evaluate the 919 momentum map at these same source points.
    m=prior.momenta(p,par)
    dictionary=max(c.maximum(m['pi']-Pi['pi']-q['dp']),c.maximum(m['E']-Pi['E']-q['dE']))
    assert dictionary<1e-13
    statistics=lambda z:dict(max=c.maximum(z),RMS=float(np.sqrt(np.mean(z*z))))
    return dict(round=920,date='2026-10-06',N=17,source_points=len(x),original_material_source_points=True,
      same912_preparation_and_104_variable_flow=True,same_time_window=bg.T,
      algebraic_checks=dict(square_identity=square_identity,actual_joint_first_variation_identity=first_variation_identity,
       auxiliary_momentum_equation_identity=mom_eq_error,A0_source_direction_identity=a0_error,independent919_momentum_dictionary=dictionary),
      actual_momentum_defects={k:statistics(q[k]) for k in ('dp','dE')},
      actual_velocity_equation_defects={k:statistics(q[k]) for k in ('velocity_phi_defect','velocity_A_defect')},
      actual_tangent_momentum_defects={k:statistics(dq[k]) for k in ('dp','dE')},
      actual_tangent_velocity_defects={k:statistics(dq[k]) for k in ('velocity_phi_defect','velocity_A_defect')},
      action_defect=dict(minimum=float(q['square'].min()),maximum=float(q['square'].max()),maximum_first_variation=c.maximum(dq['square']),
        A0_test_first_variation_max=c.maximum(actual)),
      action_defect_smallness_is_not_source_or_observable_bound=True,
      mixed_action_is_same_theory_offshell_dictionary_not_new_physical_counterterm=True,
      full_relation_observable_error_certified=False,quantum_auxiliary_measure_equivalence_proved=False,
      full872_response_computed=False,full_goal_completed=False,elapsed_seconds=round(time.time()-start,3),
      source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
       (Path(__file__),STAGE/'919/canonical_gauss_commutator.py',STAGE/'919/weak_gauss_transport.py',STAGE/'912/receiver_forced_response.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    result=run()
    if a.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
