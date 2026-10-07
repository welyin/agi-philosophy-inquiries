"""905: full action source, canonical forcing and joint Ward diagnostics on original902."""
from pathlib import Path
import argparse,json
import numpy as np
import canonical_source_bridge as b
HERE=Path(__file__).resolve().parent;TARGET=HERE/'source_bridge_checks_results.json'

def ward(model,y,h=5e-5):
    dy,z,m=model.rhs(y,aux=True);jg,jp,ja,ym=b.gauge_kinetic_source(model,y);force,const,s=b.convert(y,z,jg,jp,ja)
    mixed=np.einsum('...ma,...an->...mn',z['ig'],s);mu=z['alpha']*z['vol'];dens=mu[...,None,None]*mixed
    spatial=model.grad(dens)
    diff=np.einsum('...iin->...n',spatial[..., :,1:,:])
    gauge=np.einsum('...iia->...a',model.grad(ja[...,1:,:]))+sum(b.ev.bracket(y['A'][...,i,:],ja[...,i+1,:]) for i in range(3))
    ds=[];dj=[]
    for sign in (1,-1):
        yy=b.ev.add(y,dy,sign*h);zz=model.geometry(yy);gj,pj,aj,_=b.gauge_kinetic_source(model,yy);_,_,ss=b.convert(yy,zz,gj,pj,aj)
        ds.append((zz['alpha']*zz['vol'])[...,None]*np.einsum('...a,...an->...n',zz['ig'][...,0,:],ss))
        dj.append(aj[...,0,:])
    diff+=(ds[0]-ds[1])/(2*h);gauge+=(dj[0]-dj[1])/(2*h)
    diff-=mu[...,None]*np.einsum('...rmn,...mr->...n',z['Gamma'],mixed)
    F=np.zeros(y['A'].shape[:3]+(4,4,12));F[...,1:,1:,:]=m['curv'];F[...,0,1:,:]=dy['A'];F[...,1:,0,:]=-dy['A']
    DP=np.concatenate((dy['phi'][...,None,:],m['D']),axis=-2)
    rhs=np.einsum('...ma,...nma->...n',ja,F)+np.einsum('...A,...nA->...n',jp,DP)
    gauge+=np.einsum('...A,...aA->...a',jp,m['rp'])
    return dict(diffeomorphism_density_max=float(np.max(abs(diff-rhs))),internal_density_max=float(np.max(abs(gauge))),
        stress_divergence_max=float(np.max(abs(diff))),exchange_terms_max=float(np.max(abs(rhs))),
        constraint_load_max={k:float(np.max(abs(v))) for k,v in const.items()})

def independent_checks():
    with b.ev.ResearchRuntime(b.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=b.ev.Model(10,old);y,_=model.initial();y=b.ev.rk4(model,y,.005);dy,z,m=model.rhs(y,aux=True)
        jg,jp,ja,ym=b.gauge_kinetic_source(model,y);force,const,s=b.convert(y,z,jg,jp,ja)
        rows=[]
        for h in (.01,.005,.0025):
            vp=b.deformed_velocity_in_base_variables(model,y,h);vm=b.deformed_velocity_in_base_variables(model,y,-h)
            er={k:float(np.max(abs((vp[k]-vm[k])/(2*h)-force[k]))) for k in force}
            cp=b.deformed_constraints(model,y,h);cm=b.deformed_constraints(model,y,-h)
            ce={k:float(np.max(abs((cp[kk]-cm[kk])/(2*h)-const[k]))) for k,kk in (('H','H_extra'),('M','M_extra'),('Gauss','Gauss'))}
            rows.append(dict(step=h,force_errors=er,constraint_errors=ce))
        assert rows[-1]['force_errors']['E']<rows[0]['force_errors']['E']/8
        assert max(rows[-1]['force_errors'].values())<5e-9 and max(rows[-1]['constraint_errors'].values())<1e-10
        # Metric and scalar source checked directly against the same action density at fixed F_mu_nu.
        at=(1,2,3);g=y['g'][at];F=np.zeros((4,4,12));F[1:,1:]=m['curv'][at];F[0,1:]=dy['A'][at];F[1:,0]=-dy['A'][at]
        p=y['phi'][at][5];rng=np.random.default_rng(905);dg=rng.normal(size=(4,4))*.1;dg=(dg+dg.T)/2;dp=.13
        def lag(gg,pp):
            ig=np.linalg.inv(gg);mu=np.sqrt(-np.linalg.det(gg));return -.25*mu*pp*np.einsum('mr,ns,mna,rsa,a->',ig,ig,F,F,model.K)
        expected=np.sum(jg[at]*dg)+jp[at][5]*dp;action_errors=[]
        for h in (2e-3,1e-3,5e-4):action_errors.append(abs((lag(g+h*dg,p+h*dp)-lag(g-h*dg,p-h*dp))/(2*h)-expected))
        assert action_errors[-1]<action_errors[0]/8 and action_errors[-1]<1e-8
        # Dropping the induced scalar source or time gauge current is physically distinct.
        return dict(N=10,time=.005,canonical_and_constraint_rows=rows,
            independent_metric_scalar_action_errors=action_errors,
            dropped_scalar_force_max=float(np.max(abs(force['pi']))),
            dropped_Gauss_load_max=float(np.max(abs(const['Gauss']))),
            YM_trace_error=float(np.max(abs(np.einsum('...mn,...mn->...',z['ig'],ym['T'])))),
            finite_kinetic_deformation_is_verification_not_the_quantum_source=True)

def run():
    out=dict(round=905,date='2026-10-06',independent=independent_checks());rows=[]
    with b.ev.ResearchRuntime(b.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        for N in (12,16,24):
            model=b.ev.Model(N,old);y,_=model.initial()
            for _ in range(4):y=b.ev.rk4(model,y,.005/4)
            row=dict(N=N,time=.005,**ward(model,y));rows.append(row);print('Ward',row,flush=True)
    for key in ('diffeomorphism_density_max','internal_density_max'):
        assert all(bb[key]<aa[key] for aa,bb in zip(rows,rows[1:])),key
    out.update(rows=rows,all_checks_passed=True,formal_reports=905,fresh_numbered_groups=1,cumulative_numbered_groups=3690,
        argument_scope='Explicit Euler-density conversion to the original902/904 canonical forced system and all normal constraints, verified by a covariant pF^2 action variation and joint Ward diagnostics. Conditional retarded physical mapping for compatible zero-past sources; no numerical872 source or feedback yet.',
        forcing_and_normal_constraints_share_same_source=True,
        full_physical_retarded_inverse_numerically_solved=False,actual_receiver_or_quantum_source_computed=False,
        rigorous_source_or_time_error_bound=False,full_goal_completed=False)
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--checks-only',action='store_true');a=p.parse_args()
    r=independent_checks() if a.checks_only else run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    elif not a.checks_only:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
