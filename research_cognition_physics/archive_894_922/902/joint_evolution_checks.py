"""902: independent Hamiltonian/stress checks and short-time background diagnostics."""
from pathlib import Path
import argparse,json,time,sys
import numpy as np
import common_background_evolution as ev
HERE=Path(__file__).resolve().parent;TARGET=HERE/'joint_evolution_checks_results.json'
def metric_identity_check():
    rng=np.random.default_rng(902);errors=[]
    for _ in range(20):
        a=rng.normal(size=(4,4));g=np.diag([-1.,1.,1.,1.])+.03*(a+a.T);ig=np.linalg.inv(g)
        d=.1*rng.normal(size=(4,4,4));d=(d+np.swapaxes(d,-1,-2))/2
        dd=.2*rng.normal(size=(4,4,4,4));dd=(dd+np.swapaxes(dd,0,1))/2;dd=(dd+np.swapaxes(dd,2,3))/2
        di=-np.einsum('ab,mbc,cd->mad',ig,d,ig)
        Z=d+np.swapaxes(d,0,1)-np.moveaxis(d,0,-1);G=.5*np.einsum('rs,mns->rmn',ig,Z)
        dz=dd+np.swapaxes(dd,1,2)-np.moveaxis(dd,1,-1)
        dg=.5*np.einsum('ars,mns->armn',di,Z)+.5*np.einsum('rs,amns->armn',ig,dz)
        ric=np.einsum('rrmn->mn',dg)-np.einsum('nrmr->mn',dg)+np.einsum('rrl,lmn->mn',G,G)-np.einsum('rnl,lmr->mn',G,G)
        C=np.einsum('ab,abm->m',ig,d)-.5*np.einsum('ab,mab->m',ig,d)
        dC=np.einsum('mab,abn->mn',di,d)-.5*np.einsum('mab,nab->mn',di,d)+np.einsum('ab,mabn->mn',ig,dd)-.5*np.einsum('ab,mnab->mn',ig,dd)
        reduced=ric-.5*(dC+dC.T)+np.einsum('lmn,l->mn',G,C)
        dg0=.5*np.einsum('ars,mns->armn',di,Z)
        ric0=np.einsum('rrmn->mn',dg0)-np.einsum('nrmr->mn',dg0)+np.einsum('rrl,lmn->mn',G,G)-np.einsum('rnl,lmr->mn',G,G)
        dc0=np.einsum('mab,abn->mn',di,d)-.5*np.einsum('mab,nab->mn',di,d)
        Q=ric0-.5*(dc0+dc0.T)+np.einsum('lmn,l->mn',G,C)
        theory=-.5*np.einsum('ab,abmn->mn',ig,dd)+Q
        errors.append(float(np.max(abs(reduced-theory))))
    assert max(errors)<2e-15
    # Non-Abelian representation and adjoint bracket have one convention.
    lie=0.
    for a in range(12):
        for b in range(12):
            rhs=sum((f*ev.REP[c] for c,i,j,f in ev.TERMS if i==a and j==b),np.zeros((6,6)))
            lie=max(lie,float(np.max(abs(ev.REP[a]@ev.REP[b]-ev.REP[b]@ev.REP[a]-rhs))))
    assert lie<1e-14
    return dict(random_metric_jets=20,reduced_Ricci_wave_identity_max=max(errors),gauge_representation_bracket_error=lie)
def variational_check():
    rng=np.random.default_rng(902)
    with ev.ResearchRuntime(ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=ev.Model(12,old);y,_=model.initial();rhs,z,m=model.rhs(y,aux=True);dx=(2*np.pi/12)**3
        # Include nonzero shift: test the full canonical generator and Gauss term.
        y['g'][...,0,1]+=.017;y['g'][...,1,0]+=.017
        rhs,z,m=model.rhs(y,aux=True)
        def Ham(w):
            geo=model.geometry(w);mat=model.matter(w,geo)
            return dx*float(np.sum(geo['alpha']*geo['vol']*mat['rho']+np.einsum('...i,...i->...',geo['beta'],mat['M'])))
        rows=[]
        for key,expected_key,sign in (('phi','pi',-1),('pi','phi',1),('A','E',-1),('E','A',1)):
            d=rng.normal(size=y[key].shape);d/=np.sqrt(dx*np.sum(d*d));analytic=sign*dx*float(np.sum(d*rhs[expected_key]))
            errors=[]
            for eps in (1e-3,5e-4,2.5e-4):
                plus={k:v.copy() for k,v in y.items()};minus={k:v.copy() for k,v in y.items()};plus[key]+=eps*d;minus[key]-=eps*d
                fd=(Ham(plus)-Ham(minus))/(2*eps);errors.append(abs(fd-analytic))
            assert errors[-1]<2e-7,(key,analytic,errors)
            rows.append(dict(variable=key,canonical_derivative=analytic,steps=[1e-3,5e-4,2.5e-4],errors=errors))
        # Stress is independently the metric derivative at fixed canonical data.
        y,_=model.initial();rhs,z,m=model.rhs(y,aux=True)
        d=rng.normal(size=y['g'][...,1:,1:].shape);d=(d+np.swapaxes(d,-1,-2))/2;d/=np.sqrt(dx*np.sum(d*d))
        Sup=np.einsum('...ik,...jl,...kl->...ij',z['invgamma'],z['invgamma'],m['stress'])
        stress_derivative=-.5*dx*np.sum(z['vol']*np.einsum('...ij,...ij->...',Sup,d))
        errors=[]
        for eps in (1e-3,5e-4,2.5e-4):
            plus={k:v.copy() for k,v in y.items()};minus={k:v.copy() for k,v in y.items()};plus['g'][...,1:,1:]+=eps*d;minus['g'][...,1:,1:]-=eps*d
            errors.append(abs((Ham(plus)-Ham(minus))/(2*eps)-stress_derivative))
        assert errors[-1]<2e-7,errors
        # Actual matter energy transfer to evolving geometry on the initial slice.
        work=dx*float(np.sum(z['vol']*np.einsum('...ij,...ij->...',Sup,z['K'])))
        def matter_energy(w):
            gz=model.geometry(w);mm=model.matter(w,gz);return dx*float(np.sum(gz['vol']*mm['rho']))
        energy_errors=[]
        for eps in (1e-4,5e-5,2.5e-5):
            fd=(matter_energy(ev.add(y,rhs,eps))-matter_energy(ev.add(y,rhs,-eps)))/(2*eps)
            energy_errors.append(abs(fd-work))
        assert energy_errors[-1]<2e-5,energy_errors
    return dict(canonical_rows=rows,metric_stress_errors=errors,energy_work=work,energy_exchange_errors=energy_errors)
def point_observables(model,y):
    z=model.geometry(y);m=model.matter(y,z);i=(0,model.N//4,model.N//8)
    phi=y['phi'];h=np.linalg.norm(phi[...,:4],axis=-1);dh=np.einsum('...A,...iA->...i',phi[...,:4],m['D'][...,:,:4])/h[...,None]
    vh=np.sum(phi[...,:4]*m['v'][...,:4],axis=-1)/h
    dt=z['alpha']*vh+np.sum(z['beta']*dh,axis=-1);cov=np.r_[dt[i],dh[i]];ig=z['ig'][i]
    time_norm=-float(cov@ig@cov);assert time_norm>0
    normal=-ig@cov/np.sqrt(time_norm);P=ig+np.outer(normal,normal)
    field=np.zeros((4,4,8));field[1:,1:]=m['curv'][i][...,:8]
    F0=z['alpha'][i]*m['e'][i]+np.einsum('j,jia->ia',z['beta'][i],m['curv'][i])
    field[0,1:]=F0[...,:8];field[1:,0]=-F0[...,:8]
    Mh=.5*np.einsum('mr,ns,mna,rsa->',P,P,field,field)
    return dict(h=float(h[i]),singlet=float(phi[i][4]),probe=float(phi[i][5]),magnetic_h=float(Mh),clock_timelike_margin=time_norm,
        rho=float(m['rho'][i]),lapse=float(z['alpha'][i]),shift=z['beta'][i].tolist(),metric=y['g'][i].tolist())
def run():
    start=time.time();out=dict(round=902,status='working',formal_previous=901,cumulative_previous=3686,
        identities=metric_identity_check(),variational=variational_check())
    print('Independent action/stress checks passed.',flush=True)
    rows=[];states={}
    for N,steps in ((12,8),(16,8),(24,8),(24,16),(32,8)):
        row,model,first,last=ev.run(N,steps,.01)
        row['point_initial']=point_observables(model,first);row['point_final']=point_observables(model,last)
        rows.append(row)
        if N==24:states[steps]=last
        print('evolution',N,steps,row['final'],flush=True)
    differences={k:float(np.max(abs(states[8][k]-states[16][k]))) for k in states[8]}
    out.update(rows=rows,fixed_N24_time_refinement_differences=differences,elapsed_seconds=time.time()-start,
        original_same_classical_background=True,all_gauge_and_scalar_components_evolved=True,
        quantum_modes_evolved=False,interval_time_propagation_certificate=False,full_goal_completed=False)
    return out
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');ap.add_argument('--identities-only',action='store_true');a=ap.parse_args()
    if a.identities_only:r=dict(identities=metric_identity_check(),variational=variational_check())
    else:r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(r,ensure_ascii=False,indent=2))
