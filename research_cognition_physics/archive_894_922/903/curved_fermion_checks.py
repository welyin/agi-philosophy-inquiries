"""903 independent tetrad covariance, discrete CAR structure, and actual spatial propagation."""
from pathlib import Path
import sys,json,argparse,time
import numpy as np
import curved_fermion_propagation as f
HERE=Path(__file__).resolve().parent;TARGET=HERE/'curved_fermion_checks_results.json'

def point_geometry(g,d):
    ig=np.linalg.inv(g);gamma=g[1:,1:];alpha=1/np.sqrt(-ig[0,0]);beta=-ig[0,1:]/ig[0,0]
    Z=d+np.swapaxes(d,0,1)-np.moveaxis(d,0,-1)
    conn=.5*np.einsum('rs,mns->rmn',ig,Z)
    return dict(ig=ig,gamma=gamma,alpha=alpha,beta=beta,d=d,Gamma=conn)

def frame_checks():
    rng=np.random.default_rng(903);maximum=dict(metric=0.,lorentz=0.,root_derivative=0.,rotation=0.,hermiticity=0.)
    B=-f.vertex.GAMMA[1]@f.vertex.GAMMA[2]/2
    vals,O=np.linalg.eigh(1j*B);angle=.31;S=(O*np.exp(-1j*angle*vals))@O.conj().T
    L=np.eye(4);L[1:3,1:3]=[[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]
    Z=np.zeros((4,4));Z[1,2]=-1;Z[2,1]=1
    for _ in range(8):
        a=.025*rng.normal(size=(4,4));g=np.diag([-1.,1.,1.,1.])+a+a.T
        d=.04*rng.normal(size=(4,4,4));d=(d+np.swapaxes(d,1,2))/2
        y=dict(g=g,phi=rng.normal(size=6)*.15,A=rng.normal(size=(3,12))*.05)
        z=point_geometry(g,d);fr=f.frame_data(y,z)
        maximum['metric']=max(maximum['metric'],float(np.max(abs(fr['E'].T@f.ETA@fr['E']-g))))
        maximum['lorentz']=max(maximum['lorentz'],float(np.max(abs(fr['low']+fr['low'].transpose(0,2,1)))))
        eps=1e-5
        for mu in range(4):
            yp=y|dict(g=g+eps*d[mu]);ym=y|dict(g=g-eps*d[mu])
            ep=f.frame_data(yp,point_geometry(yp['g'],d))['E'];em=f.frame_data(ym,point_geometry(ym['g'],d))['E']
            maximum['root_derivative']=max(maximum['root_derivative'],float(np.max(abs((ep-em)/(2*eps)-fr['dE'][mu]))))
        theta=rng.normal(size=4)*.2;rot=f.frame_data(y,z,(L,theta[:,None,None]*(Z@L)))
        u=rng.normal(size=(64,2))+1j*rng.normal(size=(64,2));du=rng.normal(size=(3,64,2))+1j*rng.normal(size=(3,64,2))
        left=f.point_operator(y,z,rot,S@u,np.array([S@du[i]+theta[i+1]*B@S@u for i in range(3)]))
        right=S@f.point_operator(y,z,fr,u,du)+1j*theta[0]*B@S@u
        maximum['rotation']=max(maximum['rotation'],float(np.max(abs(left-right))))
        V=sum((c*m for c,m in f.terms(y,z,fr)),np.zeros((64,64),complex))
        maximum['hermiticity']=max(maximum['hermiticity'],float(np.max(abs(V-V.conj().T))))
    assert maximum['rotation']<1e-12,maximum
    assert maximum['root_derivative']<1e-8 and max(maximum[k] for k in ('metric','lorentz','hermiticity'))<1e-12
    lie=0.
    for a in range(12):
        for b in range(12):
            target=sum((c*1j*f.Q[j] for j,i,k,c in f.ev.TERMS if i==a and k==b),np.zeros((64,64),complex))
            lie=max(lie,float(np.max(abs((1j*f.Q[a])@(1j*f.Q[b])-(1j*f.Q[b])@(1j*f.Q[a])-target))))
    assert lie<1e-12
    return dict(random_jets=8,errors=maximum,full_internal_Lie_representation_error=lie,nonzero_spin_basis_terms=len(f.SPIN))

def discrete_checks():
    with f.ev.ResearchRuntime(f.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=f.ev.Model(9,old);y,_=model.initial()
        # Actual non-conformal future; do not test spin terms only where they vanish.
        y=f.ev.rk4(model,y,.005);z=model.geometry(y);fr=f.frame_data(y,z)
        rng=np.random.default_rng(904);u=(rng.normal(size=(9,9,9,64,2))+1j*rng.normal(size=(9,9,9,64,2)))/8
        Hu=f.apply(model,y,z,fr,u);inner=f.gram(u,Hu)
        herm=float(np.max(abs(inner-inner.conj().T)))
        Cu=u[...,f.CHARGE,:].conj();reality=float(np.max(abs(f.apply(model,y,z,fr,Cu)+Hu[...,f.CHARGE,:].conj())))
        no_spin=f.apply(model,y,z,fr,u,omit_spin=True);no_color=f.apply(model,y,z,fr,u,omit_color=True)
        spin=float(np.sqrt(np.mean(np.sum(abs(Hu-no_spin)**2,axis=-2))))
        color=float(np.sqrt(np.mean(np.sum(abs(Hu-no_color)**2,axis=-2))))
        assert herm<1e-12 and reality<1e-12 and spin>1e-8 and color>.01,(herm,reality,spin,color)
        return dict(N=9,time=.005,discrete_hermiticity_error=herm,Nambu_reality_error=reality,
            omitted_spin_action_norm=spin,omitted_color_action_norm=color,
            time_discretization_not_exactly_unitary=True)

def run():
    out=dict(round=903,date='2026-10-06',frame=frame_checks(),discrete=discrete_checks());print('Frame and CAR checks passed.',flush=True)
    rows=[];selected={}
    for N,steps in ((9,8),(13,8),(13,16),(17,8)):
        start=time.time();row,model,y,u=f.flow(N,steps)
        rows.append(row);selected[(N,steps)]=u
        print('propagation',N,steps,'Gram',row['gram_error'],'seconds',round(time.time()-start,2),flush=True)
        assert row['gram_error']<2e-10
    temporal=float(np.max(abs(selected[(13,8)]-selected[(13,16)])))
    # Compare Fourier vectors on a common low-mode set; no interpolation at unequal points.
    def band(u):
        n=len(u);c=np.fft.fftn(u,axes=(0,1,2))/n**3;idx=np.r_[np.arange(3),np.arange(n-2,n)]
        return c[np.ix_(idx,idx,idx,np.arange(64),np.arange(2))]
    spatial=[float(np.max(abs(band(selected[(a,8)])-band(selected[(b,8)])))) for a,b in ((9,13),(13,17))]
    assert temporal<1e-10 and spatial[1]<spatial[0],(temporal,spatial)
    out.update(rows=rows,time_halving_N13_max=temporal,common_band_errors_9_13_and_13_17=spatial,
        all_checks_passed=True,formal_reports=903,fresh_numbered_groups=1,cumulative_numbered_groups=3688,
        argument_scope='Full original 64 Nambu spatial Dirac/Majorana operator and propagator on the902 classical joint trajectory; tetrad covariance and CAR checks plus grid/time diagnostics. Smearing test modes are not a reselected vacuum or a computed prepared covariance.',
        same_original_mass_gauge_and_metric=True,actual_spatial_PDE_evolved=True,
        auxiliary_past_covariance_evaluated=False,quantum_source_evaluated=False,quantum_backreaction_evolved=False,
        rigorous_background_or_mode_error_enclosure=False,full_goal_completed=False)
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--checks-only',action='store_true');a=p.parse_args()
    r=dict(frame=frame_checks(),discrete=discrete_checks()) if a.checks_only else run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    elif not a.checks_only:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
