"""912 working: actual compact receiver momentum is a divergence-free spin curl.
Checks the source identity before introducing any numerical counterflow.
"""
from pathlib import Path
import sys,json,hashlib,argparse
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(HERE));sys.path.insert(0,str(STAGE/'906'))
import receiver_initial_completion as init
import compact_receiver_checks as independent
r=init.r;TARGET=HERE/'receiver_momentum_curl_results.json'

def point_jets():
    rng=np.random.default_rng(912);errs=[]
    for _ in range(24):
        conformal=rng.uniform(.85,1.3);grad=rng.normal(size=3)*.13
        kk=rng.normal(size=(3,3))*.05;kk=(kk+kk.T)/2
        y=dict(g=np.diag([-1.,*[conformal**4]*3]));dg=np.zeros((4,4,4))
        dg[0,1:,1:]=-2*kk
        dg[0,0,1:]=dg[0,1:,0]=rng.normal(size=3)*.1
        dg[0,0,0]=rng.normal()*.1
        dg[1:,1:,1:]=4*conformal**3*grad[:,None,None]*np.eye(3)
        z=independent.point_geometry(y['g'],dg);z.update(vol=conformal**6,invgamma=np.eye(3)*conformal**-4)
        fr=r.shared.frame_data(y,z)
        v=rng.normal(size=2)+1j*rng.normal(size=2);v/=np.linalg.norm(v)
        spin=np.array([np.vdot(v,ss@v).real for ss in r.SIG])
        height=rng.uniform(.1,.7);deriv=rng.normal(size=3)*.4
        chi=np.r_[height*v,[0j,0j]];dx=np.zeros((3,4),complex);dx[:,:2]=deriv[:,None]*v
        H=independent.point_action(y,z,fr,chi,dx)
        data=r.bilinear_jet(y,z,fr,chi,np.r_[(-1j*H)[None,:],dx])
        wantJ=-height*np.cross(deriv,spin)/2
        gotJ=conformal**6*data['stress'][0,1:]
        errs.append((init.norm(gotJ-wantJ),abs(conformal**6*data['stress'][0,0]-r.MASS*height**2),init.norm(data['dirac_residual'])))
    values=np.max(errs,axis=0);assert values.max()<1e-13,values
    return dict(actual_Dirac_spin_conventions=True,curved_conformal_jets=24,
        arbitrary_extrinsic_curvature_and_spin_orientation=True,
        momentum_curl_identity_error=float(values[0]),energy_identity_error=float(values[1]),Dirac_error=float(values[2]))

def actual_source(N):
    with r.ev.ResearchRuntime(r.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=r.ev.Model(N,old);y,_=model.initial();data,z,rho,J=init.initial_stress(model,y)
        points,rr=r.coordinates(N);d=(points-r.CENTER+np.pi)%(2*np.pi)-np.pi
        f=r.NORM*r.profile(rr);df=np.zeros_like(rr);sel=(rr>r.INNER)&(rr<r.OUTER)
        t=(rr[sel]-r.INNER)/(r.OUTER-r.INNER);ff=r.profile(rr[sel])
        df[sel]=-r.NORM*ff*(1-ff)*(t**-2+(1-t)**-2)/(r.OUTER-r.INNER)
        grad=df[...,None]*d/np.where(rr>0,rr,1)[...,None]
        expected=-.5*f[...,None]*np.cross(grad,np.array([0.,0.,1.]))
        # A structure-preserving finite Fourier projection: derivative acts on
        # the projected f^2, preserving the exact continuum zero mode.
        density=f*f;projected=-.25*np.cross(model.grad(density),np.array([0.,0.,1.]))
        div=sum(model.grad(projected[...,i])[...,i] for i in range(3))
        dx=(2*np.pi/N)**3
        row=dict(N=N,actual_source_momentum_identity_error=init.norm(J-expected),
            actual_source_energy_identity_error=init.norm(z['vol']*rho-r.MASS*density),
            raw_point_sample_total_momentum=(dx*J.sum(axis=(0,1,2))).tolist(),
            discrete_curl_total_momentum=(dx*projected.sum(axis=(0,1,2))).tolist(),
            discrete_curl_divergence_max=init.norm(div),
            pointwise_difference_between_two_source_approximations=init.norm(projected-J),
            source_approximation_error_certified=False)
        assert row['actual_source_momentum_identity_error']<1e-13
        assert row['actual_source_energy_identity_error']<1e-13
        assert init.norm(row['discrete_curl_total_momentum'])<1e-13
        assert row['discrete_curl_divergence_max']<1e-13
        return row

def run():
    return dict(round=912,status='working',date='2026-10-06',formal_rounds=911,cumulative_numbered_groups=3696,
        point_jets=point_jets(),actual_original_source_rows=[actual_source(N) for N in (17,25)],
        identity='For chi=f(x)(zeta,0), real f and constant normalized zeta on gamma=psi^4 I, alpha=1,beta=0: vol*rho=M f^2, J=-curl(f^2 s)/4, s=zeta^dagger sigma zeta.',
        exact_continuum_total_momentum_zero=True,extra_physical_matter_counterflow_needed_for_this_source=False,
        general_sources_still_need_inherited_counterflow_menus=True,
        full_A_initial_and_time_error_certified=False,full872_feedback_completed=False,full_goal_completed=False,
        source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (Path(__file__),HERE/'receiver_initial_completion.py',STAGE/'906/compact_receiver_checks.py',STAGE/'906/compact_receiver_propagation.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    result=run()
    if a.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
