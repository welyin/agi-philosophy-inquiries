"""906: independent Dirac source/rotation checks and actual compact-mode evolution."""
from pathlib import Path
import argparse,json,time
import numpy as np
import compact_receiver_propagation as r
from curved_fermion_checks import point_geometry
HERE=Path(__file__).resolve().parent;TARGET=HERE/'compact_receiver_checks_results.json'

def point_action(y,z,f,u,du):
    out=r.potential(y,z,f,u)
    for i in range(3):
        for A in range(4):
            out-=1j*f['W'][i+1,A]*r.mat(r.ALPHA[A],du[i])
            out-=.5j*f['dW'][i+1,i+1,A]*r.mat(r.ALPHA[A],u)
    return out

def jets():
    rng=np.random.default_rng(906);errors=dict(flat_stress=0.,flat_Dirac=0.,Clifford=0.,point_Dirac=0.,point_trace=0.,frame_covariance=0.,stress_frame_covariance=0.)
    for a in range(4):
        for b in range(4):errors['Clifford']=max(errors['Clifford'],float(np.max(abs(r.GAMMA[a]@r.GAMMA[b]+r.GAMMA[b]@r.GAMMA[a]-2*r.shared.ETA[a,b]*np.eye(4)))))
    for _ in range(12):
        k=rng.normal(size=3);E=np.sqrt(r.MASS**2+k@k);H=r.MASS*r.BETA+sum(k[i]*r.ALPHA[i+1] for i in range(3));_,O=np.linalg.eigh(H);u=O[:,-1]
        y=dict(g=np.diag([-1.,1.,1.,1.]));z=point_geometry(y['g'],np.zeros((4,4,4)));z.update(vol=np.array(1.),invgamma=np.eye(3));f=r.shared.frame_data(y,z)
        du=np.r_[(-1j*E*u)[None,:],1j*k[:,None]*u];b=r.bilinear_jet(y,z,f,u,du);p=np.r_[-E,k]
        errors['flat_stress']=max(errors['flat_stress'],float(np.max(abs(b['stress']-np.outer(p,p)/E))))
        errors['flat_Dirac']=max(errors['flat_Dirac'],float(np.max(abs(b['dirac_residual']))))
    B=-r.GAMMA[1]@r.GAMMA[2]/2;vals,O=np.linalg.eigh(1j*B);angle=.31;S=(O*np.exp(-1j*angle*vals))@O.conj().T
    L=np.eye(4);L[1:3,1:3]=[[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]
    Z=np.zeros((4,4));Z[1,2]=-1;Z[2,1]=1
    for _ in range(8):
        q=.025*rng.normal(size=(4,4));g=np.diag([-1.,1.,1.,1.])+q+q.T
        d=.04*rng.normal(size=(4,4,4));d=(d+d.transpose(0,2,1))/2;y=dict(g=g);z=point_geometry(g,d)
        z.update(vol=np.sqrt(np.linalg.det(g[1:,1:])),invgamma=np.linalg.inv(g[1:,1:]))
        f=r.shared.frame_data(y,z);u=rng.normal(size=4)+1j*rng.normal(size=4);dx=rng.normal(size=(3,4))+1j*rng.normal(size=(3,4))
        H=point_action(y,z,f,u,dx);du=np.r_[(-1j*H)[None,:],dx];data=r.bilinear_jet(y,z,f,u,du)
        errors['point_Dirac']=max(errors['point_Dirac'],float(np.max(abs(data['dirac_residual']))))
        errors['point_trace']=max(errors['point_trace'],float(abs(np.sum(z['ig']*data['stress'])+r.MASS*data['scalar'])))
        theta=rng.normal(size=4)*.2;f2=r.shared.frame_data(y,z,(L,theta[:,None,None]*(Z@L)))
        u2=S@u;du2=np.array([S@du[i]+theta[i]*B@S@u for i in range(4)])
        H2=point_action(y,z,f2,u2,du2[1:]);errors['frame_covariance']=max(errors['frame_covariance'],float(np.max(abs(H2-S@H-1j*theta[0]*B@S@u))))
        data2=r.bilinear_jet(y,z,f2,u2,du2);errors['stress_frame_covariance']=max(errors['stress_frame_covariance'],float(np.max(abs(data2['stress']-data['stress']))))
    assert max(errors.values())<2e-12,errors
    return dict(random_flat_momenta=12,random_curved_jets=8,errors=errors,normalization_128_256_difference=abs(r.normalization(128)-r.NORM),normalization_512_256_difference=abs(r.normalization(512)-r.NORM))

def offshell_action():
    # Fixed spinor jet, constant metric variation, symmetric vierbein identification.
    # Do not use Dirac equations to simplify either independent density derivative.
    rng=np.random.default_rng(908);errors=[];omission=[]
    for _ in range(10):
        g=np.diag([-1.,1.,1.,1.]);z=point_geometry(g,np.zeros((4,4,4)));z.update(vol=np.array(1.),invgamma=np.eye(3));y=dict(g=g);f=r.shared.frame_data(y,z)
        psi=rng.normal(size=4)+1j*rng.normal(size=4);dpsi=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        h=.15*rng.normal(size=(4,4));h=(h+h.T)/2
        data=r.bilinear_jet(y,z,f,psi,dpsi)
        def density(t):
            E=f['E']@r.shared.vertex.root_near_one(np.eye(4)+t*z['ig']@h)
            C=np.einsum('mA,Aab->mab',np.linalg.inv(E),r.ALPHA)
            L=-np.imag(np.einsum('a,mab,mb->',psi.conj(),C,dpsi))-r.MASS*np.vdot(psi,r.BETA@psi).real
            return abs(np.linalg.det(E))*L
        exact=.5*np.sum((z['ig']@data['stress']@z['ig'])*h)
        wrong=.5*np.sum((z['ig']@data['onshell_stress']@z['ig'])*h)
        row=[abs((density(hh)-density(-hh))/(2*hh)-exact) for hh in (2e-3,1e-3,5e-4)]
        errors.append(row);omission.append(abs(exact-wrong))
    maxima=np.max(errors,axis=0).tolist()
    assert maxima[-1]<1e-6 and maxima[-1]<maxima[0]/10 and max(omission)>.01,(maxima,omission)
    # One free incoming leg and one inhomogeneous response leg, as in872.
    y=dict(g=np.diag([-1.,1.,1.,1.]));z=point_geometry(y['g'],np.zeros((4,4,4)));z.update(vol=np.array(1.),invgamma=np.eye(3));f=r.shared.frame_data(y,z)
    u=np.array([1.,0,0,0],complex);v=np.zeros(4,complex);du=np.zeros((4,4),complex);dv=np.zeros((4,4),complex)
    du[0]=-1j*u;dv[0]=-1j*u
    plus=r.bilinear_jet(y,z,f,u+v,du+dv);minus=r.bilinear_jet(y,z,f,u-v,du-dv)
    mix=(plus['stress']-minus['stress'])/4;short=(plus['onshell_stress']-minus['onshell_stress'])/4
    assert np.max(abs(mix-np.diag([0.,.5,.5,.5])))<1e-14
    assert np.max(abs(short-np.diag([.5,0.,0.,0.])))<1e-14
    return dict(random_offshell_jets=10,density_variation_errors=maxima,omitting_lagrangian_term_max_error=max(omission),
        one_free_one_forced_leg_full_stress=mix.tolist(),one_free_one_forced_leg_on_shell_shortcut=short.tolist())

def discrete():
    with r.ev.ResearchRuntime(r.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=r.ev.Model(9,old);y,_=model.initial();y=r.ev.rk4(model,y,.005);z=model.geometry(y);f=r.shared.frame_data(y,z)
        rng=np.random.default_rng(907);u=rng.normal(size=(9,9,9,4))+1j*rng.normal(size=(9,9,9,4));v=rng.normal(size=u.shape)+1j*rng.normal(size=u.shape)
        Hu=r.apply(model,y,z,f,u);Hv=r.apply(model,y,z,f,v)
        herm=abs(np.mean(np.sum(u.conj()*Hv-Hu.conj()*v,axis=-1)))
        assert herm<1e-12
        return dict(N=9,time=.005,discrete_Hermitian_pairing_error=float(herm))

def band(u):
    N=len(u);coeff=np.fft.fftn(u,axes=(0,1,2))/N**3;idx=np.r_[np.arange(3),np.arange(N-2,N)]
    return coeff[np.ix_(idx,idx,idx,np.arange(4))]

def run():
    out=dict(round=906,date='2026-10-06',jets=jets(),offshell_action=offshell_action(),discrete=discrete());print('Independent Dirac, source and rotation checks passed.',flush=True)
    rows=[];saved={}
    for N,steps in ((17,4),(25,4),(25,8),(33,4)):
        start=time.time();row,u=r.flow(N,steps);rows.append(row);saved[N,steps]=u
        print('N',N,'steps',steps,'stress_div',row['stress_density_divergence_max'],'seconds',round(time.time()-start,2),flush=True)
    fine=[rows[0],rows[1],rows[3]]
    for key in ('stress_density_divergence_max','covariant_Dirac_residual_max','current_density_divergence_max'):
        assert fine[-1][key]<fine[0][key],(key,[q[key] for q in fine])
    assert all(q['sampled_inner_scalar_min']>0 and q['norm_drift']<1e-9 for q in rows)
    temporal=float(np.max(abs(saved[25,4]-saved[25,8])))
    spatial=[float(np.max(abs(band(saved[a,4])-band(saved[b,4])))) for a,b in ((17,25),(25,33))]
    assert temporal<1e-8 and spatial[1]<spatial[0],(temporal,spatial)
    out.update(rows=rows,preparation=dict(mass=r.MASS,center=r.CENTER.tolist(),plateau_radius=r.INNER,support_radius=r.OUTER,normalization=r.NORM,m=.5,r0=.5,r1=1.,species_count=3,identical_mode_allowed_distinct_species=True),
        N25_time_halving_max=temporal,common_band_differences_17_25_and_25_33=spatial,
        refinement_need_not_be_monotone=True,initial_nonmonotone_current_test_preserved=True,all_checks_passed=True,formal_reports=906,fresh_numbered_groups=1,cumulative_numbered_groups=3691,
        argument_scope='Original neutral doublet receiver mode with a single grid-independent compact preparation, propagated on902; scalar response and full off-shell covariant bilinear stress mapped to905 full gravity force and normal constraints. Sampled PDE, Ward and refinement diagnostics only; no full872 source or response yet.',
        absolute_reference_stress_evaluated=False,original_read_profile_b_evaluated=False,physical_source_projection_implemented=False,
        full_four_leg_source_evaluated=False,finite_time_error_certified=False,full_goal_completed=False)
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--checks-only',action='store_true');a=p.parse_args()
    out=dict(jets=jets(),offshell_action=offshell_action(),discrete=discrete()) if a.checks_only else run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    elif not a.checks_only:assert out==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(out,ensure_ascii=False,indent=2))
