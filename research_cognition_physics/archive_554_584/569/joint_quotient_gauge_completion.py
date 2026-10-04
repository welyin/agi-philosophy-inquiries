"""569: same Higgs, residual Gauss, and a faithful compact-quotient magnetic H.

Classical matching prescription; no chiral-fermion or quantum continuum claim.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_scalar_propagation_matching as scalar
import joint_gauge_matter_propagation as graph
from protected_pair_rg import inputs, gauge_squared

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_quotient_gauge_completion_results.json'


def parameters():
    old,bound=inputs()
    row=json.loads((HERE/'joint_singlet_common_mass_rg_results.json').read_text('utf8'))['examples'][2]
    gy2,gw2,gc2=gauge_squared(-row['u'],old['matched_inverse_couplings'])
    _,u,_=scalar.parameters()
    b=np.array([gc2/2,gw2/2,gy2/72])
    K=1/(2*b)  # k_h/v=c=epsilon=1, explicit classical prescription.
    wq,wl=1.,bound['r']
    indices=np.array([2*wq,(3*wq+wl)/2,66*wq+54*wl])
    delta=.5*float(np.min(K/indices))
    return dict(g2=np.array([gc2,gw2,gy2]),b=b,K=K,wq=wq,wl=wl,
        indices=indices,delta=delta,h2=float(u[0]),scale=row['mu_GeV'],u=row['u'])


def generators(n):
    out=[]
    for i in range(n):
        for j in range(i+1,n):
            a=np.zeros((n,n),complex);a[i,j]=a[j,i]=.5;out.append(a)
            a=np.zeros((n,n),complex);a[i,j]=-.5j;a[j,i]=.5j;out.append(a)
    for k in range(1,n):
        a=np.zeros((n,n),complex)
        for i in range(k):a[i,i]=1
        a[k,k]=-k;out.append(a/np.sqrt(2*k*(k+1)))
    return np.array(out)


def group_exp(coords,n):
    h=np.einsum('a,aij->ij',coords,generators(n))
    vals,vecs=np.linalg.eigh(h)
    return (vecs*np.exp(1j*vals))@vecs.conj().T


def coefficients(p,fraction=.5):
    delta=fraction*float(np.min(p['K']/p['indices']))
    residual=(p['K']-delta*p['indices'])/np.array([3.,2.,36.])
    return delta,residual


def components(C,W,z,p):
    tc,tw=np.trace(C),np.trace(W)
    chi=p['wq']*(tc*tw*z+tc.conjugate()*(z**-4+z**2))
    chi+=p['wl']*(tw*z**-3+z**6+1)
    faithful=12*p['wq']+4*p['wl']-chi.real
    return float(faithful),np.array([9-abs(tc)**2,4-abs(tw)**2,1-(z**6).real])


def potential(C,W,z,p,fraction=.5):
    delta,residual=coefficients(p,fraction)
    faithful,blind=components(C,W,z,p)
    return float(delta*faithful+residual@blind)


def units_check():
    p=parameters();gc2,gw2,gy2=p['g2'];bc,bw,bu=p['b'];h2=p['h2']
    assert np.max(abs(2*p['b']*p['K']-1))<1e-14
    H=h2/4*np.array([[gw2,-np.sqrt(gw2*gy2)],[-np.sqrt(gw2*gy2),gy2]])
    photon=np.sqrt([gy2,gw2]);photon/=np.linalg.norm(photon)
    assert np.linalg.norm(H@photon)<1e-14
    eig=np.linalg.eigvalsh(H);gap=h2*(gw2+gy2)/4
    assert np.max(abs(eig-[0,gap]))<1e-14
    delta,res=coefficients(p)
    assert min(res)>0
    # Old single-trace relation is not an all-scale equality.
    trace_defect=1/gy2-3/gw2+4/(3*gc2)
    assert abs(trace_defect)>1e-2
    return dict(same_frozen_scale_GeV=p['scale'],same_rg_u=p['u'],
        couplings_squared_colour_weak_Y=p['g2'].tolist(),
        electric_coefficients_colour_weak_integer_circle=p['b'].tolist(),
        target_magnetic_Hessian=p['K'].tolist(),Higgs_integer_charge=3,
        h_squared=h2,W_gap_squared=float(gw2*h2/4),Z_gap_squared=float(gap),
        photon_direction=photon.tolist(),positive_weight_ratio=p['wl'],
        indices=p['indices'].tolist(),delta=delta,residual_positive_coefficients=res.tolist(),
        old_single_trace_relation_defect_at_this_scale=float(trace_defect),
        MS_to_classical_matching_is_an_input=True,not_a_new_mass_prediction=True)


def neutral_symplectic_check():
    p0=parameters();_,bw,bu=p0['b'];_,kw,k0=p0['K'];h2=p0['h2'];mu=h2/4
    D,C,_,_=graph.square();E,N=D.shape;rng=np.random.default_rng(569)
    p=rng.normal(size=E);e=rng.normal(size=E);e-=D@np.linalg.lstsq(D,e,rcond=None)[0]
    a,alpha,phi=[rng.normal(size=n) for n in (E,E,N)]
    pp=-D.T@p;pu=e-6*p;theta=a+6*alpha+D@phi
    full=2/h2*(pp@pp)+bw*(p@p)+bu*(pu@pu)+mu/2*(a@a)
    full+=kw/2*np.sum((C@theta)**2)+k0/2*np.sum((C@alpha)**2)
    reduced=2/h2*np.sum((D.T@p)**2)+bw*(p@p)+bu*np.sum((e-6*p)**2)
    reduced+=mu/2*(a@a)+kw/2*np.sum((C@a+6*C@alpha)**2)+k0/2*np.sum((C@alpha)**2)
    da,dal,dphi=[rng.normal(size=n) for n in (E,E,N)]
    form=pp@dphi+p@(da+6*dal+D@dphi)+pu@dal-p@da-e@dal
    pp_dot=mu*D.T@a;pw_dot=-mu*a-kw*C.T@C@theta
    pu_dot=6*mu*a-k0*C.T@C@alpha
    gauss=max(np.max(abs(pp+D.T@p)),np.max(abs(-6*pp+D.T@pu)))
    gauss_dot=max(np.max(abs(pp_dot+D.T@pw_dot)),np.max(abs(-6*pp_dot+D.T@pu_dot)))
    actual_adot=2*bw*p-12*bu*pu-D@(4*pp/h2)
    predicted_adot=2*bw*p-12*bu*(e-6*p)+4*D@D.T@p/h2
    errs=[abs(full-reduced),abs(form),float(gauss),float(gauss_dot),float(np.max(abs(actual_adot-predicted_adot)))]
    assert max(errs)<1e-11
    return dict(energy_one_form_Gauss_and_flow_residuals=errs,
        residual_constraint='D.T @ e = 0',residual_gauge='alpha ~ alpha + D gamma')


def complete_graph_spectrum_check():
    p=parameters();_,bw,bu=p['b'];_,kw,k0=p['K'];h2=p['h2'];mu=h2/4
    mz=2*(bw+36*bu)*mu;rows=[]
    for d in (2,3):
        D,C,_,_,points=graph.torus(3,d);E,N=D.shape
        _,sing,Vh=np.linalg.svd(D.T,full_matrices=True);rank=int(np.sum(sing>1e-10))
        Q=Vh[rank:].T;n=Q.shape[1];curl=C.T@C
        A=np.block([[2*(bw+36*bu)*np.eye(E)+4*D@D.T/h2,-12*bu*Q],
                    [-12*bu*Q.T,2*bu*np.eye(n)]])
        K=np.block([[mu*np.eye(E)+kw*curl,6*kw*curl@Q],
                    [6*kw*Q.T@curl,(36*kw+k0)*Q.T@curl@Q]])
        ev,U=np.linalg.eigh(A);assert min(ev)>0
        root=(U*np.sqrt(ev))@U.T;got=np.linalg.eigvalsh(root@K@root)
        lambdas=np.linalg.eigvalsh(D.T@D);expected=[0.]*d+[mz]*d
        for lam in lambdas[1:]:expected += [float(lam)]*(d-1)+[float(lam+mz)]*d
        error=float(np.max(abs(got-np.sort(expected))))
        assert len(expected)==E+n and error<2e-10 and rank==N-1
        rows.append(dict(d=d,sites=N,neutral_physical_linear_coordinates=E+n,
            harmonic_photon_zero_frequencies=d,nonzero_momentum_bosonic_modes_excluding_colour=4*d+1,
            spectrum_error=error,minimum_kinetic_eigenvalue=float(min(ev))))
    return dict(rows=rows,no_extra_photon_longitudinal_mode=True,
        scope='linear tangent phase space; global nonlinear total-charge constraint not counted')


def nonlinear_Higgs_hessian_check():
    p=parameters();D,C,edges,faces=graph.square();E,N=D.shape;h=np.sqrt(p['h2'])
    def fun(q):
        phi=q[:N];theta=q[N:N+E];alpha=q[N+E:]
        X=h*np.exp(-.5j*phi);val=0.
        for e,(i,j) in enumerate(edges):val+=abs(X[i]-np.exp(-.5j*theta[e]+3j*alpha[e])*X[j])**2/2
        for f in range(len(C)):
            W=group_exp([0,0,float(C[f]@theta)],2)
            val+=potential(np.eye(3),W,np.exp(1j*(C[f]@alpha)),p)
        return float(val)
    B=np.column_stack((-D,np.eye(E),-6*np.eye(E)))
    expected=p['h2']/4*(B.T@B)
    expected[N:N+E,N:N+E]+=p['K'][1]*(C.T@C)
    expected[N+E:,N+E:]+=p['K'][2]*(C.T@C)
    step=2e-4;found=scalar.hessian(fun,np.zeros(N+2*E),step)
    error=float(np.max(abs(found-expected)))
    assert error<5e-4
    return dict(variables=N+2*E,finite_difference_step=step,full_potential_Hessian_max_error=error)


def centre_obstruction_and_zeros_check():
    p=parameters();omega=np.exp(2j*np.pi/3);z=np.exp(1j*np.pi/3)
    beta_w=4/p['g2'][1];beta_u=36/p['g2'][2]
    cover_at_quotient_identity=2*beta_w+.5*beta_u
    assert cover_at_quotient_identity>1
    zero=[];classes=set();worst_blind=0.
    for a,b,k in itertools.product(range(3),range(2),range(6)):
        C=omega**a*np.eye(3);W=(-1.)**b*np.eye(2);zz=z**k
        faithful,blind=components(C,W,zz,p);worst_blind=max(worst_blind,float(np.max(abs(blind))))
        val=potential(C,W,zz,p)
        is_kernel=(a==k%3 and b==k%2)
        assert (val<1e-9)==is_kernel
        if is_kernel:zero.append([a,b,k])
        classes.add(min(((a+t)%3,(b+t)%2,(k+t)%6) for t in range(6)))
    assert len(zero)==6 and len(classes)==6 and worst_blind<1e-12
    return dict(cover_fundamental_magnetic_energy_for_identity_lift=0.,
        cover_energy_for_equivalent_central_lift=float(cover_at_quotient_identity),
        blind_zero_points=36,blind_zero_classes_in_quotient=len(classes),
        faithful_completion_zero_lifts=zero,faithful_completion_zero_classes=1,
        identity_face_does_not_fix_global_flat_holonomy=True)


def Lie_hessian_check():
    p=parameters();step=3e-4
    def fun(x):return potential(group_exp(x[:8],3),group_exp(x[8:11],2),np.exp(1j*x[11]),p)
    target=np.diag(np.r_[np.repeat(p['K'][0],8),np.repeat(p['K'][1],3),p['K'][2]])
    found=scalar.hessian(fun,np.zeros(12),step)
    err=float(np.max(abs(found-target)))
    off=found-np.diag(np.diag(found))
    assert err<2e-4 and np.max(abs(off))<2e-5
    return dict(generators=12,step=step,expected_diagonal=target.diagonal().tolist(),
        finite_difference_max_error=err,cross_block_max=float(np.max(abs(off))))


def nonlinear_gauge_and_lifts_check():
    p=parameters();rng=np.random.default_rng(5691);D,C,edges,faces=graph.square()
    X=rng.normal(size=(4,2))+1j*rng.normal(size=(4,2))
    Cs=[group_exp(rng.normal(size=8),3) for _ in edges]
    Ws=[group_exp(rng.normal(size=3),2) for _ in edges]
    zs=np.exp(1j*rng.normal(size=4))
    def energy(X,Cs,Ws,zs):
        val=sum(np.linalg.norm(X[i]-zs[e]**3*Ws[e]@X[j])**2/2 for e,(i,j) in enumerate(edges))
        for face in faces:
            pc=np.eye(3,dtype=complex);pw=np.eye(2,dtype=complex);pz=1.+0j
            for e,sign in face:
                pc=pc@(Cs[e] if sign==1 else Cs[e].conj().T)
                pw=pw@(Ws[e] if sign==1 else Ws[e].conj().T)
                pz*=zs[e]**sign
            val+=potential(pc,pw,pz,p)
        return float(val)
    gc=[group_exp(rng.normal(size=8),3) for _ in range(4)]
    gw=[group_exp(rng.normal(size=3),2) for _ in range(4)];gz=np.exp(1j*rng.normal(size=4))
    transformed=energy(np.array([gz[i]**3*gw[i]@X[i] for i in range(4)]),
        [gc[i]@Cs[e]@gc[j].conj().T for e,(i,j) in enumerate(edges)],
        [gw[i]@Ws[e]@gw[j].conj().T for e,(i,j) in enumerate(edges)],
        np.array([gz[i]*zs[e]/gz[j] for e,(i,j) in enumerate(edges)]))
    original=energy(X,Cs,Ws,zs);gauge_err=abs(transformed-original)
    lift_errors=[];minimum=float('inf')
    for _ in range(12):
        lifts=rng.integers(0,6,size=4)
        value=energy(X,[a*np.exp(2j*np.pi*m/3) for a,m in zip(Cs,lifts)],
            [a*(-1)**int(m) for a,m in zip(Ws,lifts)],zs*np.exp(1j*np.pi*lifts/3))
        lift_errors.append(abs(value-original))
        cc=group_exp(rng.normal(size=8),3);ww=group_exp(rng.normal(size=3),2);zz=np.exp(1j*rng.normal())
        minimum=min(minimum,potential(cc,ww,zz,p))
    assert max([gauge_err]+lift_errors)<1e-10 and minimum>0
    return dict(node_gauge_invariance_error=gauge_err,independent_link_lift_error=max(lift_errors),
        minimum_sampled_face_potential=minimum,samples=12,
        positivity_and_zero_set_proven_analytically_not_from_sampling=True)


def remaining_freedom_check():
    p=parameters();C=group_exp(np.array([.2,-.4,.7,.1,.3,.2,-.2,.6]),3)
    W=group_exp([.3,-.5,.7],2);z=np.exp(.37j)
    rows=[]
    for fraction in (.25,.5,.75):
        delta,rest=coefficients(p,fraction)
        residual=float(np.max(abs(delta*p['indices']+rest*np.array([3,2,36])-p['K'])))
        assert min(rest)>0 and residual<1e-12
        rows.append(dict(delta=delta,magnetic_Hessian_residual=residual,
            finite_holonomy_potential=potential(C,W,z,p,fraction)))
    difference=abs(rows[0]['finite_holonomy_potential']-rows[-1]['finite_holonomy_potential'])
    assert difference>1e-3
    return dict(rows=rows,finite_holonomy_difference=difference,
        equal_quadratic_limit_does_not_imply_same_full_quantum_theory=True)


def run():
    checks=[units_check,neutral_symplectic_check,complete_graph_spectrum_check,
        nonlinear_Higgs_hessian_check,centre_obstruction_and_zeros_check,Lie_hessian_check,
        nonlinear_gauge_and_lifts_check,remaining_freedom_check]
    evidence={f.__name__:f() for f in checks}
    deps=['joint_scalar_propagation_matching.py','joint_gauge_matter_propagation.py',
        'protected_pair_rg.py','joint_singlet_common_mass_rg_results.json',
        'joint_gauge_matter_constraints_results.json','joint_yukawa_higgs_matching_results.json',
        'research_note_531.md','research_note_547.md','research_note_553.md','research_note_568.md']
    return dict(round=569,tests_run=len(checks),failures=0,errors=0,checks=[f.__name__ for f in checks],
        evidence=evidence,dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(same_Higgs_vacuum_and_radial_potential=True,residual_photon_Gauss_preserved=True,
            faithful_quotient_magnetic_completion=True,extra_magnetic_function_and_matching_inputs=True,
            no_new_Z_mass_or_hypercharge_theorem=True,no_chiral_fermion_dynamics_implemented=True,
            no_quantum_continuum_dimension_GR_or_unification_completion=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False))
