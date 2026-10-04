"""570: nonlinear compatibility of the same finite Higgs-gauge phase space.

Small Gauss residual is not a uniform distance certificate near a stabilizer.
The paired-loop construction satisfies exact Gauss with its full internal source.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_quotient_gauge_completion as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauss_nonlinear_compatibility_results.json'
T=old.generators(2)


def realify(x):return np.r_[x.real,x.imag]


def complexify(x):return x[:2]+1j*x[2:]


def adjoint(U):
    return np.array([[2*np.trace(a@U@b@U.conj().T).real for b in T] for a in T])


def cycle(n=4):return [(i,(i+1)%n) for i in range(n)], [[(i,1) for i in range(n)]]


def paired_graph():
    edges,faces=cycle();edges=edges+[(i+4,j+4) for i,j in edges]+[(0,4)]
    return edges,faces+[[(i+4,1) for i in range(4)]]


def configuration(eta,paired=False):
    p=old.parameters();edges,faces=paired_graph() if paired else cycle()
    N=max(max(e) for e in edges)+1;E=len(edges)
    angles=np.r_[np.ones(4),-np.ones(4),0.] if paired else np.ones(E)
    U=np.array([old.group_exp([eta*t,0,0],2) for t in angles])
    X=np.tile(np.array([0,np.sqrt(p['h2'])],complex),(N,1))
    return dict(par=p,edges=edges,faces=faces,N=N,E=E,X=X,U=U,z=np.ones(E,complex))


def incidence(q):
    D=np.zeros((q['E'],q['N']))
    for k,(i,j) in enumerate(q['edges']):D[k,i]=1;D[k,j]=-1
    return D


def generator_matrix(q):
    N,E=q['N'],q['E'];B=np.zeros((4*N+4*E,4*N))
    for i,X in enumerate(q['X']):
        for a in range(3):B[4*i:4*i+4,4*i+a]=realify(1j*T[a]@X)
        B[4*i:4*i+4,4*i+3]=realify(3j*X)
    for k,(i,j) in enumerate(q['edges']):
        B[4*N+3*k:4*N+3*k+3,4*i:4*i+3]=np.eye(3)
        B[4*N+3*k:4*N+3*k+3,4*j:4*j+3]=-adjoint(q['U'][k])
        B[4*N+3*E+k,4*i+3]=1;B[4*N+3*E+k,4*j+3]=-1
    return B


def kinetic_diagonal(q):
    _,bw,b0=q['par']['b'];N,E=q['N'],q['E']
    return np.r_[np.ones(4*N),np.full(3*E,2*bw),np.full(E,2*b0)]


def pack(q,P,weak,circle):return np.r_[np.concatenate([realify(x) for x in P]),weak.ravel(),circle]


def unpack(q,p):
    N,E=q['N'],q['E']
    return np.array([complexify(x) for x in p[:4*N].reshape(N,4)]),p[4*N:4*N+3*E].reshape(E,3),p[4*N+3*E:]


def direct_gauss(q,p):
    P,weak,circle=unpack(q,p);G=np.zeros((q['N'],4))
    for i in range(q['N']):
        for a in range(3):G[i,a]=np.vdot(P[i],1j*T[a]@q['X'][i]).real
        G[i,3]=np.vdot(P[i],3j*q['X'][i]).real
    for e,(i,j) in enumerate(q['edges']):
        G[i,:3]+=weak[e];G[j,:3]-=adjoint(q['U'][e]).T@weak[e]
        G[i,3]+=circle[e];G[j,3]-=circle[e]
    return G


def trial(q,amplitude):
    weak=np.zeros((q['E'],3));weak[:8 if q['N']==8 else q['E'],1]=amplitude
    return pack(q,np.zeros_like(q['X']),weak,np.zeros(q['E']))


def project(q,p):
    root=np.sqrt(kinetic_diagonal(q));z=root*p
    A=generator_matrix(q).T/root[None,:]
    u,s,vh=np.linalg.svd(A,full_matrices=False);active=s>1e-12
    vertical=vh[active].T@(vh[active]@z);zp=z-vertical
    corrected=zp/root
    return corrected,dict(residual=float(np.linalg.norm(A@zp)),
        removed_norm=float(np.linalg.norm(vertical)),initial_norm=float(np.linalg.norm(z)),
        singular_values=s,rank=int(np.sum(active)),orthogonality=float(abs(zp@vertical)),
        energy_before=float(z@z/2),energy_after=float(zp@zp/2))


def potential(q):
    val=0.
    for e,(i,j) in enumerate(q['edges']):
        val+=np.linalg.norm(q['X'][i]-q['z'][e]**3*q['U'][e]@q['X'][j])**2/2
    for face in q['faces']:
        W=np.eye(2,dtype=complex);z=1.+0j
        for e,sign in face:
            W=W@(q['U'][e] if sign==1 else q['U'][e].conj().T);z*=q['z'][e]**sign
        val+=old.potential(np.eye(3),W,z,q['par'])
    # All configurations below retain |X|=h_star and s=s_star, so V and s gradients vanish.
    return float(val)


def corrected_pair(eta):
    q=configuration(eta,True);p=trial(q,eta);r=direct_gauss(q,p)[:,:3]
    P=np.array([sum(-4*r[i,a]/q['par']['h2']*(1j*T[a]@q['X'][i]) for a in range(3)) for i in range(q['N'])])
    f=np.array([-9,-3,3,9,9,3,-3,-9,24],float)
    circle=eta*np.sin(eta)*f
    _,weak,_=unpack(q,p)
    return q,p,pack(q,P,weak,circle),r


def generator_and_momentum_check():
    q=configuration(.31,True);B=generator_matrix(q);rng=np.random.default_rng(570)
    p=rng.normal(size=B.shape[0]);xi=rng.normal(size=B.shape[1]).reshape(q['N'],4)
    direct=direct_gauss(q,p).ravel();err=float(np.max(abs(direct-B.T@p)))
    step=1e-6
    def gauge_at(t):
        gs=[old.group_exp(t*x[:3],2) for x in xi];zs=np.exp(1j*t*xi[:,3])
        X=np.array([zs[i]**3*gs[i]@q['X'][i] for i in range(q['N'])])
        U=np.array([gs[i]@q['U'][e]@gs[j].conj().T for e,(i,j) in enumerate(q['edges'])])
        alpha=np.array([t*(xi[i,3]-xi[j,3]) for i,j in q['edges']])
        return X,U,alpha
    xp,up,ap=gauge_at(step);xm,um,am=gauge_at(-step)
    dx=np.concatenate([realify(x) for x in (xp-xm)/(2*step)])
    du=[]
    for e in range(q['E']):
        left=(up[e]-um[e])/(2*step)@q['U'][e].conj().T
        du.extend([(-2j*np.trace(t@left)).real for t in T])
    fd=np.r_[dx,du,(ap-am)/(2*step)]
    fd_error=float(np.max(abs(fd-B@xi.ravel())))
    assert err<1e-12 and fd_error<2e-9
    return dict(moment_map_independent_error=err,finite_gauge_derivative_error=fd_error)


def weighted_projection_check():
    rng=np.random.default_rng(5701);rows=[]
    for eta in (0.,.17):
        q=configuration(eta);p=rng.normal(size=4*q['N']+4*q['E'])
        corrected,info=project(q,p);res=float(np.max(abs(direct_gauss(q,corrected))))
        correction=float(np.sum(kinetic_diagonal(q)*(p-corrected)**2)/2)
        identity=abs(info['energy_before']-info['energy_after']-correction)
        assert res<1e-12 and identity<2e-12 and info['energy_after']<=info['energy_before']
        rows.append(dict(eta=eta,rank=info['rank'],constraint_residual=res,
            energy_drop=correction,orthogonal_energy_identity_error=identity))
    assert rows[0]['rank']==15 and rows[1]['rank']==16
    return dict(rows=rows,scope='fixed-configuration mathematical projection, not a preparation protocol')


def quadratic_obstruction_check():
    rows=[];q0=configuration(0.);p1=trial(q0,1.)
    assert np.linalg.norm(direct_gauss(q0,p1))==0
    xi=np.tile([0.,0.,6.,1.],q0['N']);assert np.linalg.norm(generator_matrix(q0)@xi)<1e-13
    for eta in (.2,.1,.03,.01):
        q=configuration(eta);p=trial(q,eta);G=direct_gauss(q,p)
        expected=np.tile([0.,eta*(1-np.cos(eta)),-eta*np.sin(eta),0.],(4,1))
        Q=float(xi@G.ravel());formula=-24*eta*np.sin(eta)
        assert np.max(abs(G-expected))<1e-14 and abs(Q-formula)<1e-14
        rows.append(dict(eta=eta,Gauss_norm=float(np.linalg.norm(G)),
            stabilizer_charge=Q,charge_over_eta_squared=Q/eta**2))
    step=1e-5
    derivative=(generator_matrix(configuration(step))-generator_matrix(configuration(-step)))/(2*step)
    coefficient=float(p1@derivative@xi)
    assert abs(coefficient+24)<1e-8
    return dict(first_order_Gauss_residual=0.,quadratic_compatibility_coefficient=coefficient,
        exact_coefficient=-24,rows=rows,nonzero_coefficient_precludes_exact_C1_curve_with_this_tangent=True)


def small_residual_distance_check():
    rows=[]
    for eta in (.2,.05,.01,.003):
        q=configuration(eta)
        for amp in (eta,1.):
            p=trial(q,amp);_,info=project(q,p)
            relative=info['removed_norm']/info['initial_norm'];bound=np.cos(eta/2)
            residual=float(np.linalg.norm(direct_gauss(q,p)))
            assert relative>=bound-1e-10 and relative<=1+1e-12 and info['residual']<1e-11
            rows.append(dict(eta=eta,amplitude=amp,Gauss_norm=residual,
                minimum_correction_norm=info['removed_norm'],relative_correction=relative,
                proven_relative_lower_bound=float(bound),full_trial_energy=info['energy_before']+potential(q)))
    return dict(rows=rows,no_uniform_residual_to_distance_bound_across_stabilizer_change=True,
        not_a_counterexample_to_a_specific_constraint_preserving_continuum_sampling=True)


def exact_paired_completion_check():
    rows=[]
    for eta in (.2,.1,.05,.01):
        q,p,pc,r=corrected_pair(eta);res=float(np.max(abs(direct_gauss(q,pc))))
        A=kinetic_diagonal(q);extra=float(np.sum(A*(pc-p)**2)/2)
        b0=q['par']['b'][2];h2=q['par']['h2']
        expected=64/h2*eta**2*np.sin(eta/2)**2+936*b0*eta**2*np.sin(eta)**2
        delta_norm=np.sqrt(2*extra);bridge=float(unpack(q,pc)[2][-1])
        assert res<2e-14 and abs(extra-expected)<1e-13
        assert abs(bridge-24*eta*np.sin(eta))<1e-14
        assert abs(np.sum(A*pc**2)/2-np.sum(A*p**2)/2-extra)<1e-14
        rows.append(dict(eta=eta,full_Gauss_residual=res,extra_kinetic_energy=extra,
            extra_energy_over_eta_fourth=extra/eta**4,correction_norm_over_eta_squared=delta_norm/eta**2,
            internal_bridge_circle_momentum=bridge,total_energy=float(np.sum(A*pc**2)/2+potential(q))))
    return dict(rows=rows,exact_circle_flow_coefficients=[-9,-3,3,9,9,3,-3,-9,24],
        same_first_order_tangent_preserved=True,all_countercharge_and_flux_internal=True,
        extra_energy_is_initial_data_energy_not_a_protocol_cost=True)


def gauge_covariance_check():
    q,_,p,_=corrected_pair(.13);rng=np.random.default_rng(5702)
    gs=[old.group_exp(rng.normal(size=3),2) for _ in range(q['N'])]
    zs=np.exp(1j*rng.normal(size=q['N']));P,weak,circle=unpack(q,p)
    qg=dict(q);qg['X']=np.array([zs[i]**3*gs[i]@q['X'][i] for i in range(q['N'])])
    qg['U']=np.array([gs[i]@q['U'][e]@gs[j].conj().T for e,(i,j) in enumerate(q['edges'])])
    qg['z']=np.array([zs[i]*q['z'][e]/zs[j] for e,(i,j) in enumerate(q['edges'])])
    pg=pack(qg,np.array([zs[i]**3*gs[i]@P[i] for i in range(q['N'])]),
        np.array([adjoint(gs[i])@weak[e] for e,(i,j) in enumerate(q['edges'])]),circle)
    err=float(np.max(abs(direct_gauss(qg,pg))))
    H=lambda a,b:float(np.sum(kinetic_diagonal(a)*b*b)/2+potential(a))
    energy_err=abs(H(qg,pg)-H(q,p))
    assert err<1e-13 and energy_err<1e-11
    return dict(full_Gauss_after_arbitrary_node_gauge_error=err,full_H_energy_invariance_error=energy_err)


def stabilizer_gap_check():
    rows=[];q0=configuration(0.)
    for eta in (.2,.05,.01,.003):
        q=configuration(eta);root=np.sqrt(kinetic_diagonal(q));A=generator_matrix(q).T/root
        singular=np.linalg.svd(A,compute_uv=False);sigma=float(singular[-1]);bw=q['par']['b'][1]
        upper=12*np.sin(eta/2)/np.sqrt(74*bw)
        assert sigma>0 and sigma<=upper+1e-13
        rows.append(dict(eta=eta,least_singular_value=sigma,explicit_upper_bound=float(upper),sigma_over_eta=sigma/eta))
    return dict(vacuum_rank=int(np.linalg.matrix_rank(generator_matrix(q0))),
        perturbed_rank=4*q0['N'],rows=rows,
        condition_is_a_model_error_bound_not_a_new_cognitive_axiom=True)


def commuting_positive_control_check():
    rows=[]
    for eta in (.07,.2):
        q=configuration(eta);weak=np.tile([eta,0,0],(q['E'],1))
        p=pack(q,np.zeros_like(q['X']),weak,np.zeros(q['E']))
        res=float(np.max(abs(direct_gauss(q,p))))
        assert res<1e-14
        rows.append(dict(eta=eta,exact_Gauss_residual=res,total_energy=float(np.sum(kinetic_diagonal(q)*p*p)/2+potential(q))))
    return dict(rows=rows,single_loop_is_not_generically_forbidden=True)


def run():
    checks=[generator_and_momentum_check,weighted_projection_check,quadratic_obstruction_check,
        small_residual_distance_check,exact_paired_completion_check,gauge_covariance_check,
        stabilizer_gap_check,commuting_positive_control_check]
    evidence={f.__name__:f() for f in checks}
    deps=['joint_quotient_gauge_completion.py','joint_quotient_gauge_completion_results.json',
        'research_note_549.md','research_note_550.md','research_note_568.md','research_note_569.md']
    return dict(round=570,tests_run=len(checks),failures=0,errors=0,checks=[f.__name__ for f in checks],
        evidence=evidence,dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(same_finite_Higgs_hypercharge_Hamiltonian=True,
            quadratic_Gauss_compatibility_is_needed_beyond_linear_modes=True,
            exact_internal_paired_completion=True,known_linearization_instability_tools_not_claimed_new=True,
            no_specific_continuum_sampling_failure_proven=True,
            no_spacetime_dimension_GR_chiral_matter_or_quantum_continuum_completion=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False))
