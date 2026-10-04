"""567: match the existing graph matter to the two-scalar continuum candidate.

Classical fixed-background quadratic interface, not a quantum continuum limit.
All parameters are inherited; the extra singlet gradient is explicitly an input.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_finite_time_gauge_probe as inherited

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_scalar_propagation_matching_results.json'


def parameters():
    L, u, _ = inherited.constants()
    D = np.diag(np.sqrt(u))
    return L, u, 2 * D @ L @ D


def hessian(fun, x, step=1e-4):
    ans = np.zeros((len(x), len(x)))
    e = np.eye(len(x)) * step
    for i in range(len(x)):
        ans[i, i] = (fun(x+e[i])-2*fun(x)+fun(x-e[i])) / step**2
        for j in range(i):
            ans[i, j] = ans[j, i] = (fun(x+e[i]+e[j])-fun(x+e[i]-e[j])
                -fun(x-e[i]+e[j])+fun(x-e[i]-e[j])) / (4*step**2)
    return ans


def sector_check():
    L, u, _ = parameters()
    # Left multiplication by i,j,k on real quaternions; twice the Lie generators.
    generators = np.array([
        [[0,-1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,1,0]],
        [[0,0,-1,0],[0,0,0,1],[1,0,0,0],[0,-1,0,0]],
        [[0,0,0,-1],[0,0,-1,0],[0,1,0,0],[1,0,0,0]]], dtype=float)/2
    n = np.array([1.,2.,3.,4.]); n /= np.linalg.norm(n)
    h = np.array([.71,.82,.95,.64]); s = np.array([.4,.6,.51,.7])
    X = h[:,None]*n; pX = np.array([.2,-.1,.7,.3])[:,None]*n
    edges = [(0,1),(0,2),(2,3)]; k = .8; volume = 1.3
    d = np.column_stack((h*h-u[0], s*s-u[1]))
    onsite = d @ L.T
    grad = volume*onsite[:,0,None]*X
    link_torques = []
    for i,j in edges:
        diff = X[i]-X[j]
        grad[i] += k*diff; grad[j] -= k*diff
        link_torques.extend(-k*diff@(T@X[j]) for T in generators)
    projected = grad-(grad@n)[:,None]*n
    charges = np.array([[pX[i]@(T@X[i]) for T in generators] for i in range(4)])
    def energy(z):
        x=z.reshape(4,4)
        ds=np.column_stack((np.sum(x*x,axis=1)-u[0],s*s-u[1]))
        return volume*np.einsum('ni,ij,nj->',ds,L,ds)/4 + sum(k*np.sum((x[i]-x[j])**2)/2 for i,j in edges)
    step=1e-5; flat=X.ravel(); fd=np.zeros(16)
    for a in range(16):
        e=np.eye(16)[a]*step; fd[a]=(energy(flat+e)-energy(flat-e))/(2*step)
    residual=float(np.max(abs(fd-grad.ravel())))
    zeros=max(float(np.max(abs(projected))),float(np.max(abs(charges))),float(np.max(np.abs(link_torques))))
    assert residual<1e-9 and zeros<1e-14
    return dict(Cartesian_force_finite_difference_residual=residual,
        transverse_force_charge_link_torque_max=zeros,edges=[list(edge) for edge in edges],
        scope='classical invariant sector; not a normalizable sharp-link quantum state')


def graph_check():
    L,u,mass=parameters(); N=5; k=.8; volume=1.3
    edges=[(i,(i+1)%N) for i in range(N)]
    lap=np.zeros((N,N))
    for i,j in edges:
        lap[i,i]+=1; lap[j,j]+=1; lap[i,j]-=1; lap[j,i]-=1
    centre=np.tile(np.sqrt(u),N)
    errors=[]; spectral=[]
    for ks in (0.,k):
        def potential(z):
            phi=z.reshape(N,2); d=phi*phi-u
            return np.einsum('ni,ij,nj->',d,L,d)/4 + sum(
                (k*(phi[i,0]-phi[j,0])**2+ks*(phi[i,1]-phi[j,1])**2)/(2*volume) for i,j in edges)
        predicted=np.kron(np.eye(N),mass)+np.kron(lap,np.diag([k,ks])/volume)
        errors.append(float(np.max(abs(hessian(potential,centre)-predicted))))
        fourier=np.sort(np.concatenate([np.linalg.eigvalsh(mass+4*np.sin(np.pi*j/N)**2*np.diag([k,ks])/volume) for j in range(N)]))
        spectral.append(float(np.max(abs(fourier-np.linalg.eigvalsh(predicted)))))
    assert max(errors)<2e-8 and max(spectral)<2e-14
    return dict(sites=N,full_potential_Hessian_errors=errors,full_matrix_Fourier_errors=spectral,
        original_and_repaired_common_vacuum=np.sqrt(u).tolist())


def scales_check():
    _,_,mass=parameters(); c2=.8; p=np.array([.3,.7,1.1]); p2=float(p@p)
    rows=[]
    for eps in (1.,.5,.25,.125):
        lam=float(4*np.sum(np.sin(eps*p/2)**2)); q2=lam/eps**2
        bound=eps**2*float(np.sum(p**4))/12
        old=np.linalg.eigvalsh(mass+np.diag([c2*q2,0.]))
        limit=np.linalg.eigvalsh(mass+np.diag([c2*p2,0.]))
        err=float(np.max(abs(old-limit)))
        assert 0<=p2-q2<=bound+1e-14 and err<=c2*bound+1e-14
        volume=eps**3
        rows.append(dict(epsilon=eps,lattice_q_squared=q2,continuum_q_squared=p2,
            symbol_error=p2-q2,symbol_error_bound=bound,dispersion_error=err,
            fixed_kappa_speed_squared=c2*eps**2,
            fixed_speed_required_k=c2*volume/eps**2))
    assert rows[-1]['dispersion_error']<rows[0]['dispersion_error']/50
    return dict(rows=rows,dimension_three_is_input=True,
        statement='bounded physical momentum symbol convergence only; no quantum limit')


def principal_check():
    _,_,mass=parameters(); c2=.8; P=np.diag([c2,0.])
    eig,basis=np.linalg.eigh(mass); slopes=c2*basis[0,:]**2
    A,B,C=mass[0,0],mass[0,1],mass[1,1]
    exact_slopes=c2/2*(1+np.array([-1.,1.])*(A-C)/np.sqrt((A-C)**2+4*B*B))
    assert np.max(abs(slopes-exact_slopes))<1e-14
    congruence_errors=[]
    for J in (np.array([[1.,.3],[.2,1.4]]),np.array([[2.,-1.],[.4,.9]])):
        K=J.T@J; G=J.T@P@J
        speeds=np.sort(np.linalg.eigvals(np.linalg.solve(K,G)).real)
        congruence_errors.append(float(np.max(abs(speeds-[0.,c2]))))
        assert np.linalg.matrix_rank(G,tol=1e-12)==1
    sym=np.array([[2.,.5],[.5,2.]])
    vals,vec=np.linalg.eigh(sym); symmetric_slopes=c2*vec[0,:]**2
    nonlinear=np.linalg.eigvalsh(sym+P)-(vals+c2/2)
    assert np.max(abs(symmetric_slopes-c2/2))<1e-14 and max(abs(nonlinear))>.1
    assert max(congruence_errors)<1e-14
    rows=[]
    for q2 in (0.,.1,1.,100.):
        formula=(A+C+c2*q2+np.array([-1.,1.])*np.sqrt((A-C+c2*q2)**2+4*B*B))/2
        direct=np.linalg.eigvalsh(mass+q2*P)
        assert max(abs(formula-direct))<2e-14
        rows.append(dict(q_squared=q2,omega_squared=direct.tolist()))
    return dict(mass_squared_matrix=mass.tolist(),mass_squared_eigenvalues=eig.tolist(),
        d_omega_squared_d_q_squared_at_zero=slopes.tolist(),
        generalized_speed_congruence_errors=congruence_errors,
        equal_diagonal_mass_slopes=symmetric_slopes.tolist(),
        equal_slope_example_departure_at_q_squared_one=nonlinear.tolist(),rows=rows,
        massive_group_velocity_at_zero=0.)


def repair_check():
    _,_,mass=parameters(); c2=.8; m2=np.linalg.eigvalsh(mass)
    error=0.
    for q2 in (0.,.01,1.,100.):
        error=max(error,float(np.max(abs(np.linalg.eigvalsh(mass+c2*q2*np.eye(2))-(m2+c2*q2)))))
    assert error<2e-14
    # The two eigenvalues of K^{-1}G, unlike the mass basis, determine common principal speed.
    branches=[]
    for cs2 in (0.,.4,.8):
        G=np.diag([c2,cs2]); speeds=np.linalg.eigvalsh(G)
        branches.append(dict(s_speed_squared=cs2,spatial_rank=int(np.linalg.matrix_rank(G)),
            shared_nonzero_principal_speed=bool(cs2==c2 and cs2>0),speeds_squared=speeds.tolist()))
    return dict(exact_repaired_spectral_translation_residual=error,branches=branches,
        singlet_edge_rule_and_speed_matching_are_added_inputs=True,
        preserved_on_site_mass_and_vacuum=True)


def single_mode_check():
    _,_,mass=parameters(); A,B,C=mass[0,0],mass[0,1],mass[1,1]; c2=.8
    delta=.25; rows=[]; identity_error=0.
    for q2,w2 in ((0.,.02),(.02,0.),(.04,.01),(.01,.04)):
        x=c2*q2-w2; assert abs(x)<=delta*A
        exact=C-w2-B*B/(A+x)
        approx=C-B*B/A+B*B*c2*q2/A**2-(1+B*B/A**2)*w2
        remainder=-B*B*x*x/(A*A*(A+x))
        bound=B*B*x*x/(A**3*(1-delta))
        assert abs(exact-approx-remainder)<1e-15 and abs(remainder)<=bound+1e-15
        det=np.linalg.det(mass+np.diag([c2*q2,0.])-w2*np.eye(2))
        identity_error=max(identity_error,float(abs((A+x)*exact-det)))
        rows.append(dict(q_squared=q2,omega_squared=w2,exact_inverse=exact,
            derivative_approximation=approx,remainder=remainder,remainder_bound=bound))
    ratio=float(np.linalg.eigvalsh(mass)[0]/A)
    scaling=[float(np.linalg.eigvalsh(scale*mass)[0]/(scale*A)) for scale in (.1,1.,10.)]
    assert identity_error<1e-15 and max(abs(ratio-np.array(scaling)))<1e-14
    return dict(rows=rows,Schur_determinant_residual=identity_error,
        induced_spatial_coefficient=B*B*c2/A**2,time_coefficient=1+B*B/A**2,
        effective_single_mode_speed_squared=(B*B*c2/A**2)/(1+B*B/A**2),
        light_mass_squared_over_eliminated_diagonal=ratio,common_mass_scaling_ratios=scaling,
        elimination_requires_homogeneous_source_or_explicit_initial_data=True,
        not_a_two_independent_field_completion=True)


def run():
    checks=[sector_check,graph_check,scales_check,principal_check,repair_check,single_mode_check]
    evidence={f.__name__:f() for f in checks}
    deps=['joint_finite_time_gauge_probe.py','joint_matter_energy_moment_control.py',
          'research_note_548.md','research_note_549.md','research_note_566.md']
    return dict(round=567,tests_run=len(checks),failures=0,errors=0,
        checks=[f.__name__ for f in checks],evidence=evidence,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        scope=dict(classical_fixed_background_quadratic_interface=True,
            not_a_quantum_Gauss_source_or_continuum_limit=True,
            extra_s_gradient_and_lattice_scaling_are_explicit_inputs=True,
            single_low_energy_induced_propagation_not_excluded=True,
            loops_nonminimal_gravity_and_other_backgrounds_not_excluded=True,
            no_dimension_GR_or_unified_completion=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true'); args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False))
