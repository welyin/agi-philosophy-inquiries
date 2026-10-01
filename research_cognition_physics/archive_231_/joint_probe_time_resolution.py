"""560: full-source time-translation response bound for the same slow probe.

The full unbounded H theorem and Gaussian-polynomial physical clock witness
are analytic. Laguerre matrices are exact Galerkin forms of that H, but their
subspace is NOT invariant; their finite-time output is a diagnostic, not a
converged full-matter simulation. No new Hamiltonian is substituted in proof.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import joint_finite_time_gauge_probe as old
import joint_matter_energy_moment_control as poly
import joint_gauge_link_reference as link

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_probe_time_resolution_results.json'


def response_bound(E,tau,s,c,g=1.,sigma=1.):
    m2=sigma**2; m3=2*math.sqrt(2/math.pi)*sigma**3; m4=3*sigma**4
    first=m2*(c['B0']+c['B1']*E)
    second=m2*(c['B0']+2*c['B1']*E)+2*c['B1']*(E+c['A'])*g*m3/tau+c['B1']*c['B_squared']*g*g*m4/tau**2
    return abs(s)*g/tau*(math.sqrt(first)+math.sqrt(second))


def required_source_energy(tau,s,error,c,g=1.,sigma=1.):
    assert s!=0 and g>0 and sigma>0 and 0<=error<.5
    m2=sigma**2;m3=2*math.sqrt(2/math.pi)*sigma**3;m4=3*sigma**4
    numerator=(1-2*error)**2*tau*tau/(2*s*s*g*g)-2*c['B0']*m2-2*c['B1']*c['A']*g*m3/tau-c['B1']*c['B_squared']*g*g*m4/tau**2
    denominator=c['B1']*(3*m2+2*g*m3/tau)
    return max(0.,numerator/denominator)


def laguerre(x,N,alpha):
    vals=[np.ones_like(x)]
    if N:
        vals.append(1+alpha-x)
        for n in range(1,N):
            vals.append(((2*n+1+alpha-x)*vals[-1]-(n+alpha)*vals[-2])/(n+1))
    return np.array(vals).T


def gamma_quadrature(n):
    j=np.arange(n,dtype=float)
    J=np.diag(2*j+2)+np.diag(np.sqrt((j[1:])*(j[1:]+1)),1)+np.diag(np.sqrt((j[1:])*(j[1:]+1)),-1)
    x,V=np.linalg.eigh(J)
    return x,V[0]**2


def quadrature_forms(N=8,nq=32,b=1.,omega=1.):
    L,u,c=old.constants();a=1/(2*c['v']);t=1/(2*omega);k=c['kappa']*c['v']
    scale=2*k*t;x,w=gamma_quadrature(nq)
    P=laguerre(x,N,1)/np.sqrt(np.arange(N+1)+1)
    dP=np.zeros_like(P)
    dP[:,1:]=-laguerre(x,N-1,2)/np.sqrt(np.arange(1,N+1)+1)
    z2=2*t*x
    Rmean=(z2+4*t)/2
    Rfour=z2*z2/4+3*t*z2+6*t*t
    W=c['v']/2*(L[0,0]*(Rfour-2*u[0]*Rmean+u[0]**2)
        +2*L[0,1]*(Rmean-u[0])*(t-u[1])+L[1,1]*(3*t*t-2*u[1]*t+u[1]**2))
    local_H=10*a*omega-a*omega*omega*(z2+6*t)+W+scale*x
    D=k*(4*a+3*b*t/4)
    H=P.T@((w*local_H)[:,None]*P)+(D/scale)*dP.T@((w*x)[:,None]*dP)
    H=(H+H.T)/2
    F=scale*P.T@((w*x)[:,None]*P);F=(F+F.T)/2
    G=P.T@(w[:,None]*P)
    return H,F,G


def compression(N=8,nq=32,b=1.,omega=1.,eta=.3):
    L,u,c=old.constants();a=1/(2*c['v']);t=1/(2*omega);k=c['kappa']*c['v'];scale=2*k*t
    n=np.arange(N+1,dtype=float)
    X=np.diag(2*n+2)-np.diag(np.sqrt(n[1:]*(n[1:]+1)),1)-np.diag(np.sqrt(n[1:]*(n[1:]+1)),-1)
    X2=X@X;X2[-1,-1]+=(N+1)*(N+2)
    c2=c['v']/2*L[0,0]*t*t
    c1=c['v']/2*(L[0,0]*(6*t*t-2*u[0]*t)+2*L[0,1]*t*(t-u[1]))-2*a*omega*omega*t+scale
    c0=c['v']/2*(L[0,0]*(6*t*t-4*u[0]*t+u[0]**2)+2*L[0,1]*(2*t-u[0])*(t-u[1])+L[1,1]*(3*t*t-2*u[1]*t+u[1]**2))+10*a*omega-6*a*omega*omega*t
    D=k*(4*a+3*b*t/4)
    H=c0*np.eye(N+1)+c1*X+c2*X2+D/scale*np.diag(n)
    F=scale*X
    assert np.linalg.eigvalsh(H)[0]>0 and np.linalg.eigvalsh(F)[0]>0
    assert np.linalg.eigvalsh(c['B0']*np.eye(N+1)+c['B1']*H-F@F)[0]>0
    psi=np.zeros(N+1,complex);psi[0]=1;psi[1]=-1j*eta*scale*math.sqrt(2)
    psi/=np.linalg.norm(psi)
    return H,F,psi,c,dict(omega=omega,eta=eta,b=b,scale=scale,
        ground_amplitude_energy=float(H[0,0]),witness_energy=float(np.vdot(psi,H@psi).real),
        witness_F_velocity=float((1j*np.vdot(psi,(H@F-F@H)@psi)).real),
        exact_Laguerre_multiplication_with_X_squared_boundary=True)


def full_source_witness():
    L,u,c=old.constants();d=10;omega=1.;t=.5;b=1.;eta=.3;a=.5;k=c['kappa']*c['v']
    q=[poly.var(i,d) for i in range(d)]
    radii=[poly.add(*(poly.mul(q[j],q[j]) for j in range(i,i+4))) for i in (0,5)]
    W={}
    for i,o in enumerate((0,5)):
        ds=[poly.add(radii[i],poly.const(-u[0],d)),poly.add(poly.mul(q[o+4],q[o+4]),poly.const(-u[1],d))]
        W=poly.add(W,poly.scale(poly.quadratic(ds,L),c['v']/4))
    diff=[poly.add(q[i],poly.scale(q[5+i],-1)) for i in range(4)]
    F=poly.scale(poly.add(*(poly.mul(z,z) for z in diff)),k/2)
    f=poly.add(F,poly.const(-4*k*t,d))
    norm=1+eta*eta*8*k*k*t*t
    P=poly.scale(poly.add(poly.const(1,d),poly.scale(f,1j*eta)),1/math.sqrt(norm))
    assert abs(poly.norm2(P,t)-1)<1e-12
    derivatives=[poly.add(poly.deriv(P,i),poly.scale(poly.mul(q[i],P),-omega)) for i in range(d)]
    df=[poly.deriv(F,i) for i in range(d)]
    dg=[poly.scale(poly.add(*(poly.scale(poly.mul(q[i],q[5+j]),gen[i,j])
        for i in range(4) for j in range(4) if gen[i,j])),-k) for gen in link.J]
    dgP=[poly.scale(z,1j*eta/math.sqrt(norm)) for z in dg]
    density=poly.mul(poly.conj(P),P)
    E=a*sum(poly.norm2(z,t) for z in derivatives)+b*sum(poly.norm2(z,t) for z in dgP)
    E+=poly.real(poly.expect(poly.mul(density,poly.add(W,F)),t))
    velocity=2*a*sum(poly.expect(poly.mul(poly.mul(poly.conj(P),dp),dz),t).imag for dp,dz in zip(derivatives,df))
    velocity+=2*b*sum(poly.expect(poly.mul(poly.mul(poly.conj(P),dp),dz),t).imag for dp,dz in zip(dgP,dg))
    G=16*a*k*k*t+3*b*k*k*t*t
    closed_velocity=2*eta*G/norm
    H,Fm,ps,c,info=compression(4,24,b,omega,eta)
    assert abs(velocity-closed_velocity)<1e-11
    assert abs(E-info['witness_energy'])<1e-10
    assert abs(velocity-info['witness_F_velocity'])<1e-10
    # Source is phi*(1+i eta(F-meanF)), with real gauge-invariant F and phi.
    # No exact time propagation is needed for its full-H instantaneous current.
    return dict(full_source_E1=E,full_source_initial_F_velocity=velocity,
        analytic_velocity=closed_velocity,base_F_mean=4*k*t,base_F_variance=8*k*k*t*t,
        base_gradient_expectation=G,normalization=norm,eta=eta,
        polynomial_vs_Galerkin_energy_residual=abs(E-info['witness_energy']),
        all_source_coordinates_and_group_gradients_included=True,
        witness_is_normalized_Gauss_physical_and_nonstationary=True)


def evolve(H,v,time):
    E,V=np.linalg.eigh(H)
    return V@(np.exp(-1j*time*E)*(V.conj().T@v))


def output_case(tau,N=8,nq=48,s=.4):
    H,F,psi,c,info=compression(N,max(32,2*N+8))
    shifted=evolve(H,psi,s); initial_gap=float(np.vdot(shifted,F@shifted).real-np.vdot(psi,F@psi).real)
    nodes,w=np.polynomial.hermite.hermgauss(nq);w/=math.sqrt(math.pi)
    wave0=[];wave1=[];joint2=0.;means=np.zeros(2);cost=np.zeros(2);comm_res=0.
    for p,weight in zip(math.sqrt(2)*nodes,w):
        U0,d0=old.probe.propagate(p,tau,H,F,psi)
        U1,d1=old.probe.propagate(p,tau,H,F,shifted)
        aligned=evolve(H,U0,s)
        joint2+=weight*np.linalg.norm(U1-aligned)**2
        wave0.append(math.sqrt(weight)*U0);wave1.append(math.sqrt(weight)*U1)
        for j,(v,dv) in enumerate(((U0,d0),(U1,d1))):
            means[j]+=weight*(1j*np.vdot(v,dv)).real
            cost[j]+=weight*(np.vdot(v,H@v).real-info['witness_energy'])
        # Independent commutator identity, checked on the complete finite basis.
        ev,V=np.linalg.eigh(H+p*F/tau);U=(V*np.exp(-1j*tau*ev))@V.T
        comm_res=max(comm_res,float(np.linalg.norm(H@U-U@H+p/tau*(F@U-U@F))))
    A=np.array(wave0);B=np.array(wave1)
    rho0=A@A.conj().T;rho1=B@B.conj().T
    distance=float(np.sum(abs(np.linalg.eigvalsh(rho1-rho0)))/2)
    bound=response_bound(info['witness_energy'],tau,s,c)
    assert comm_res<1e-9 and distance<=math.sqrt(joint2)+1e-10
    assert math.sqrt(joint2)<=bound+1e-10
    assert abs(np.trace(rho0)-1)<1e-12 and abs(np.trace(rho1)-1)<1e-12
    return dict(duration=tau,basis_degree=N,pointer_quadrature=nq,source_time_shift=s,
        initial_source_F_difference=initial_gap,source_E1=info['witness_energy'],
        pointer_state_half_trace_distance_quadrature=distance,
        aligned_joint_state_norm=math.sqrt(joint2),full_model_response_upper_bound=bound,
        any_pointer_readout_binary_error_lower_bound=(1-min(1.,bound))/2,
        actual_compensated_pointer_means=means.tolist(),source_energy_changes=cost.tolist(),
        commutator_identity_residual=comm_res)


def run():
    witness=full_source_witness()
    H,F,psi,c,info=compression(8,32)
    form_errors=[]
    for nq in (24,32):
        H2,F2,G=quadrature_forms(6,nq)
        H6,F6,_,_,_=compression(6)
        form_errors.extend([np.max(abs(H6-H2)),np.max(abs(F6-F2)),np.max(abs(G-np.eye(7)))])
    quad_difference=float(max(form_errors))
    assert quad_difference<1e-8
    rows=[output_case(t) for t in (10.,100.,1000.)]
    errors=[]
    for row in rows:
        other=output_case(row['duration'],nq=64)
        errors.append(abs(row['pointer_state_half_trace_distance_quadrature']-other['pointer_state_half_trace_distance_quadrature']))
    assert max(errors)<2e-7
    assert all(abs(row['initial_source_F_difference'])>.05 for row in rows)
    assert rows[-1]['pointer_state_half_trace_distance_quadrature']<.001
    # Changing basis is only a Galerkin diagnostic: never assert convergence of
    # this incomplete invariant-variable family to the full source dynamics.
    other_basis=output_case(100.,N=12,nq=64)
    energy_rows=[]
    for tau in (100.,1000.,10000.):
        minimum=required_source_energy(tau,.4,.1,c)
        energy_rows.append(dict(duration=tau,necessary_E1_at_binary_error_point_one=minimum,
            bound_at_fixed_E1_18=response_bound(18.,tau,.4,c)))
        if minimum>0:
            m3=2*math.sqrt(2/math.pi)
            bracket=c['B1']*minimum*(3+2*m3/tau)+2*c['B0']+2*c['B1']*c['A']*m3/tau+3*c['B1']*c['B_squared']/tau**2
            assert abs(2*.4**2/tau**2*bracket-.8**2)<1e-12
    for row in rows:
        assert row['duration']>math.sqrt(c['B1']/2)
    checks=['normalized_full_Gauss_clock_source_and_exact_probability_current',
        'independent_polynomial_and_conditional_Gamma_energy_forms',
        'same_H_compressed_probe_joint_evolution_and_commutator',
        'actual_pointer_density_response_and_independent_momentum_quadrature',
        'finite_energy_time_response_bound_without_ground_state_or_gap',
        'binary_error_required_source_budget_and_stable_pointer_duration']
    deps=('joint_finite_time_gauge_probe.py','joint_finite_time_gauge_probe_results.json',
        'joint_gauge_link_reference.py','joint_matter_energy_moment_control.py','joint_singlet_common_mass_rg_results.json')
    return dict(round=560,tests_run=len(checks),failures=0,errors=0,checks=checks,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        full_source_clock_witness=witness,conditional_form_quadrature_difference=quad_difference,
        compressed_probe_rows=rows,independent_pointer_quadrature_max_difference=max(errors),
        changed_Galerkin_degree_diagnostic=other_basis,necessary_source_budget=energy_rows,
        scope=dict(full_source_theorem_analytic_and_energy_uniform=True,
            time_shift_fixed_or_in_fixed_bounded_interval=True,
            source_may_be_nonstationary_and_correlated_with_old_reference=True,
            only_pointer_and_old_reference_readout_not_direct_source=True,
            fixed_single_probe_family_gain_and_pointer_moments=True,
            Galerkin_subspace_not_invariant_or_full_dynamics_convergence=True,
            budget_lower_bound_necessary_not_sufficient_or_universal_clock_bound=True,
            no_spacetime_gravity_or_unified_completion=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else: assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=560,tests=result['tests_run'],witness=result['full_source_clock_witness'],
        rows=result['compressed_probe_rows'],budget=result['necessary_source_budget']),ensure_ascii=False))
