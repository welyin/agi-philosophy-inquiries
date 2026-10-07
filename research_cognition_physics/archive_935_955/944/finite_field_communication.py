"""944: finite-time scalar-field communication, recoil, and finite-mass control.

Smooth massive Nelson-type effective matter is an explicit physical input.
This is not a simulation or a certified matching of full Einstein--Dirac QFT.
No Fock occupation truncation or point-source/UV limit is taken.
"""
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'finite_field_communication_results.json'
I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1.,-1.])
PARAM=dict(scalar_mass=1.,profile_sigma=1.,separation=2.,duration=30.,
           sender_coupling=.5,receiver_coupling=.5,position_std=1e-4,mass_lower_bound=1e14)

def norm(a):return float(np.linalg.norm(a,2))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def nodes(n,a,b):
    x,w=np.polynomial.legendre.leggauss(n)
    return a+(x+1)*(b-a)/2,w*(b-a)/2
def sinc(x):return np.sinc(x/np.pi)
def sincprime(x):return (x*np.cos(x)-np.sin(x))/x**2

def kernels(n,t=30.,r=2.):
    k,w=nodes(n,0.,12.);om=np.sqrt(k*k+1.)
    rad=w*k*k*np.exp(-k*k)/(2*np.pi**2);s=sinc(k*r)
    kp=float(np.sum(rad*s/(om*om)))
    theta=float(.25*np.sum(rad*s*(om*t-np.sin(om*t))/om**3))
    theta_r=float(.25*np.sum(rad*k*sincprime(k*r)*(om*t-np.sin(om*t))/om**3))
    n0=float(np.sum(rad*(1-np.cos(om*t))/om**3))
    nr=float(np.sum(rad*s*(1-np.cos(om*t))/om**3))
    noise=.25*np.array([[n0,nr],[nr,n0]])
    field0=float(np.sum(rad*(1-np.cos(om*t))/om**2))
    fieldr=float(np.sum(rad*s*(1-np.cos(om*t))/om**2))
    field=.25*np.array([[field0,fieldr],[fieldr,field0]])
    return dict(static_kernel=kp,phase=theta,phase_gradient=theta_r,noise=noise,field=field)

def run():
    p=PARAM;t=p['duration'];r=p['separation']
    a=kernels(256,t,r);b=kernels(384,t,r)
    quadrature_error=max(float(np.max(abs(np.asarray(a[k])-np.asarray(b[k])))) for k in a)
    assert quadrature_error<2e-11
    closed_k=(math.exp(-1)-math.exp(3)*math.erfc(2))/(16*np.pi)
    static_error=abs(a['static_kernel']-closed_k)
    assert static_error<1e-13
    k,w=nodes(256,0.,12.);u,wu=nodes(48,-1.,1.)
    om=np.sqrt(k*k+1.)[:,None]
    measure=(w*k*k/(4*np.pi**2))[:,None]*wu[None,:]
    lam_s=.5*np.exp(-k*k/2)[:,None]/np.sqrt(2*om)*np.ones((1,len(u)))
    lam_r=lam_s*np.exp(-1j*k[:,None]*r*u[None,:])
    sectors=np.array([[z,y] for z in (1.,-1.) for y in (-1.,1.)])
    alphas=[];phases=[];energies=[];forces=[]
    for z,y in sectors:
        lam=z*lam_s+y*lam_r
        alpha=lam*(np.exp(-1j*om*t)-1)/om
        phase=float(np.sum(measure*abs(lam)**2*(om*t-np.sin(om*t))/om**2))
        ef=float(np.sum(measure*om*abs(alpha)**2))
        ei=float(2*np.real(np.sum(measure*lam.conj()*alpha)))
        derivative_r=(-1j*k[:,None]*u[None,:])*y*lam_r
        force=float(-2*np.real(np.sum(measure*derivative_r.conj()*alpha)))
        alphas.append(alpha);phases.append(phase);energies.append([ef,ei]);forces.append(force)
    # Independent coherent-state overlaps versus the reduced Gaussian Schur channel.
    gram=np.empty((4,4),complex);closed=np.empty_like(gram)
    for i in range(4):
        for j in range(4):
            overlap=np.sum(measure*alphas[j].conj()*alphas[i])
            distance=float(np.sum(measure*abs(alphas[i]-alphas[j])**2))
            gram[i,j]=np.exp(1j*(phases[i]-phases[j]+overlap.imag)-distance/2)
            d=sectors[i]-sectors[j]
            closed[i,j]=np.exp(1j*a['phase']*(np.prod(sectors[i])-np.prod(sectors[j]))-.5*d@a['noise']@d)
    overlap_error=norm(gram-closed)
    assert overlap_error<1e-12
    positivity=float(np.linalg.eigvalsh(gram).min())
    assert positivity>-1e-12 and norm(np.diag(np.diag(gram))-np.eye(4))<1e-12
    # Reuse the actual invariant two-mode reference effect of round 939.
    y,vy=np.linalg.eigh(Y);rotate=np.kron(vy,vy)
    px=np.array([1,1],complex)/np.sqrt(2)
    prep=rotate.conj().T@np.kron(px,px)
    rho=np.outer(prep,prep.conj());charges=(y[:,None]+y[None,:]).reshape(-1)
    rho*=np.isclose(charges[:,None],charges[None,:])
    relational=rotate.conj().T@((np.kron(Z,X)-np.kron(X,Z))/2)@rotate
    read=(np.eye(4)+relational)/2
    probabilities=[]
    for i in (0,1):
        rho_s=np.zeros((2,2));rho_s[i,i]=1
        joint=np.kron(rho_s,rho).reshape(4,2,4,2)
        joint=(joint*gram[:,None,:,None]).reshape(8,8)
        probabilities.append(float(np.real(np.trace(joint@np.kron(I,read)))))
    exact_gap=.5*np.exp(-2*a['noise'][1,1])*np.sin(2*a['phase'])
    probability_error=abs(probabilities[0]-probabilities[1]-exact_gap)
    assert probability_error<1e-13 and exact_gap>0
    # Direct time integration of the microscopic force equals the phase gradient.
    times,tw=nodes(128,0.,t)
    rad=w*k*k*np.exp(-k*k)/(2*np.pi**2)
    force_time=.25*np.sum((rad*k*sincprime(k*r)/(k*k+1))[:,None]
                          *(1-np.cos(np.sqrt(k*k+1)[:,None]*times[None,:])),axis=0)
    impulse=float(force_time@tw)
    force_identity_error=abs(impulse-a['phase_gradient'])
    fd=(kernels(256,t,r+1e-4)['phase']-kernels(256,t,r-1e-4)['phase'])/(2e-4)
    derivative_error=abs(fd-a['phase_gradient'])
    assert force_identity_error<1e-12 and derivative_error<1e-9
    radial_force=.25*np.sum(rad*k*sincprime(k*r)*(1-np.cos(np.sqrt(k*k+1)*t))/(k*k+1))
    force_mode_error=max(abs(f-z*y*radial_force) for f,(z,y) in zip(forces,sectors))
    energy_error=float(np.max(abs(np.sum(energies,axis=1))))
    energy_formula_error=max(abs(en[0]-s@a['field']@s) for en,s in zip(energies,sectors))
    assert force_mode_error<1e-12 and energy_error<1e-12 and energy_formula_error<1e-12
    # Analytic, deliberately loose bounds: no quadrature convergence assumption is
    # used to assert a strictly positive signal in the finite-mass dynamical model.
    i3_upper=1/(8*np.pi**1.5)
    c_upper=3/(32*np.pi**1.5)
    k_lower=math.exp(-1)/(16*np.pi)*(1-1/(2*np.sqrt(np.pi)))
    k_upper=math.exp(-1)/(16*np.pi)
    theta_lower=.25*(t*k_lower-i3_upper)
    theta_upper=.25*(t*k_upper+i3_upper)
    assert 0<theta_lower<theta_upper<np.pi/4
    static_gap_lower=.5*math.exp(-i3_upper)*math.sin(2*theta_lower)
    a2_upper=2*i3_upper  # 4 (g_s^2+g_r^2) I_3, uniformly at all intermediate times.
    delta=p['position_std'];mass=p['mass_lower_bound']
    kinetic_error=t*2*np.sqrt(15)/(8*mass*delta**2)
    position_error=t*delta*np.sqrt(c_upper)*np.sqrt(1+4*a2_upper)
    unitary_error_bound=float(kinetic_error+position_error)
    moving_gap_lower=float(static_gap_lower-2*unitary_error_bound)
    assert moving_gap_lower>.03
    # Gaussian tail: all powers used above are bounded by (T+2) k^4 e^-k^2
    # for k>=12 (also covers sine/cosine spatial derivatives at r=2).
    lower=12.;tail0=np.sqrt(np.pi)*math.erfc(lower)/2
    tail2=lower*math.exp(-lower**2)/2+tail0/2
    tail4=lower**3*math.exp(-lower**2)/2+3*tail2/2
    tail_envelope=(t+2)*tail4/(2*np.pi**2)
    files=[Path(__file__),STAGE/'943/material_content_access_results.json',
           STAGE/'939/receiver_relational_access_results.json',
           STAGE/'853/dirac_receiver_record_results.json']
    # The old 853 result filename is discovered without silently substituting an input.
    if not files[-1].exists():
        candidates=sorted((STAGE/'853').glob('*results.json'))
        assert len(candidates)==1,candidates
        files[-1]=candidates[0]
    return dict(round=944,date='2026-10-07',all_scientific_checks_passed=True,parameters=p,
        finite_field=dict(phase=a['phase'],noise_matrix=a['noise'].tolist(),static_kernel=a['static_kernel'],
            relational_plus_probabilities=probabilities,relational_probability_gap=float(exact_gap),
            independent_CAR_effect_probability_error=probability_error,
            coherent_overlap_channel_error=overlap_error,channel_gram_minimum_eigenvalue=positivity,
            radial_quadrature_comparison_error=quadrature_error,static_closed_form_error=static_error,
            analytic_radial_tail_envelope=float(tail_envelope),quadrature_agreement_is_not_interval_certificate=True),
        same_source_resources=dict(sector_field_and_interaction_energies=energies,
            total_energy_balance_error=energy_error,independent_field_energy_error=energy_formula_error,
            phase_gradient=a['phase_gradient'],integrated_receiver_force=impulse,
            phase_gradient_equals_impulse_error=force_identity_error,
            force_from_displacement_vs_radial_error=force_mode_error,finite_difference_derivative_error=derivative_error,
            finite_mass_H_has_exact_total_momentum_conservation=True,
            fixed_center_force_is_not_automatic_full_moving_force_error_bound=True),
        finite_mass_control=dict(I3_upper=float(i3_upper),gradient_form_factor_upper=float(c_upper),
            static_kernel_lower=float(k_lower),static_kernel_upper=float(k_upper),
            phase_lower=float(theta_lower),phase_upper=float(theta_upper),
            fixed_center_relational_gap_lower=float(static_gap_lower),
            kinetic_Duhamel_part=float(kinetic_error),position_Duhamel_part=float(position_error),
            global_state_norm_error_bound=unitary_error_bound,
            actual_moving_model_probability_gap_lower=moving_gap_lower,
            bound_uniform_over_unknown_internal_input_and_passive_reference=True),
        scope=dict(new_heavy_particle_effective_realization_is_explicit_physical_input=True,
            profile_width_not_claimed_physical_minimum_scale=True,point_source_limit_taken=False,
            finite_mass_and_quantum_field_both_retained=True,
            Gaussian_profiles_not_claimed_compact_causal_support=True,
            bounded_probability_error_does_not_certify_unbounded_stress_error=True,
            original_Dirac_GR_SM_joint_matching_certified=False,
            autonomous_reference_preparation_and_readout_certified=False,
            full_six_protocol_organization_completed=False,full_goal_completed=False),
        source_hashes={str(q.relative_to(ROOT)):sha(q) for q in files})

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert abs(a-b)<1e-10,(a,b)
    else:assert a==b,(a,b)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,json.loads(TARGET.read_text('utf-8')))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
