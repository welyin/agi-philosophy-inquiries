"""945: one scalar-field record process with a Newton constraint interaction.

The weak-field coupling and short-distance effective kernel are physical inputs.
This verifies joint bounded tasks, not full Einstein/SM matching or stress bounds.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'joint_record_newton_results.json'
I=np.eye(2,dtype=complex);X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1.,-1.])

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a,2))
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def run():
    old=load('field944',STAGE/'944/finite_field_communication.py')
    saved=json.loads((STAGE/'944/finite_field_communication_results.json').read_text('utf-8'))
    spectral=json.loads((STAGE/'929/joint_protocol_selection_results.json').read_text('utf-8'))
    assert saved['all_scientific_checks_passed'] and spectral['full_H_maximum']==30
    t=30.;radii=np.array([2.,3.]);mass=1e14;receiver_mass=1e14;newton=1e-30;core=1e-3;delta=1e-4
    u=1/np.sqrt(radii*radii+core*core)
    weights=np.array([math.comb(30,j)/2**30 for j in range(31)])
    masses=mass+np.arange(31);kappas=newton*receiver_mass*masses
    kappa_min=newton*receiver_mass*mass;kappa_max=float(kappas.max())
    k,w=old.nodes(256,0.,12.);q,wq=old.nodes(64,-1.,1.)
    om=np.sqrt(k*k+1.)[:,None]
    measure=(w*k*k/(4*np.pi**2))[:,None]*wq[None,:]
    ls=.5*np.exp(-k*k/2)[:,None]/np.sqrt(2*om)*np.ones((1,len(q)))
    sectors=[(l,b,y) for l in (0,1) for b in (1.,-1.) for y in (-1.,1.)]
    alpha=[];phase=[]
    for l,b,y in sectors:
        lr=ls*np.exp(-1j*k[:,None]*radii[l]*q[None,:]);lam=b*ls+y*lr
        alpha.append(lam*(np.exp(-1j*om*t)-1)/om)
        phase.append(float(np.sum(measure*abs(lam)**2*(om*t-np.sin(om*t))/om**2)))
    amplitude=np.asarray(alpha).reshape(8,-1)*np.sqrt(measure.reshape(-1))[None,:]
    overlap=amplitude@amplitude.conj().T
    squared=np.diag(overlap).real
    distance=squared[:,None]+squared[None,:]-2*overlap.real
    # overlap[i,j]=<alpha_j,alpha_i>, in the density-matrix convention used here.
    phase=np.asarray(phase)
    field_gram=np.exp(1j*(phase[:,None]-phase[None,:]+overlap.imag)-distance/2)
    gravity_gram=np.array([[np.sum(weights*np.exp(1j*kappas*t*(u[l]-u[ll])))
                           for ll,bb,yy in sectors] for l,b,y in sectors])
    joint_gram=field_gram*gravity_gram
    minimum=float(np.linalg.eigvalsh(joint_gram).min())
    assert minimum>-1e-12 and np.max(abs(np.diag(joint_gram)-1))<1e-12

    # Round 939 invariant relational reference; plus an explicitly declared
    # coherent path-label/recombination contract, not a claimed built apparatus.
    yy,vy=np.linalg.eigh(Y);rot=np.kron(vy,vy);plus=np.array([1,1],complex)/np.sqrt(2)
    prep=rot.conj().T@np.kron(plus,plus);rho_ra=np.outer(prep,prep.conj())
    total_y=(yy[:,None]+yy[None,:]).reshape(-1)
    rho_ra*=np.isclose(total_y[:,None],total_y[None,:])
    relative=rot.conj().T@((np.kron(Z,X)-np.kron(X,Z))/2)@rot
    read=(np.eye(4)+relative)/2
    path_plus=np.outer(plus,plus.conj());path0=np.diag([1.,0.]);path_read=(I-Y)/2
    effects=dict(path=np.kron(path_read,np.eye(8)),
        branch0_record=np.kron(path0,np.kron(I,read)),
        joint_path_record=np.kron(path_read,np.kron(I,read)))
    outputs={}
    for label,gram in [('G0',field_gram),('GN',joint_gram)]:
        rows=[]
        for b in (0,1):
            rho_s=np.zeros((2,2));rho_s[b,b]=1
            rho=np.kron(path_plus,np.kron(rho_s,rho_ra)).reshape(8,2,8,2)
            rho=(rho*gram[:,None,:,None]).reshape(16,16)
            assert abs(np.trace(rho)-1)<1e-12 and np.linalg.eigvalsh(rho).min()>-1e-12
            rows.append({name:float(np.real(np.trace(rho@effect))) for name,effect in effects.items()})
        outputs[label]=rows
    fields=[old.kernels(256,t,r) for r in radii]
    kernel_delta=old.kernels(256,t,float(radii[1]-radii[0]))
    dephasing=kernel_delta['noise'][1,1]-kernel_delta['noise'][0,1]
    field_visibility=float(np.exp(-dephasing)*np.cos(fields[0]['phase']-fields[1]['phase']))
    mass_coherence=np.sum(weights*np.exp(1j*kappas*t*(u[0]-u[1])))
    path_contrast=float(.5*field_visibility*mass_coherence.imag)
    actual_contrast=outputs['GN'][0]['path']-outputs['G0'][0]['path']
    path_error=abs(actual_contrast-path_contrast)
    record_gap=outputs['GN'][0]['branch0_record']-outputs['GN'][1]['branch0_record']
    record_error=abs(record_gap-.5*saved['finite_field']['relational_probability_gap'])
    assert path_error<1e-12 and record_error<1e-12
    assert actual_contrast>.024 and record_gap>.020
    assert abs(outputs['GN'][0]['joint_path_record']-outputs['G0'][0]['joint_path_record'])>.01

    # Uniform finite-position bound for the actual moving H, conditional on each
    # orthogonal path label.  The tail bound is deliberately rounded UP, not to zero.
    base_error=saved['finite_mass_control']['global_state_norm_error_bound']
    rmin=float(radii.min())
    assert rmin*rmin/(32*delta*delta)>100*np.log(10)
    tail_bound=2**1.5*1e-100
    shape_variance=24*delta**2/rmin**4+(1/core+1/rmin)**2*tail_bound
    gravity_position_error=float(t*kappa_max*np.sqrt(shape_variance))
    moving_error=base_error+gravity_position_error
    i3u=saved['finite_mass_control']['I3_upper']
    # |Theta(r)| <= g_s g_r (T K(0)+I3), with K(0)<=1/(8 pi^(3/2)).
    theta_abs_bound=.25*(t*i3u+i3u)
    visibility_lower=float(np.exp(-i3u)*np.cos(2*theta_abs_bound))
    phase_min=float(kappa_min*t*(u[0]-u[1]));phase_max=float(kappa_max*t*(u[0]-u[1]))
    assert 0<phase_min<=phase_max<np.pi/2 and 2*theta_abs_bound<np.pi/2
    static_path_contrast_lower=.5*visibility_lower*np.sin(phase_min)
    moving_path_contrast_lower=float(static_path_contrast_lower-moving_error-base_error)
    moving_record_gap_lower=float(.5*saved['finite_mass_control']['fixed_center_relational_gap_lower']-2*moving_error)
    assert moving_path_contrast_lower>.022 and moving_record_gap_lower>.015

    # Only the long-distance Newton term is transported from the Einstein input.
    point_phases=kappas[0]*t/radii;soft_phases=kappas[0]*t*u
    smoothing_bound=float(t*kappa_max*core**2/(2*rmin**3))
    assert float(np.max(abs(point_phases-soft_phases)))<=smoothing_bound*1.000001
    weak_potential_upper=float(newton*float(masses.max())/core)
    initial_rms_velocity=float(np.sqrt(3)/(2*mass*delta))
    # Recoil and source derivatives are taken from this same pair Hamiltonian.
    force=-kappas[0]*radii/(radii*radii+core*core)**1.5
    finite_difference=[]
    for r in radii:
        step=1e-5
        phase_plus=kappas[0]*t/np.sqrt((r+step)**2+core**2)
        phase_minus=kappas[0]*t/np.sqrt((r-step)**2+core**2)
        finite_difference.append((phase_plus-phase_minus)/(2*step))
    force_phase_error=float(np.max(abs(np.array(finite_difference)-t*force)))
    assert force_phase_error<1e-10
    # A changed G is one changed Hamiltonian, not a separately fitted phase.
    dg=1e-4
    def contrast_at(scale):
        z=np.sum(weights*np.exp(1j*scale*kappas*t*(u[0]-u[1])))
        return .5*field_visibility*z.imag
    response=(contrast_at(1+dg)-contrast_at(1-dg))/(2*dg)
    expected_response=.5*field_visibility*np.real(np.sum(weights*kappas*t*(u[0]-u[1])
                                                      *np.exp(1j*kappas*t*(u[0]-u[1]))))
    source_response_error=abs(response-expected_response)
    assert source_response_error<1e-10
    files=[Path(__file__),STAGE/'944/finite_field_communication.py',
        STAGE/'944/finite_field_communication_results.json',STAGE/'929/joint_protocol_selection_results.json',
        STAGE/'940/native_exchange_bridge_results.json',STAGE/'941/composite_operation_bridge_results.json']
    return dict(round=945,date='2026-10-07',all_scientific_checks_passed=True,
       parameters=dict(G=newton,base_mass=mass,receiver_mass=receiver_mass,internal_energy_range=[0,30],
            radii=radii.tolist(),duration=t,soft_core=core,position_std=delta),
       joint_process=dict(gram_minimum_eigenvalue=minimum,outputs=outputs,
            scalar_field_path_visibility=field_visibility,gravity_path_phase_minimum=phase_min,
            gravity_path_probability_contrast=actual_contrast,conditional_record_joint_event_gap=record_gap,
            inherited_record_probability_error=record_error,independent_path_formula_error=path_error,
            same_H_gravity_source_response_error=float(source_response_error)),
       finite_moving_control=dict(inherited_field_and_kinetic_bound=base_error,
            added_gravity_position_bound=gravity_position_error,full_state_error_bound=moving_error,
            conservative_field_visibility_lower=visibility_lower,
            gravity_path_probability_contrast_lower=moving_path_contrast_lower,
            conditional_record_joint_event_gap_lower=moving_record_gap_lower,
            uniform_unknown_input_and_old_reference=True,Gaussian_tail_probability_bound=tail_bound),
       leading_gravity_match=dict(long_distance_Newton_phase_softening_error_upper=smoothing_bound,
            global_dimensionless_potential_upper=weak_potential_upper,initial_rms_velocity=initial_rms_velocity,
            same_pair_energy_force_phase_derivative_error=force_phase_error,
            internal_mass_spectrum_retained=True,Newton_constraint_normalization_reused_from_940=True,
            full_scalar_stress_postNewton_radiative_or_loop_matching=False),
       scope=dict(coherent_path_preparation_and_recombination_are_declared_access_inputs=True,
            short_distance_kernel_is_effective_choice_not_universal_minimum_scale=True,
            joint_record_and_weak_gravity_are_one_H_not_two_output_tables=True,
            six_protocol_organization_and_full_SM_recovery_completed=False,
            full_goal_completed=False),
       source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

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
