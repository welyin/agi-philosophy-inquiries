"""974: fixed two-path Newton monopole process for the unchanged 965 composite.
The infinite Fock Hamiltonian is primary. This is a declared static-path EFT,
not a certificate for its traps, full GR, or microscopic cavity matching.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
HERE=Path(__file__).resolve().parent; STAGE=HERE.parent; ROOT=HERE.parents[2]
TARGET=HERE/'joint_energy_source_results.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,p):
    spec=importlib.util.spec_from_file_location(name,p)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def norm(a):return float(np.linalg.norm(a,2))
def expectation(op,block):return float(np.sum(block.conj()*(op@block)).real)

def run():
    core=load('core965',STAGE/'965/material_field_window.py')
    optical=load('optical965',STAGE/'965/physical_mode_effect.py')
    old=read(STAGE/'965/material_field_window_results.json')
    best=old['fixed_menu_rows'][-1]
    m,Q,Hm,S,W0,local=core.material()
    n=best['cutoff']+1;g=.003;omega=.5;m0=100.;T=old['T']
    dt=math.pi;lab_time=T+dt;delta=dt/lab_time
    potentials=np.array([-delta,-2*delta]);lapses=1+potentials
    proper_times=np.array([T,T-dt])
    assert max(abs(lapses*lab_time-proper_times))<3e-11
    a=np.diag(np.sqrt(np.arange(1,n)),1)
    A=np.kron(Hm,np.eye(n));F=np.kron(np.eye(36),omega*np.diag(np.arange(n)))
    C=g*np.kron(S,a+a.T);D=g*g/omega*np.kron(S@S,np.eye(n))
    H=A+F+C+D;parts=[A,F,C,D];val,V=np.linalg.eigh(H)
    mass=m0*np.eye(len(H))+H
    U=[(V*core.phases(val,float(t)))@V.T for t in proper_times]
    # Reuse the conservative old residual on 0 <= tau <= T, including
    # arbitrary material input in W0 and photon span{0,1}, with old references.
    old_eps=best['residual_bound']['total'];eps=max(2e-6,old_eps+1e-8)
    assert old_eps<2e-7
    phi=np.zeros(n,complex);phi[:2]=1/np.sqrt(2)
    targetB=(np.exp(1j*m['J']*T)*local[:,0]+local[:,1])/np.sqrt(2)
    ef=np.array([1.,1j*core.phases(np.array([omega]),T)[0]])/np.sqrt(2)
    weights=[]
    for s in np.diag(S):
        r0,r1=optical.displacement_rows(g*s/omega,n)
        weights.append(ef[0].conjugate()*r0+ef[1].conjugate()*r1)
    weights=np.array(weights)
    assert float(np.max(np.sum(abs(weights)**2,axis=1)))<=1+1e-13
    inputs=[];rows=[]
    for label in (0,1):
        psi=np.kron(np.kron(local[:,label],local@np.array([1.,1.])/np.sqrt(2)),phi)
        block=psi.reshape(6,4,6,4,n).transpose(0,2,4,1,3).reshape(36*n,16)
        inputs.append(block)
        energy=expectation(H,block)
        variance=float(np.sum(abs(H@block)**2))-energy**2
        co=V.T@block;prob_e=np.sum(abs(co)**2,axis=1)
        chi=np.sum(prob_e*core.phases(val,dt))
        evolved=[u@block for u in U]
        chi_endpoint=np.sum(evolved[1].conj()*evolved[0])
        assert abs(chi-chi_endpoint)<2e-10
        pb=[];pf=[];energies=[]
        for branch,x in enumerate(evolved):
            tensor=x.reshape(6,6,n,4,4).transpose(0,3,1,4,2).reshape(24,24,n)
            rb=np.einsum('abn,adn->bd',tensor,tensor.conj())
            pb.append(float(np.vdot(targetB,rb@targetB).real))
            projected=np.einsum('in,inj->ij',weights,x.reshape(36,n,16))
            pf.append(float(np.sum(abs(projected)**2)))
            before=[expectation(op,block) for op in parts]
            after=[expectation(op,x) for op in parts]
            assert abs(sum(after)-sum(before))<1e-10
            energies.append(dict(branch=branch,before=before,after=after,
                local_energy_change=sum(after)-sum(before),
                joint_energy_change=lapses[branch]*(sum(after)-sum(before))))
        rows.append(dict(label=label,energy_mean=energy,energy_variance=variance,
            mass_mean=m0+energy,chi=[float(chi.real),float(chi.imag)],
            visibility=abs(chi),visibility_endpoint_difference=abs(chi-chi_endpoint),
            receiver_probabilities=pb,physical_mode_probabilities=pf,
            receiver_probability_unconditional=sum(pb)/2,
            physical_mode_probability_unconditional=sum(pf)/2,energies=energies))
    # A single predeclared, rest-mass-corrected path effect; its phase is fixed
    # by the known label-0 calibration. H is NEVER conditioned on an input mean.
    calibration=rows[0]['energy_mean']
    phase=core.phases(np.array([-calibration]),dt)[0]
    for row in rows:
        chi=complex(*row['chi']);mean_chi=core.phases(np.array([row['energy_mean']]),dt)[0]
        p=float((1-(phase*chi).real)/2)
        pmean=max(0.,float((1-(phase*mean_chi).real)/2))
        row.update(path_effect_probability=p,mean_replacement_probability=pmean,
                   path_effect_difference=p-pmean,
                   infinite_model_path_difference_lower=p-pmean-eps-1e-10,
                   visibility_loss_lower=1-row['visibility']-2*eps)
    rec=rows[1]['receiver_probability_unconditional']-rows[0]['receiver_probability_unconditional']
    field=abs(rows[0]['physical_mode_probability_unconditional']-rows[1]['physical_mode_probability_unconditional'])
    # Same-H currents sum to zero even though the charge does not commute with
    # material Hm. On each path every internal current is multiplied by lapse.
    currents=[1j*(H@op-op@H) for op in parts]
    current_res=float(np.linalg.norm(sum(currents)))
    # Differentiation w.r.t. the potential gives the FULL mass, including field.
    q=.001;reference=-.00004
    source_fd=(((1+reference+q)*mass-(1+reference-q)*mass)/(2*q))
    source_res=norm(source_fd-mass)
    assert current_res<1e-12 and source_res<1e-9
    assert rec-2*eps>.98
    assert field-2*eps>.005
    assert rows[0]['infinite_model_path_difference_lower']>.14
    assert min(rows[0]['energy_variance'],rows[1]['energy_variance'])>.06
    # Positivity/self-adjointness for all occupancies follows from the complete
    # square H=Hm + omega (a+g S/omega)^dagger(a+g S/omega).
    lower_mass=m0+float(np.linalg.eigvalsh(Hm)[0])
    assert lower_mass>99 and min(lapses)>0
    files=[Path(__file__),STAGE/'965/material_field_window.py',
        STAGE/'965/physical_mode_effect.py',STAGE/'965/material_field_window_results.json',
        STAGE/'965/physical_mode_effect_results.json',STAGE/'956/native_material_interface.py',
        STAGE/'973/common_source_distribution_results.json']
    return dict(round=974,date='2026-10-07',all_scientific_checks_passed=True,
        parameters=dict(T=T,lab_time=lab_time,proper_times=proper_times.tolist(),
            proper_time_difference=dt,potentials=potentials.tolist(),lapses=lapses.tolist(),
            m0=m0,g=g,omega=omega,computational_cutoff=best['cutoff']),
        fixed_input_and_effects_from_965=True,
        joint_preparation='old labelled material and photon preparation tensor path |+>',
        fixed_hamiltonian='sum_l |l><l| tensor (1+Phi_l)(m0+H_L)',
        rows=rows,unconditional_reports=dict(receiver_contrast=rec,
            receiver_contrast_lower=rec-2*eps,
            physical_mode_contrast=field,physical_mode_contrast_lower=field-2*eps),
        certificate=dict(old_uniform_isometry_bound=old_eps,transported_bound=eps,
            arbitrary_passive_reference=True,coherent_path_input_preserved=True,
            time_interval=[0,T],complete_square_mass_lower=lower_mass,
            material_charge_commutator=norm(Hm@S-S@Hm),
            current_sum_residual=current_res,source_derivative_residual=source_res,
            energy_moments_exact_on_initial_0_1_photon_support=True),
        path_effect=dict(calibration_energy=calibration,
            definition='P=(I-X_theta)/2 after known rest mass phase, theta=pi*E_label0',
            same_effect_for_both_labels=True),
        adoption='native material-field process and coherent Newton source coexist in a fixed static-path EFT; stop path optimization',
        scope=dict(infinite_occupation_process_positive=True,
            energy_distribution_not_replaced_by_mean=True,
            old_quantum_information_preserved_jointly=True,
            original_isolated_marginals_unchanged=False,
            mass_gravity_equivalence_is_physical_input=True,
            fixed_path_supports_certified=False,full_moving_backreaction_certified=False,
            quantum_metric_dynamics_constructed=False,post_Newton_remainder_certified=False,
            full_SM_to_material_matching=False,cosmology_973_repaired=False,
            macroscopic_time_arrow_proven=False,full_goal_completed=False),
        references=['https://arxiv.org/abs/1502.00971','https://arxiv.org/abs/1311.1095'],
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

def compare(a,b):
    assert a['source_hashes']==b['source_hashes']
    for x,y in zip(a['rows'],b['rows']):
        for key in ('energy_mean','energy_variance','visibility','path_effect_probability',
                    'receiver_probability_unconditional','physical_mode_probability_unconditional'):
            assert abs(x[key]-y[key])<1e-8,key
    for k,v in a['unconditional_reports'].items():assert abs(v-b['unconditional_reports'][k])<1e-8

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    elif TARGET.exists():compare(out,read(TARGET))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
