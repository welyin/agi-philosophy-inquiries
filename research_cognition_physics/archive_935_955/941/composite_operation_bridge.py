"""941: shared mass/clock/source extension of the actual 925 internal process.

Composite binding, universal mass-energy coupling and background are inputs.
No claim that the 925 apparatus is already realized by the original field action.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json
import numpy as np

HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'composite_operation_bridge_results.json'
OLD=STAGE/'925/record_resource_exchange.py'
spec=importlib.util.spec_from_file_location('old925',OLD)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
def norm(a):return float(np.linalg.norm(a,2))
def function(h,f):
    e,v=np.linalg.eigh(h);return (v*f(e))@v.conj().T
def exp_h(h,t):return function(h,lambda e:np.exp(-1j*t*e))
def sqroot(h):return function(h,np.sqrt)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run():
    z=old.model();h=z['h'];h0=z['h0'];hint=z['hint'];d=z['d'];eye=np.eye(d)
    m=100.;mass=m*eye+h;eps,v=np.linalg.eigh(h);masses=m+eps
    T=np.pi/(2*z['g']);N=.97;phi=N-1
    initial=[];target=[];w=[]
    for q in range(3):
        src=z['idx'](0,q,3);dst=z['idx'](1,q,3+q)
        a=eye[:,src];b=eye[:,dst]
        initial.append(a);target.append(-1j*np.exp(-1j*T*h0[src,src])*b)
        w.extend([a,b])
    initial=np.column_stack(initial);target=np.column_stack(target);w=np.column_stack(w)
    assert norm(exp_h(h,T)@initial-target)<5e-14
    # Common redshift rescales the actual interaction as well as bare energies.
    Hrest=N*mass
    ur=exp_h(Hrest,T/N)*np.exp(1j*m*T)
    exact_record_error=norm(ur@initial-target)
    bad=m*N*eye+h+phi*h0
    ub=exp_h(bad,T/N)*np.exp(1j*m*T)
    accepted=target@target.conj().T
    pr=lambda u:float(np.trace(initial.conj().T@u.conj().T@accepted@u@initial).real/3)
    p_good=pr(ur);p_bad=pr(ub)
    bad_formula=np.sin(z['g']*T/N)**2
    missing_source=norm(mass-(m*eye+h0))
    assert exact_record_error<2e-12 and abs(p_bad-bad_formula)<1e-12
    assert 1-p_bad>1e-3 and abs(missing_source-.4)<1e-13
    # A common-momentum preparation has a controlled, not exact, common proper time.
    p=20.;Ebase=np.sqrt(p*p+m*m);rate=N*m/Ebase;t=T/rate
    evals=N*np.sqrt(p*p+masses*masses)
    Hfixed=(v*evals)@v.conj().T
    linear=N*Ebase*eye+rate*h
    # h>=0 here; f'' decreases as positive mass increases, uniformly on the spectrum.
    hmax=float(eps.max());gen_bound=N*p*p*hmax*hmax/(2*(p*p+m*m)**1.5)
    gen_error=norm(Hfixed-linear)
    Uf=exp_h(Hfixed,t)*np.exp(1j*N*Ebase*t)
    Ulin=exp_h(h,T)
    unitary_error=norm(Uf-Ulin)
    p_fixed=pr(Uf);channel_bound=min(1.,t*gen_bound)
    assert gen_error<=gen_bound+1e-12 and unitary_error<=t*gen_bound+1e-12
    assert abs(p_fixed-1)<=channel_bound+1e-12
    # Same velocity means different momenta for different mass-energy eigenstates.
    vel=.3;gamma=1/np.sqrt(1-vel*vel)
    E=gamma*mass;P=gamma*vel*mass;tboost=gamma*T;xboost=vel*tboost
    phase=E*tboost-P*xboost
    boost_error=norm(phase-mass*T)
    boost_u=exp_h(phase,1.)*np.exp(1j*m*T)
    boost_record=norm(boost_u@initial-target)
    shell_error=norm(E@E-P@P-mass@mass)
    # These are phase comparisons on a specified path, not construction of a material boost device.
    assert boost_error<2e-12 and boost_record<3e-12 and shell_error<2e-11
    # Same source under noncommuting detuning: solve the square-root derivative Sylvester equation.
    Dh=z['nref']@z['nref'];ee=np.sqrt(p*p+masses*masses)
    dmat=v.conj().T@Dh@v
    derivative=N*v@(((masses[:,None]+masses[None,:])/(ee[:,None]+ee[None,:]))*dmat)@v.conj().T
    delta=2e-4
    # 5-point central difference: independent square-root evaluations.
    def at(a):
        mm=mass+a*Dh
        return N*sqroot(p*p*eye+mm@mm)
    diff=(-at(2*delta)+8*at(delta)-8*at(-delta)+at(-2*delta))/(12*delta)
    deriv_error=norm(diff-derivative)
    sylv=norm(Hfixed@derivative+derivative@Hfixed-N*N*(mass@Dh+Dh@mass))
    naive=N*function(h,lambda e:(m+e)/np.sqrt(p*p+(m+e)**2))@Dh
    naive_gap=norm(naive-derivative);naive_nonhermitian=norm(naive-naive.conj().T)
    assert deriv_error<1e-6 and sylv<1e-8 and naive_gap>1e-5 and naive_nonhermitian>1e-5
    # Exact regrouping retains the whole pair, including interaction and source.
    he=w.conj().T@h@w;me=m*np.eye(6)+he
    He=N*sqroot(p*p*np.eye(6)+me@me)
    dh_small=w.conj().T@Dh@w
    e_small,s=np.linalg.eigh(he);ms=m+e_small;es=np.sqrt(p*p+ms*ms)
    de=N*s@(((ms[:,None]+ms[None,:])/(es[:,None]+es[None,:]))*(s.conj().T@dh_small@s))@s.conj().T
    regroup=dict(internal_intertwiner_error=norm(h@w-w@he),
        mass_shell_intertwiner_error=norm(Hfixed@w-w@He),
        source_intertwiner_error=norm(derivative@w-w@de))
    assert max(regroup.values())<2e-12
    # Finite joint energy/source check on every input column; no external mean-field substitution.
    initial_energy=np.diag(initial.conj().T@Hrest@initial).real
    final_energy=np.diag(initial.conj().T@ur.conj().T@Hrest@ur@initial).real
    energy_error=float(np.max(abs(initial_energy-final_energy)))
    source_difference_on_initial=np.diag(initial.conj().T@hint@initial).real
    assert energy_error<2e-12 and np.max(abs(source_difference_on_initial))==0
    paths=[OLD,STAGE/'925/record_resource_exchange_results.json',
           STAGE/'935/drafts/unified_operation_hypotheses_v0_1.md',
           STAGE/'939/drafts/common_model_recovery_map.md',Path(__file__)]
    return dict(round=941,date='2026-10-07',all_scientific_checks_passed=True,
        inherited_internal_dimension=d,retained_pair_dimension=6,base_mass=m,
        internal_energy_minimum=float(eps.min()),internal_energy_maximum=hmax,
        redshift=dict(N=N,proper_transfer_time=T,coordinate_transfer_time=T/N,
            all_unknown_input_transfer_error=exact_record_error,probability=p_good,
            omitted_interaction_source_norm=missing_source,
            omission_invisible_to_initial_diagonal_energy_means=True,
            bad_bare_mass_only_probability=p_bad,bad_probability_formula=bad_formula,
            finite_record_probability_gap=p_good-p_bad,energy_error=energy_error),
        fixed_momentum=dict(momentum=p,coordinate_time=t,proper_rate_linear=rate,
            generator_error=gen_error,generator_error_bound=gen_bound,
            full_unitary_error=unitary_error,unitary_bound=t*gen_bound,
            every_input_and_ancilla_trace_distance_bound=channel_bound,
            record_probability=p_fixed),
        common_velocity=dict(velocity=vel,mass_shell_error=shell_error,
            path_phase_error=boost_error,all_unknown_input_transfer_error=boost_record,
            comparison_includes_spatial_phase=True,material_boost_device_constructed=False),
        source=dict(noncommuting_h_source_commutator_norm=norm(h@Dh-Dh@h),
            frechet_vs_independent_five_point_error=deriv_error,sylvester_identity_error=sylv,
            scalar_chain_rule_error=naive_gap,scalar_chain_rule_nonhermitian_norm=naive_nonhermitian),
        exact_regrouping=regroup,
        shared_internal_process_clock_mass_and_source_bridge=True,
        mass_energy_and_quantum_equivalence_coupling_are_physical_inputs=True,
        original_field_material_binding_or_full_worldline_matching_proved=False,
        dynamical_gravitational_backreaction_or_full_SM_proved=False,
        new_clock_precision_optimization_performed=False,full_goal_completed=False,
        source_hashes={str(q.relative_to(ROOT)):sha(q) for q in paths})

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args();r=run()
    if a.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        before=json.loads(TARGET.read_text('utf-8'))
        def compare(x,y):
            if isinstance(x,dict):
                assert x.keys()==y.keys()
                for k in x:compare(x[k],y[k])
            elif isinstance(x,list):
                assert len(x)==len(y)
                for v,w in zip(x,y):compare(v,w)
            elif isinstance(x,float):assert abs(x-y)<1e-9
            else:assert x==y,(x,y)
        compare(r,before)
    print(json.dumps({k:v for k,v in r.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
