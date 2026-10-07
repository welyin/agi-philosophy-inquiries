"""957: audit the task-domain budget, with actual saved 955/956 constants.

No new cutoff, parent field solution, communication run, or precision scan.
The number-mode example separates hypotheses; it is not a new physics model.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/"task_domain_budget_results.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a,2))
def module(path):
    spec=importlib.util.spec_from_file_location("native956",path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def run():
    r=read(STAGE/"955/joint_effective_window_results.json")
    pars=read(STAGE/"946/portal_common_process_results.json")["parameters"]["inherited_gravity_and_moving_parameters"]
    T=r["parameters"]["duration"]
    m0=math.sqrt(r["parameters"]["scalar_masses_squared"][0])
    mass=pars["base_mass"];mr=pars["receiver_mass"];G=pars["G"];core=pars["soft_core"]
    delta=pars["position_std"];I0=1/(8*math.pi**1.5);hmax=math.pi
    vmax=G*mr*(mass+hmax)/core
    b=math.sqrt(I0/(2*m0));q2=2*I0/m0**2
    Hnorm=hmax+2*math.sqrt(15)/(8*mass*delta**2)+b+vmax
    assert abs(Hnorm-r["bounds"]["global_time_window"]/r["parameters"]["readout_half_window"])<1e-12
    # ||V psi|| <= sqrt(q2)*sqrt(<A>) + b, A=Tkin+Hfield.
    # Young at 1/2 gives ||A psi|| <= 2(||H psi||+||h||+||VN||+b+q2/2).
    Anorm=2*(Hnorm+hmax+vmax+b+q2/2)
    field_mean=r["bounds"]["field_energy_upper"]
    V_each=.5*math.sqrt(4*I0/m0**2*field_mean+I0/m0)
    rates={"internal_h":hmax,"source_kinetic":Anorm,"receiver_kinetic":Anorm,
        "field_H":Anorm,"source_field_interaction":V_each,
        "receiver_field_interaction":V_each,"soft_Newton":vmax}
    bounds955=dict(total_H_rms_upper=Hnorm,positive_kinetic_plus_field_rms_upper=Anorm,
        same_saved_field_energy_upper=field_mean,term_rms_upper=rates,
        term_integrated_action_upper={k:T*v for k,v in rates.items()},
        sum_integrated_action_upper=T*sum(rates.values()))
    assert all(math.isfinite(x) for x in bounds955["term_integrated_action_upper"].values())
    # Check the algebraic Young step independently over a broad positive range.
    xs=np.concatenate(([0.],np.geomspace(1e-8,1e5,200)))
    assert np.max(np.sqrt(q2)*np.sqrt(xs)+b-(.5*xs+q2/2+b))<1e-13
    # A positive unbounded interaction: H=(omega+g Z)N, omega>|g|.
    omega=1.;g=.2;M2=4.;rows=[]
    for n in (8,32,128,512):
        time=math.pi/(2*g*n)
        unrestricted_spin_overlap=abs(math.cos(g*n*time))
        full_distance=math.sqrt(max(0.,1-unrestricted_spin_overlap**2))
        p=M2/(n*n)
        overlap=1-p+p*math.cos(g*n*time)
        limited_distance=math.sqrt(max(0.,1-overlap**2))
        limit=g*time*math.sqrt(M2)
        assert full_distance>1-1e-13 and limited_distance<=limit+1e-13
        # Fixed first moment alone does not bound the chosen RMS action.
        mean_energy=1.
        first_moment_only_rms=g*math.sqrt(mean_energy*n)
        rows.append(dict(n=n,time=time,unrestricted_trace_distance=full_distance,
            second_moment_budget=M2,tail_probability=p,task_trace_distance=limited_distance,
            task_trace_distance_upper=limit,first_moment_only_interaction_rms=first_moment_only_rms))
    # Independent entangled finite-support checks of the SAME infinite-mode law.
    rng=np.random.default_rng(957);N=np.arange(9,dtype=float)
    energies=(omega+g*np.array([1.,-1.])[:,None])*N[None,:]
    Delta=g*np.array([1.,-1.])[:,None]*N[None,:]
    worst_margin=float("inf")
    for t in (.001,.07,.4):
        psi=rng.normal(size=(2,9,3))+1j*rng.normal(size=(2,9,3))
        psi/=np.linalg.norm(psi)
        # Remove the shared free oscillator evolution; ancilla is retained.
        changed=np.exp(-1j*t*Delta)[:,:,None]*psi
        second=float(np.sum(abs(psi)**2*N[None,:,None]**2))
        actual=np.sqrt(max(0.,1-abs(np.vdot(psi,changed))**2))
        upper=abs(g*t)*math.sqrt(second)
        assert actual<=upper+1e-12
        worst_margin=min(worst_margin,upper-actual)
        assert np.max(abs(np.sum(abs(changed)**2,axis=(0,2))-np.sum(abs(psi)**2,axis=(0,2))))<1e-14
    # 956 finite electronic material already meets the strong all-state version.
    dimer=module(STAGE/"956/native_material_interface.py").material(1.,.1)
    hopping=norm(-.1*dimer["hop"]);coulomb=norm(dimer["D"])
    duration=(math.pi/2)/.04
    bounds956=dict(hopping_operator_norm=hopping,coulomb_operator_norm=coulomb,
        duration=duration,all_state_integrated_action_upper=duration*(hopping+coulomb),
        all_six_dimensional_postmeasurement_states_included=True)
    assert abs(hopping-.2)<1e-12 and abs(coulomb-1)<1e-12
    files=[Path(__file__),STAGE/"955/joint_effective_window_results.json",
        STAGE/"946/portal_common_process_results.json",STAGE/"956/native_material_interface.py",
        STAGE/"956/native_material_interface_results.json",
        STAGE/"935/drafts/unified_operation_hypotheses_v0_1.md"]
    return dict(round=957,date="2026-10-07",all_scientific_checks_passed=True,
        actual_955_budget=bounds955,actual_956_budget=bounds956,
        separating_example=dict(omega=omega,g=g,positive_energy_coefficient=omega-abs(g),
            full_interaction_operator_norm="infinite",fixed_second_moment_rows=rows,
            entangled_reference_minimum_bound_margin=worst_margin,
            instrument_contract="spin-only instruments preserve oscillator N moments"),
        scope=dict(A1_D_is_declared_revision_not_old_A1_proved=True,
            all_955_H_terms_controlled_during_saved_free_run=True,
            arbitrary_955_terminal_instrument_moment_closure_certified=False,
            finite_956_full_material_strong_budget_preserved=True,
            Duhamel_controls_bounded_tasks_not_all_unbounded_sources=True,
            old_propagation_or_dimension_theorems_automatically_inherited=False,
            complete_native_parent_model_certified=False,full_goal_completed=False),
        references=["https://arxiv.org/abs/1712.10267","https://arxiv.org/abs/1812.07447",
            "https://doi.org/10.1007/s00220-025-05282-w"],
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--write",action="store_true");a=p.parse_args()
    if a.write:assert not TARGET.exists()
    out=run()
    if a.write:
        with TARGET.open("x",encoding="utf-8") as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write("\n")
    else:
        old=read(TARGET)
        assert old["source_hashes"]==out["source_hashes"] and old["scope"]==out["scope"]
        assert abs(old["actual_955_budget"]["sum_integrated_action_upper"]-
                   out["actual_955_budget"]["sum_integrated_action_upper"])<1e-10
    print(json.dumps({k:v for k,v in out.items() if k!="source_hashes"},ensure_ascii=False,indent=2))
