"""958: same physical density interaction writes relational material content.

One fixed finite Hubbard candidate. No coupling scan, control compiler, or
claim of full QED/GR matching. Analytic proofs are in research_note_958.md.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/"capacitive_material_write_results.json"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text("utf-8-sig"))
def norm(a):return float(np.linalg.norm(a,2))
def load956():
    spec=importlib.util.spec_from_file_location("material956",STAGE/"956/native_material_interface.py")
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def comm(a,b):return a@b-b@a
def run():
    old=load956();U=1.;v=.1;k=.2
    m=old.material(U,v);h=m["H"];W=m["W"];J=m["J"];d=m["occupancy"]
    cs=[old.annihilate(i) for i in range(4)]
    ns=[c.T@c for c in cs]
    sector=np.eye(16)[:,[s for s in range(16) if s.bit_count()==2]]
    Q=sector.T@(ns[0]+ns[1]-ns[2]-ns[3])@sector/2
    assert norm(Q@Q-m["D"])<1e-14 and norm(W.T@Q@W)<1e-14
    HA=np.kron(h,np.eye(6));HB=np.kron(np.eye(6),h);V=np.kron(Q,Q)
    H0=HA+HB;H=H0+k*V;W0=np.kron(W,W)
    g=W@np.array([0.,1.,-1.,0.])/np.sqrt(2)
    psi0=np.kron(g,g);E0=-2*J
    ss=np.stack([m["s"],m["d"],Q@m["d"]],axis=1)
    B=np.kron(ss,ss);HS=B.T@H@B
    energies,vec=np.linalg.eigh(HS);E=float(energies[0]);phi=B@vec[:,0]
    if np.vdot(psi0,phi).real<0:phi=-phi
    chi=E0-E
    delta=U+J;eta=abs(k)*d/(delta-abs(k));error_upper=2*eta
    trial_lower=2*k*k*d*d/(math.sqrt((2*delta)**2+4*k*k*d*d)+2*delta)
    chi_upper=k*k*d*d/(delta-abs(k))
    chi2=k*k*d*d/(2*delta)
    assert k<delta and trial_lower<=chi<=chi_upper and chi>0
    difference=norm(phi-psi0)
    assert difference<=eta+1e-14
    sing=np.array([0.,1.,-1.,0.])/np.sqrt(2);joint_sing=np.kron(sing,sing)
    Wk=W0+np.outer(phi-psi0,joint_sing)
    K=np.kron(m["low"],np.eye(4))+np.kron(np.eye(4),m["low"])-chi*np.outer(joint_sing,joint_sing)
    identity_error=norm(H@Wk-Wk@K);iso_error=norm(Wk.T@Wk-np.eye(16))
    assert max(identity_error,iso_error)<2e-13
    # Exact full CAR/number-sector calculation independent of singlet reduction.
    ew,ev=np.linalg.eigh(H)
    evolve=lambda t:(ev*np.exp(-1j*t*ew))@ev.conj().T
    T=math.pi/chi;half_window=.01/(J+chi)
    rows=[]
    for time in (T-half_window,T,T+half_window):
        err=norm(evolve(time)@W0-W0@old.evolution(K,time))
        assert err<=error_upper+2e-9
        rows.append(dict(time=time,unknown_input_isometry_error=err))
    # Independent use of old 434 total-spin-zero relation code; no absolute axis.
    eye=np.eye(16)
    code=np.stack([(eye[5]-eye[6]-eye[9]+eye[10])/2,
        (2*eye[3]+2*eye[12]-eye[5]-eye[6]-eye[9]-eye[10])/np.sqrt(12)],axis=1)
    local=np.kron(W,np.eye(4))@code
    record=np.kron(local,local) # local ordering: active material, spectator pair
    def act(u,x):
        a=x.reshape(6,4,6,4,-1).transpose(0,2,1,3,4)
        a=(u@a.reshape(36,-1)).reshape(6,6,4,4,-1)
        return a.transpose(0,2,1,3,4).reshape(576,-1)
    assert norm(local.T@local-np.eye(2))<1e-14
    target=(np.exp(1j*J*T)*local[:,0]+local[:,1])/np.sqrt(2)
    plus=np.array([1.,1.])/np.sqrt(2)
    declared_contrast_lower=1-2*error_upper-.02
    energy_errors=[];source_values=[]
    for row in rows:
        u=evolve(row["time"]);probs=[]
        for a in (0,1):
            initial=record@np.kron(np.eye(2)[:,a],plus)
            evolved=act(u,initial).reshape(576)
            probs.append(float(np.linalg.norm(evolved.reshape(24,24)@target.conj())**2))
            def expect(op,x):return float(np.vdot(x,act(op,x).reshape(576)).real)
            energy_errors.append(abs(expect(H,evolved)-expect(H,initial)))
            if row["time"]==T and a==0:
                source_values.append(dict(actual_initial_code_pair=True,
                    charge_cross=expect(V,evolved),
                    charge_A=expect(np.kron(Q,np.eye(6)),evolved),
                    charge_B=expect(np.kron(np.eye(6),Q),evolved)))
        row["receiver_probabilities_for_logical_0_1"]=probs
        row["receiver_contrast"]=probs[1]-probs[0]
        assert row["receiver_contrast"]>=declared_contrast_lower
    # Same source operator, not a separately fitted force.
    derivative=float(phi@V@phi);dk=2e-5
    eval_at=lambda x:float(np.linalg.eigvalsh(B.T@(H0+x*V)@B)[0])
    fd=(eval_at(k+dk)-eval_at(k-dk))/(2*dk)
    source_error=norm(Wk.T@V@Wk-derivative*np.outer(joint_sing,joint_sing))
    assert abs(fd-derivative)<1e-9 and source_error<1e-13 and derivative<0
    assert abs(phi@np.kron(Q,np.eye(6))@phi)<1e-14
    assert abs(phi@np.kron(np.eye(6),Q)@phi)<1e-14
    currents=[1j*comm(H,A) for A in (HA,HB,k*V)]
    current_balance=norm(sum(currents));assert current_balance<1e-13
    assert max(energy_errors)<1e-12
    # Positive bare Coulomb density coefficients realize this operator up to
    # a state-independent term in the fixed particle sector; retain that term
    # for physical sources, or explicitly include neutralizing material.
    vc=.3;couplings=np.array([[vc+k/4,vc-k/4],[vc-k/4,vc+k/4]])
    nsites=[np.eye(6)+Q,np.eye(6)-Q]
    coulomb=sum(couplings[i,j]*np.kron(nsites[i],nsites[j]) for i in range(2) for j in range(2))
    density_error=norm(coulomb-(4*vc*np.eye(36)+k*V));assert density_error<1e-14
    # One static geometry dictionary, explicitly extra point-orbital input.
    R=3.;a=1.;coef=k/(2*(1/R-1/math.sqrt(R*R+a*a)))
    kprime=2*coef*(-1/R**2+R/(R*R+a*a)**1.5)
    ground_force=-kprime*derivative
    files=[Path(__file__),STAGE/"956/native_material_interface.py",
        STAGE/"956/native_material_interface_results.json",
        STAGE/"957/task_domain_budget_results.json",
        STAGE.parent/"archive_429_466/434/encoded_exchange_response_audit.py"]
    return dict(round=958,date="2026-10-07",all_scientific_checks_passed=True,
        parameters=dict(U=U,hopping=v,kappa=k,active_number_sector_dimension=36,
            relational_material_with_spectator_dimension=576,
            physical_input="fixed Hubbard plus electrostatic density interaction"),
        exact_interface=dict(one_body_charge_projection=norm(W.T@Q@W),
            double_occupancy_probability=d,first_order_joint_projection=norm(W0.T@V@W0),
            conditional_energy_chi=chi,second_order_chi=chi2,
            variational_chi_lower=trial_lower,resolvent_chi_upper=chi_upper,
            singlet_sector_gap_at_zero=delta,dressed_product_distance=difference,
            distance_upper_eta=eta,all_time_unknown_input_error_upper=error_upper,
            controlled_phase_time=T,read_window_half_width=half_window,
            analytic_receiver_contrast_lower=declared_contrast_lower,finite_time_rows=rows),
        sources=dict(ground_energy_kappa_derivative=derivative,
            finite_difference_derivative=fd,ground_marginal_charges=[0.,0.],
            actual_initial_code_at_read_time=source_values,
            extra_static_geometry=dict(R=R,a=a,coulomb_coefficient=coef,kappa_prime=kprime,
                conditional_ground_force=ground_force,
                moving_geometry_or_full_stress_certified=False)),
        verification=dict(exact_intertwining_error=identity_error,isometry_error=iso_error,
            same_source_intertwining_error=source_error,current_balance_error=current_balance,
            energy_conservation_max_error=max(energy_errors),
            positive_density_representation_error=density_error,
            all_density_coefficients_positive=bool(np.min(couplings)>0)),
        scope=dict(same_electrostatic_interaction_writes_relational_content=True,
            arbitrary_unknown_input_and_passive_reference_bound=True,
            controller_pulse_schedule_used=False,
            mean_charge_product_would_miss_cross_source=True,
            physical_preparation_and_terminal_readout_are_declared_inputs=True,
            full_six_protocols_implemented=False,full_QED_GR_matching_certified=False,
            physical_lifetime_at_calculated_read_time_certified=False,
            three_dimensions_derived=False,full_goal_completed=False),
        references=["https://arxiv.org/abs/1408.4740","https://arxiv.org/abs/1404.5420",
            "https://arxiv.org/abs/1105.0675"],
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})
if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true");args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    out=run()
    if args.write:
        with TARGET.open("x",encoding="utf-8") as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write("\n")
    else:
        saved=read(TARGET);assert saved["source_hashes"]==out["source_hashes"] and saved["scope"]==out["scope"]
        assert abs(saved["exact_interface"]["conditional_energy_chi"]-out["exact_interface"]["conditional_energy_chi"])<1e-14
    print(json.dumps({k:v for k,v in out.items() if k!="source_hashes"},ensure_ascii=False,indent=2))
