"""959: a Gauss-constrained electrostatic parent of the actual 958 process.

The full link Hilbert spaces are rotors. The finite diagnostic block contains
the entire declared neutral physical sector; no UV cutoff is asserted.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/"gauss_material_parent_results.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a,2))
def module(path):
    spec=importlib.util.spec_from_file_location("old_material",path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def comm(a,b):return a@b-b@a
def run():
    old=module(STAGE/"956/native_material_interface.py")
    saved=read(STAGE/"958/capacitive_material_write_results.json")
    U=saved["parameters"]["U"];v=saved["parameters"]["hopping"];k=saved["parameters"]["kappa"]
    m=old.material(U,v);h=m["H"];W=m["W"];d=m["occupancy"];J=m["J"]
    cs=[old.annihilate(i) for i in range(4)];ns=[c.T@c for c in cs]
    sector=np.eye(16)[:,[s for s in range(16) if s.bit_count()==2]]
    Q=sector.T@(ns[0]+ns[1]-ns[2]-ns[3])@sector/2
    tp=sector.T@sum(cs[i].T@cs[i+2] for i in (0,1))@sector
    assert norm(comm(Q,tp)-tp)<1e-14
    E=np.diag([-1.,0.,1.]);L=np.diag(np.ones(2),-1)
    assert norm(comm(E,L)-L)<1e-14
    Qloc=np.kron(Q,np.eye(3));Eloc=np.kron(np.eye(6),E)
    G=Eloc-Qloc
    Tloc=np.kron(tp,L);hop=-v*(Tloc+Tloc.T)
    j=np.zeros((18,6))
    for a,q in enumerate(np.diag(Q)):j[3*a+int(q)+1,a]=1.
    D=np.kron(j,j);P=D@D.T
    I=np.eye(18)
    EA=np.kron(Eloc,I);EB=np.kron(I,Eloc)
    QA=np.kron(Qloc,I);QB=np.kron(I,Qloc)
    GA=EA-QA;GB=EB-QB
    kinetic=np.kron(hop,I)+np.kron(I,hop)
    electric=U*(EA@EA+EB@EB)+k*EA@EB
    H=kinetic+electric
    oldHC=U*(np.kron(m["D"],np.eye(6))+np.kron(np.eye(6),m["D"]))+k*np.kron(Q,Q)
    oldH=np.kron(h,np.eye(6))+np.kron(np.eye(6),h)+k*np.kron(Q,Q)
    intertwining=norm(H@D-D@oldH)
    gauss_error=max(norm(GA@D),norm(GB@D),norm(comm(H,GA)),norm(comm(H,GB)))
    physical_count=int(np.count_nonzero((np.diag(GA)==0)&(np.diag(GB)==0)))
    assert physical_count==36 and intertwining<1e-13 and gauss_error<1e-13
    assert norm((np.eye(324)-P)@H@D)<1e-13
    K=np.array([[2*U,k],[k,2*U]]);C=np.linalg.inv(K)
    cg=1/(2*U+k);cm=k/(4*U*U-k*k)
    assert norm(C-np.array([[cg+cm,-cm],[-cm,cg+cm]]))<1e-14
    assert np.linalg.eigvalsh(K)[0]>0 and cg>0 and cm>0
    sources=dict(
        U=norm((EA@EA+EB@EB)@D-D@(np.kron(m["D"],np.eye(6))+np.kron(np.eye(6),m["D"]))),
        kappa=norm((EA@EB)@D-D@np.kron(Q,Q)),
        hopping=norm((kinetic/v)@D-D@(-np.kron(m["hop"],np.eye(6))-np.kron(np.eye(6),m["hop"]))))
    assert max(sources.values())<1e-13
    # Unknown state and arbitrary reference: a map identity, also sampled here.
    rng=np.random.default_rng(959)
    psi=rng.normal(size=(36,3))+1j*rng.normal(size=(36,3));psi/=np.linalg.norm(psi)
    out=D@psi
    assert abs(np.linalg.norm(out)-1)<1e-14
    # A local gauge-invariant, complete two-result physical instrument.
    lowplus=W@np.array([1.,1.,1.,1.])/2
    M=np.outer(lowplus,lowplus);Mg=j@M@j.T;Ploc=j@j.T
    instrument_error=norm(Mg@Mg+(Ploc-Mg)@(Ploc-Mg)-Ploc)
    assert instrument_error<1e-13 and norm(comm(Mg,G))<1e-13
    # Do NOT reinterpret tracing out constrained flux as the old matter state.
    singlet=np.array([0.,1.,-1.,0.])/np.sqrt(2);g=W@singlet
    charged=j@g;bare=charged.reshape(6,3)@charged.reshape(6,3).T
    old_hop=-v*m["hop"]
    expected_dressed=float(g@old_hop@g);bare_hop=float(np.trace(bare@old_hop))
    assert abs(charged@hop@charged-expected_dressed)<1e-14 and abs(bare_hop)<1e-14
    # The unrelated E=0 vacuum is not the Gauss dressing of the same unknown input.
    vac=np.array([0.,1.,0.]);bad_local=np.kron(g,vac);bad=np.kron(bad_local,bad_local)
    violation=float(np.vdot(bad,(GA@GA+GB@GB)@bad).real)
    admissible=float(np.linalg.norm(D.T@bad)**2)
    assert abs(violation-2*d)<1e-14 and abs(admissible-(1-d)**2)<1e-14
    # Actual 958 readout, original preparation and clock value; no retuning.
    b=np.eye(16);code=np.stack([(b[5]-b[6]-b[9]+b[10])/2,
        (2*b[3]+2*b[12]-b[5]-b[6]-b[9]-b[10])/np.sqrt(12)],axis=1)
    local=np.kron(W,np.eye(4))@code;record=np.kron(local,local)
    T=saved["exact_interface"]["controlled_phase_time"]
    target=(np.exp(1j*J*T)*local[:,0]+local[:,1])/np.sqrt(2);plus=np.array([1.,1.])/np.sqrt(2)
    def act(u,x):
        a=x.reshape(6,4,6,4,-1).transpose(0,2,1,3,4)
        a=(u@a.reshape(36,-1)).reshape(6,6,4,4,-1)
        return a.transpose(0,2,1,3,4).reshape(576,-1)
    def receiver(Hm):
        u=old.evolution(Hm,T);ps=[]
        for a in (0,1):
            initial=record@np.kron(np.eye(2)[:,a],plus)
            final=act(u,initial).reshape(24,24)
            ps.append(float(np.linalg.norm(final@target.conj())**2))
        return ps
    correct=receiver(oldH);doubled=receiver(oldH+oldHC)
    old_probs=saved["exact_interface"]["finite_time_rows"][1]["receiver_probabilities_for_logical_0_1"]
    assert max(abs(a-b) for a,b in zip(correct,old_probs))<1e-12
    double_difference=max(abs(a-b) for a,b in zip(correct,doubled))
    assert double_difference>.01
    # Full parent propagator check on its invariant physical sector at finite time.
    tcheck=.7;up=old.evolution(H,tcheck);um=old.evolution(oldH,tcheck)
    dynamic_error=norm(up@D-D@um);assert dynamic_error<1e-12
    files=[Path(__file__),STAGE/"956/native_material_interface.py",
        STAGE/"958/capacitive_material_write.py",STAGE/"958/capacitive_material_write_results.json",
        STAGE/"957/drafts/unified_operation_hypotheses_v0_2.md"]
    return dict(round=959,date="2026-10-07",all_scientific_checks_passed=True,
        parameters=dict(U=U,hopping=v,kappa=k,full_link_space="l2(Z) tensor l2(Z)",
            matrix_block_dimension=324,entire_declared_physical_dimension=physical_count,
            electric_quadratic_K=K.tolist(),capacitance_K_inverse=C.tolist(),
            ground_capacitance=cg,cross_capacitance=cm),
        exact_transport=dict(generator_error=intertwining,Gauss_error=gauss_error,
            finite_time_propagator_error=dynamic_error,all_source_errors=sources,
            instrument_completeness_error=instrument_error,
            inherits_958_unknown_input_upper=saved["exact_interface"]["all_time_unknown_input_error_upper"],
            inherits_958_read_window_contrast_lower=saved["exact_interface"]["analytic_receiver_contrast_lower"]),
        negative_controls=dict(unrelated_zero_flux_preparation=dict(
                sum_squared_Gauss=violation,probability_in_physical_sector=admissible),
            incorrect_bare_partial_trace=dict(material_purity=float(np.trace(bare@bare)),
                correct_single_material_hopping_energy=expected_dressed,
                bare_matter_hopping_energy=bare_hop),
            double_counted_electrostatics=dict(original_read_time=T,
                correct_receiver_probabilities=correct,duplicate_receiver_probabilities=doubled,
                max_probability_change=double_difference,extra_generator_norm=norm(oldHC))),
        scope=dict(exact_parent_for_958_static_constraint_sector=True,
            dressing_is_fixed_independent_of_unknown_input=True,
            no_duplicate_Coulomb_energy=True,physical_flux_bound_is_not_UV_cutoff=True,
            full_parent_SM_Einstein_certified=False,transverse_photon_recovery_certified=False,
            preparation_or_detector_manufactured=False,six_protocols_jointly_certified=False,
            three_dimensional_geometry_derived=False,full_goal_completed=False),
        references=["https://arxiv.org/abs/1807.01294","https://arxiv.org/abs/1409.3085",
            "https://doi.org/10.1103/PhysRevD.11.395"],
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})
if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true");args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    out=run()
    if args.write:
        with TARGET.open("x",encoding="utf-8") as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write("\n")
    else:
        saved=read(TARGET)
        assert saved["source_hashes"]==out["source_hashes"] and saved["scope"]==out["scope"]
        assert abs(saved["negative_controls"]["double_counted_electrostatics"]["max_probability_change"]-
            out["negative_controls"]["double_counted_electrostatics"]["max_probability_change"])<1e-10
    print(json.dumps({k:v for k,v in out.items() if k!="source_hashes"},ensure_ascii=False,indent=2))
