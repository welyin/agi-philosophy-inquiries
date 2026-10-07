"""971: shared native mass gap, covariant dipole and recoil kinematics.
No scattering rate, real molecule selection rules or full common EFT is certified.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
OUT=HERE/"native_covariant_source_results.json"
ETA=np.diag([-1.,1.,1.,1.])
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mdot(a,b):return float(a@ETA@b)
def boost(v):
    v=np.asarray(v,float);v2=float(v@v)
    B=np.eye(4)
    if v2:
        g=1/math.sqrt(1-v2);B[0,0]=g;B[0,1:]=g*v;B[1:,0]=g*v
        B[1:,1:]+=(g-1)*np.outer(v,v)/v2
    assert np.max(abs(B.T@ETA@B-ETA))<1e-14
    return B
def run():
    spec=importlib.util.spec_from_file_location("old956",STAGE/"956/native_material_interface.py")
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    m=old.material(1.,.1);J=m["J"];d=m["occupancy"]
    cs=[old.annihilate(i) for i in range(4)];ns=[c.T@c for c in cs]
    sec=np.eye(16)[:,[s for s in range(16) if s.bit_count()==2]]
    Q=sec.T@(ns[0]+ns[1]-ns[2]-ns[3])@sec/2
    g=math.sqrt(1-d)*m["s"]+math.sqrt(d)*m["d"];e=Q@m["d"]
    vertex=float(e@Q@g);Delta=float(e@m["H"]@e-g@m["H"]@g)
    assert abs(vertex-math.sqrt(d))<1e-14 and abs(Delta-(1+J))<1e-14
    m0=100.;Mg=m0-J;Me=m0+1.
    # Fixed allowed-transition mass input for the kinematic audit.
    # Real rotational/structural energies would have to be included in these masses.
    wabs=Delta+Delta*Delta/(2*Mg)
    wem=Delta-Delta*Delta/(2*Me)
    shift=wabs-Delta
    assert Mg>0 and wem>0
    # Rest absorption P_e=P_g+k, and emission P_g=P_e-k.
    n=np.array([1.,0.,0.])
    ka=np.r_[wabs,wabs*n];ke=np.r_[wem,wem*n]
    Pg=np.array([Mg,0.,0.,0.]);Pe=Pg+ka
    Pe_rest=np.array([Me,0.,0.,0.]);Pg_final=Pe_rest-ke
    assert abs(mdot(ka,ka))<1e-13 and abs(mdot(ke,ke))<1e-13
    assert abs(mdot(Pe,Pe)+Me*Me)<5e-12
    assert abs(mdot(Pg_final,Pg_final)+Mg*Mg)<5e-12
    # Polarization tensor D=u wedge d, with d spacelike, u.d=0.
    u=np.array([1.,0.,0.,0.]);dip=np.array([0.,0.,0.,vertex])
    pol=np.array([0.,0.,0.,1.])
    def amplitude(U,D,K,eps):
        tensor=np.outer(U,D)-np.outer(D,U)
        kl=ETA@K;el=ETA@eps
        F=np.outer(kl,el)-np.outer(el,kl)
        full=float(np.sum(tensor*F)/2)
        # Deliberately omit spatial-spatial entries in a moving frame.
        truncated=tensor.copy();truncated[1:,1:]=0
        naive=float(np.sum(truncated*F)/2)
        current=kl@tensor
        ward=float(kl@current)
        return full,naive,ward
    base=amplitude(u,dip,ka,pol)[0];rows=[]
    max_ward=0.;max_cov=0.;max_gauge=0.;max_shell=0.
    for v in ([0.,0.,0.],[1/3,0.,0.],[.2,-.1,.15],[-.4,.1,0.]):
        B=boost(v)
        U,D,K,eps=(B@x for x in (u,dip,ka,pol))
        full,naive,ward=amplitude(U,D,K,eps)
        shifted=amplitude(U,D,K,eps+.7*K)[0]
        max_ward=max(max_ward,abs(ward))
        max_cov=max(max_cov,abs(full-base))
        max_gauge=max(max_gauge,abs(full-shifted))
        incoming=B@Pg;outgoing=B@Pe
        shell=abs(mdot(outgoing,outgoing)+Me*Me)
        max_shell=max(max_shell,shell)
        assert np.max(abs(outgoing-incoming-K))<2e-14
        comoving=-mdot(U,K)
        assert abs(comoving-wabs)<1e-14
        rows.append(dict(velocity=list(v),lab_photon_energy=float(K[0]),
            comoving_photon_energy=comoving,full_vertex=full,
            electric_only_vertex=naive,
            electric_only_relative_error=abs(naive/base-1),
            gauge_shift_error=abs(full-shifted),ward_residual=abs(ward),
            final_mass_shell_residual=shell))
    assert max(max_ward,max_cov,max_gauge)<1e-13 and max_shell<2e-11
    assert abs(rows[1]["electric_only_relative_error"]-.5)<1e-12
    # Rest vertices and currents are from the same operator, not independently fit.
    # 970 width is used ONLY as a concrete precision benchmark, not its vacuum momentum.
    old970=read(STAGE/"970/native_photon_scattering_results.json")
    delta=old970["parameters"]["half_packet_energy_width"]
    wrong_k=np.array([Delta,Delta,0.,0.]);wrong_final=Pg+wrong_k
    mismatch=mdot(wrong_final,wrong_final)+Me*Me
    # Me=Mg+Delta => the unshifted photon leaves a mass-shell deficit Delta^2.
    assert abs(mismatch-Delta*Delta)<5e-12
    assert shift/delta>100
    target_mass=Delta*Delta/(2*delta)
    # Bare geometric inertia m0 differs from Mg; using it changes the same finite recoil.
    wrong_recoil_with_bare_mass=Delta*Delta/(2*m0)
    # Matching audit: 960's resolved charge polarizability is not an additional
    # independent contact term when the SAME transitions are still explicit.
    alpha_resolved=2*d/Delta
    field=1e-4
    E0=float(np.linalg.eigvalsh(m["H"])[0])
    Em=float(np.linalg.eigvalsh(m["H"]-field*Q)[0])
    Ep=float(np.linalg.eigvalsh(m["H"]+field*Q)[0])
    alpha_direct=-(Ep+Em-2*E0)/(field*field)
    alpha_duplicated=-(Ep+Em-2*E0-alpha_resolved*field*field)/(field*field)
    assert abs(alpha_direct-alpha_resolved)<1e-7
    assert abs(alpha_duplicated-2*alpha_resolved)<1e-7
    files=[STAGE/"956/native_material_interface.py",
        STAGE/"970/native_photon_scattering_results.json",
        STAGE/"research_note_941.md",STAGE/"research_note_953.md",
        STAGE/"957/drafts/unified_operation_hypotheses_v0_2.md",
        HERE/"drafts/native_parent_decision.md"]
    return dict(round=971,all_scientific_checks_passed=True,
        fixed_native_inputs=dict(U=1.,v=.1,m0=m0,J=J,d=d,
            electronic_mass_gap=Delta,electric_dipole_matrix_element=vertex,
            ground_mass=Mg,excited_mass=Me),
        recoil=dict(absorption_frequency_rest=wabs,emission_frequency_rest=wem,
            absorption_shift=shift,emission_shift=Delta-wem,
            no_recoil_mass_shell_residual=mismatch,
            benchmark_970_half_width=delta,
            recoil_shift_over_benchmark_width=shift/delta,
            sufficient_ground_mass_for_recoil_below_width=target_mass,
            bare_mass_recoil_difference=shift-wrong_recoil_with_bare_mass,
            benchmark_is_not_a_waveguide_to_vacuum_identification=True),
        covariance=dict(rest_vertex=base,rows=rows,
            maximum_ward_residual=max_ward,maximum_lorentz_vertex_error=max_cov,
            maximum_gauge_shift_error=max_gauge,maximum_mass_shell_residual=max_shell),
        matching=dict(static_field_step=field,resolved_polarizability=alpha_resolved,
            actual_CAR_energy_curvature=alpha_direct,
            wrong_explicit_plus_full_contact_curvature=alpha_duplicated,
            residual_contact_for_this_finite_material=0.,
            self_polarization_of_965_removed=False,
            full_SM_polarizability_claimed=False),
        scope=dict(covariant_dipole_current_identity=True,
            same_internal_energy_in_mass_and_recoil=True,
            standard_physics_inputs_explicit=True,
            actual_transition_rates_certified=False,
            real_rotational_selection_rules_matched=False,
            full_SM_to_material_matching=False,
            full_quantum_positive_parent_constructed=False,
            round970_disproved=False,full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files},
        references=["https://arxiv.org/abs/1806.00234",
                    "https://arxiv.org/abs/hep-th/0409156",
                    "https://arxiv.org/abs/hep-th/0511133"])
def compare(a,b,path=""):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+"/"+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for k,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f"/{k}")
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-9,abs_tol=2e-11),(path,a,b)
    else:assert a==b,(path,a,b)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--write",action="store_true");a=p.parse_args()
    result=run()
    if a.write:
        with OUT.open("x",encoding="utf-8") as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else:compare(result,read(OUT))
    print(json.dumps({k:result[k] for k in ("round","all_scientific_checks_passed",
          "fixed_native_inputs","recoil","covariance")},ensure_ascii=False,indent=2))
