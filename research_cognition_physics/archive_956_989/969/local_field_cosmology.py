"""969: one explicit monopole/mean-field extension of the frozen 964/965 tasks.
The bound local mode is not relabelled as a free comoving radiation mode.
No real cavity matching, full quantum gravity or thermodynamic arrow is claimed.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
OUT=HERE/"local_field_cosmology_results.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,p):
    s=importlib.util.spec_from_file_location(name,p)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def comm(a,b):return a@b-b@a
def run():
    core=load("core965",STAGE/"965/material_field_window.py")
    m,Q,Hm,S,W0,local=core.material()
    old=read(STAGE/"965/material_field_window_results.json")
    phys=read(STAGE/"965/physical_mode_effect_results.json")
    cosmos=read(STAGE/"964/material_cosmology_results.json")["cosmology"]
    T=old["T"];J=m["J"];d=m["occupancy"];g=.003;omega=.5;m0=100.
    # SAME 964 gravitational coefficients and free-radiation preparation.
    mu=cosmos["mu"];R=cosmos["photon_comoving_energy"]
    V0=cosmos["coordinate_volume"];Mp2=cosmos["reduced_planck_mass_squared"]
    M=m0-J+omega/2+g*g/omega*d
    M_old=m0-J
    sources=old["fixed_menu_rows"][0]["sources"]["energy_rows"]
    initial=np.mean([r["before"] for r in sources],axis=0)
    final=np.mean([r["after"] for r in sources],axis=0)
    assert abs(m0+sum(initial)-M)<1e-12
    assert abs(sum(initial)-sum(final))<1e-11
    assert abs(initial[3]-g*g/omega*d)<1e-15
    assert abs(mu-3*Mp2*V0)<1e-3  # very large absolute coefficients
    # Nonselected quantum preparation has ONE background.
    def primitive(a,mass):
        return 2*(mass*a-2*R)*math.sqrt(mass*a+R)/(3*mass*mass)
    def time(a,mass=M):
        return math.sqrt(mu)*(primitive(a,mass)-primitive(1,mass))
    def scale(t,mass=M):
        lo,hi=1.,4.
        while time(hi,mass)<t:hi*=2
        for _ in range(65):
            mid=(lo+hi)/2
            if time(mid,mass)<t:lo=mid
            else:hi=mid
        return (lo+hi)/2
    def eta(a):
        return 2*math.sqrt(mu)/M*(math.sqrt(M*a+R)-math.sqrt(M+R))
    af=scale(T)
    assert abs(scale(T,M_old)-2)<1e-12
    assert af>2
    # Independent Raychaudhuri evolution, no constraint reset.
    def rhs(y):
        a,B,conformal,work=y
        return np.array([a*B,-3*T*T/(2*mu)*(M/a**3+4*R/(3*a**4)),
                         1/a,R*B/a])
    def solve(n):
        y=np.array([1.,T*math.sqrt((M+R)/mu),0.,0.])
        step=1/n;constraint=0.
        for _ in range(n):
            a=rhs(y);b=rhs(y+step*a/2);c=rhs(y+step*b/2);e=rhs(y+step*c)
            y+=step*(a+2*b+2*c+e)/6
            constraint=max(constraint,abs(y[1]**2-T*T*(M/y[0]**3+R/y[0]**4)/mu))
        return y,constraint
    yc,ec=solve(600);yf,ef=solve(1200)
    exact=np.array([af,T*math.sqrt((M*af+R)/(mu*af**4)),eta(af)/T,R*(1-1/af)])
    error=float(np.max(abs(yf-exact)))
    assert error<2e-9 and ef<2e-10
    rows=[]
    for a in (1.,1.5,af):
        hub=math.sqrt((M*a+R)/(mu*a**4))
        rho=(M/a**3+R/a**4)/V0;P=R/(3*V0*a**4)
        continuity=-hub*(3*M/a**3+4*R/a**4)/V0+3*hub*(rho+P)
        pa=-2*mu*a*a*hub
        constraint=-pa*pa/(4*mu*a)+M+R/a
        # For H_L independent of a, ONLY the free radiation contributes volume work.
        rows.append(dict(a=a,time_over_T=time(a)/T,
            free_photon_to_bare_material_gap_ratio=1/a,
            free_to_local_mode_frequency_ratio=(J/a)/omega,
            photon_temperature=cosmos["initial_temperature"]/a,
            local_complete_mass=M,total_comoving_energy=M+R/a,
            pressure_work=R*(1-1/a),
            continuity_residual=abs(continuity),constraint_residual=abs(constraint)))
        assert abs(continuity)<1e-15 and abs(constraint)<1e-10
    # A genuine finite interaction: local photon energy is NOT separately conserved.
    # The first two derivatives at t=0 have finite exact support in Fock number.
    derivatives=[]
    for n in (4,5):
        ann=np.diag(np.sqrt(np.arange(1,n)),1);N=np.diag(np.arange(n,dtype=float))
        A=np.kron(Hm,np.eye(n));F=np.kron(np.eye(36),omega*N)
        C=g*np.kron(S,ann+ann.T);D=g*g/omega*np.kron(S@S,np.eye(n))
        H=A+F+C+D
        current=1j*comm(H,F);acc=-comm(H,comm(H,F))
        phi=np.zeros(n);phi[:2]=1/np.sqrt(2)
        inputs=[]
        for label in (0,1):
            x=np.kron(np.kron(local[:,label],local@np.array([1.,1.])/np.sqrt(2)),phi)
            inputs.append(x.reshape(6,4,6,4,n).transpose(0,2,4,1,3).reshape(36*n,16))
        def exp(op):
            return float(sum(np.sum(x.conj()*(op@x)).real for x in inputs)/2)
        first=exp(current);second=exp(acc)
        assert abs(first)<1e-15 and abs(second-2*g*g*omega*d)<1e-15
        derivatives.append(dict(fock_dimension=n,first=first,second=second,
            current_operator_norm=float(np.linalg.norm(current,2))))
    H0=math.sqrt((M+R)/mu)
    # Wrong assignment: E_wrong=m0+<H_L-F>+<F>/a+R/a,
    # P_wrong=(<F>+R)/(3 V0 a^4) with the ORIGINAL bound H_L evolution.
    # Then dE_wrong/dt + P_wrong*dV/dt=(1/a-1)*d<F>/dt.
    # Leading t^2 coefficient is strictly negative, analytically.
    wrong_leading=-H0*2*g*g*omega*d
    assert wrong_leading<0
    wrong_mass_offset_at_T=(1/af-1)*float(final[1])
    assert wrong_mass_offset_at_T<-.12
    lower_record=old["fixed_menu_rows"][0]["probabilities"]["receiver_contrast"]-4e-6
    lower_field=phys["rows"][0]["conservative_reported_lower"]
    assert lower_record>.9987 and lower_field>.0498
    # Exact thermal entropy for the 12 unchanged free modes.
    nbar=cosmos["mean_occupation"];modes=cosmos["photon_modes"]
    sg=modes*((1+nbar)*math.log1p(nbar)-nbar*math.log(nbar))
    full_entropy=math.log(2)+sg
    files=[STAGE/"965/material_field_window.py",STAGE/"965/material_field_window_results.json",
           STAGE/"965/physical_mode_effect.py",STAGE/"965/physical_mode_effect_results.json",
           STAGE/"964/material_cosmology_results.json",STAGE/"research_note_968.md",
           HERE/"drafts/common_background_decision.md"]
    return dict(round=969,all_scientific_checks_passed=True,
        preparation="old 965 balanced material labels and coherent local mode, tensor old 964 free Gibbs modes",
        added_input="stable local support; disjoint bound and free mode menus; common monopole mean-field",
        original_material_parameters=old["old_material_parameters"],
        same_local_interaction_parameters=dict(g=g,omega=omega),
        shared_background_parameters=dict(mu=mu,Mp2=Mp2,V0=V0,R=R,T=T,m0=m0),
        local_sources=dict(mean_initial_parts=initial.tolist(),mean_final_parts=final.tolist(),
            complete_mass=M,previous_964_mass=M_old,
            mass_increment=M-M_old,analytic_mass_formula="-J + omega/2 + (g*g/omega)*d + m0",
            minimum_mass_lower=m0+float(np.linalg.eigvalsh(Hm)[0]),
            energy_change_error=float(abs(sum(final)-sum(initial)))),
        background=dict(rows=rows,scale_at_T=af,previous_scale_at_T=2.,
            exact_endpoint=exact.tolist(),integrated_endpoint=yf.tolist(),
            integration_error=error,step_difference=float(np.max(abs(yc-yf))),
            Raychaudhuri_constraint_max=ef),
        joint_task_transport=dict(exact_factorization=True,original_local_time=T,
            inherited_receiver_contrast_lower=lower_record,
            inherited_polarization_effect_contrast_lower=lower_field,
            inherited_local_residual_bound=old["fixed_menu_rows"][0]["residual_bound"]["total"],
            infinite_occupation_primary=True,
            photon_temperature_ratio=1/af,
            global_entropy_initial=full_entropy,global_entropy_final=full_entropy),
        wrong_dictionary=dict(energy_balance_identity="(1/a-1)*d<F>/dt",
            initial_field_derivatives=derivatives,nonzero_t_squared_coefficient=wrong_leading,
            source_energy_offset_at_T=wrong_mass_offset_at_T,
            applies_to="bound-mode dynamics retained but its energy assigned a free-radiation redshift"),
        scope=dict(one_nonselected_mean_background=True,
            microscopic_cavity_matching=False,microscopic_Einstein_derivation=False,
            full_quantum_geometry=False,external_radiation_absorption_recovered=False,
            irreversible_arrow_derived=False,round968_exchange_automatically_included=False,
            full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files},
        references=["https://arxiv.org/abs/0810.2712","https://arxiv.org/abs/1911.08427"])
def compare(a,b,path=""):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+"/"+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for k,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f"/{k}")
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-8,abs_tol=2e-11),(path,a,b)
    else:assert a==b,(path,a,b)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--write",action="store_true");a=p.parse_args()
    result=run()
    if a.write:
        with OUT.open("x",encoding="utf-8") as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else:compare(result,read(OUT))
    print(json.dumps({k:result[k] for k in ("round","all_scientific_checks_passed","local_sources",
        "background","joint_task_transport","wrong_dictionary")},ensure_ascii=False,indent=2))
