"""Match a real aligned, vector-like down mass row at the existing upper neutral pole."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
def read(name):
    return json.loads((HERE/name).read_bytes())
CFG = read("inputs.json")
def module(name):
    return runpy.run_path(str(HERE/CFG[name]))
JOINT = module("charged_neutral")
MATCH = module("low_matching")
CP = module("CP_readings")
TABLE = read(CFG["weak_point"])["table_IV"]
WEAK = read(CFG["weak_point"])["B_K_CP_weak_offsets_from_independent_inverse"]
DIRECT = JOINT["DIRECT"]
IDENTITY = np.eye(3, dtype=complex)


def project_shared(arm):
    n = len(arm)
    cov = np.array(TABLE["correlation"])[:n,:n]*np.outer(TABLE["sigmas"][:n],TABLE["sigmas"][:n])
    avg = np.array(TABLE["means"][:n])
    prec = np.linalg.inv(cov)
    norm = (arm@prec@avg)/(arm@prec@arm)
    diff = norm*arm-avg
    return dict(radius=math.sqrt(diff@prec@diff), common_s_over_c=norm, pulls=diff/np.sqrt(np.diag(cov)))


def absolute_table(arm, uc, prof):
    # With Kuc=(u/c)_MZ/(u/c)_lattice, the four predicted ratios are arm*n/Kuc.
    norm=uc/arm[3]
    center=np.array(TABLE["means"])
    cov=np.array(TABLE["correlation"])*np.outer(TABLE["sigmas"],TABLE["sigmas"])
    radius=lambda x: np.sqrt((x-center)@np.linalg.solve(cov,x-center))
    precision=arm@np.linalg.solve(cov,arm)
    best=prof["common_s_over_c"]
    allowed=None
    if prof["radius"]<=2:
        delta=math.sqrt((4-prof["radius"]**2)/precision)
        allowed=[norm/(best+delta),norm/(best-delta)]
    return dict(u_over_c_at_MZ=uc,absolute_at_Kuc_one=norm*arm,
                radius_at_Kuc_one=radius(norm*arm),Kuc_interval_for_radius_two=allowed,
                min_Kuc_relative_change_for_two=(max(0.,allowed[0]-1.,1.-allowed[1]) if allowed else None))

def charged_endpoint(count, step, weak=None, up_shift=0., g3_shift=0., top_shift=0., charm_shift=0.):
    # Both charged shooting and the neutral frame reconstruction must see the same spectrum.
    engine = module("charged_neutral")
    original_load = engine["module"]
    saved = {}
    def add_spectrum(flavour):
        flow = flavour["FLOW"]
        flow["low"]["yukawa"]["up"][0][0] = engine["HAD"]["UP"][0]+up_shift*engine["HAD"]["UP"][1]
        top=flow["low"]["yukawa"]["up"][2]
        top[0]+=top_shift*top[1]
        charm=flow["low"]["yukawa"]["up"][1]
        charm[0]+=charm_shift*charm[1]
        gauge=flow["physics"]["weak"]["g3"]
        gauge[0]+=g3_shift*gauge[1]
        mass = flow["fields"]["neutral_poles_GeV"][2]
        soft_mass = math.hypot(mass,flow["MS"])
        extra = [(mass,count*np.array(CFG["new_Dirac_triplet_beta"])),
                 (soft_mass,count*np.array(CFG["new_two_scalar_triplet_beta"]))]
        if count:
            flow["physics"]["thresholds"].extend(extra)
            flow["segments"].cache_clear()
        return flavour
    def load(which):
        result = original_load(which)
        if which=="flavour":
            return add_spectrum(result)
        if which=="neutral":
            original=result["load"]
            def same_spectrum(path):
                return add_spectrum(original(path))
            result["load"]=same_spectrum
            # shoot wraps the endpoint loader once more to install its chosen u mass.
            endpoint=result["endpoint"]
            def record(*args,**kwargs):
                obj=endpoint(*args,**kwargs)
                saved["endpoint"]=obj
                saved["neutral"]=result
                return obj
            result["endpoint"]=record
            endpoint.__globals__["load"]=same_spectrum
        return result
    engine["shoot"].__globals__["module"]=load
    summary=engine["shoot"](WEAK if weak is None else weak,up_shift,step)
    return summary,saved["endpoint"],saved["neutral"]


def initial_state(pt):
    f=pt["flow"]
    y=np.zeros((4,3,4),complex)
    y[0,:,:3]=np.diag(f["central"]("up"))
    y[1,:,:3]=f["low_ckm"](pt["row"]["gamma_MZ_deg"],f["angles"]("direct"))@np.diag(f["central"]("down"))
    y[2,:,:3]=pt["low_E"]@np.diag(f["central"]("charged_leptons"))
    return y


def advance(data, transport, dt, grid, flow):
    def rhs(y, s, g):
        e, nu=(part@part.conj().T for part in y[2:])
        return flow["beta"](y,g,False), -(-1.5*e+.5*nu)@s
    y,s=data,transport
    for j in range(0,len(grid)-1,2):
        a,k=rhs(y,s,grid[j])
        b,l=rhs(y+dt*a/2,s+dt*k/2,grid[j+1])
        c,m=rhs(y+dt*b/2,s+dt*l/2,grid[j+1])
        d,n=rhs(y+dt*c,s+dt*m,grid[j+2])
        y=y+dt*(a+2*b+2*c+d)/6
        s=s+dt*(k+2*l+2*m+n)/6
    return y,s


def ew_projection(y, mass, f):
    m=y*CFG["VEV_GeV"]/math.sqrt(2)
    mix=1-f*f
    gap=mass*mass-m*m
    disc=math.hypot(gap,2*m*mass*math.sqrt(mix))
    heavy=math.sqrt((m*m+mass*mass+disc)/2)
    return dict(heavy_GeV=heavy, light_over_original=f*mass/heavy,
                weak_left_missing=2*(m*mass)**2*mix/(disc*(disc+gap)),
                colored_Higgs_Yukawa_moment_GeV2=3*(mass*y)**2*mix/(8*math.pi**2))


def physical(summary, pt, nu, rule, step):
    f=pt["flow"]
    poles=f["fields"]["neutral_poles_GeV"]
    mass=poles[2]
    low=initial_state(pt)
    y=low.copy()
    sv=IDENTITY.copy()
    lower_segments=[segment for segment in f["segments"](step) if segment[0]<mass]
    for start,susy,dt,g in lower_segments:
        assert not susy
        for j,m in enumerate(poles):
            if start==m:
                y[3,:,j]=f["physics"]["Y0"][:,j]
        y,sv=advance(y,sv,dt,g,f)
    up,_=f["basis"](y[0])
    down,sd,right=np.linalg.svd(y[1,:,:3])
    down,sd,right=down[:,::-1],sd[::-1],right.conj().T[:,::-1]
    p=abs(up[:,0].conj()@down)**2
    c=math.sqrt(p[0])
    factors=np.array([1. if rule=="unit" else c, 1., {"unit":1.,"d_only":1.,"full_full":c,
                      "d_and_b_bisector":math.sqrt((1+c)/2)}[rule]])
    light=y.copy()
    light[1,:,:3]=y[1,:,:3]@(right@np.diag(factors)@right.conj().T)
    # Sum over the heavy subspace (the individual assignment of degenerate poles is free).
    missing_source_gram=(down@np.diag(sd**2*(1-factors**2))@down.conj().T)
    labelled=np.real(np.diag(down.conj().T@missing_source_gram@down))
    tagged=None
    if rule in ("full_full","d_and_b_bisector"):
        normalized=(labelled[2]/(sd[2]*factors[2])**2)/(labelled[0]/(sd[0]*factors[0])**2)
        expected=(c/(1+c))**2 if rule=="d_and_b_bisector" else 1.
        assert abs(normalized-expected)<1e-9
        tagged=dict(stripped_bH_over_dH=labelled[2]/labelled[0], divided_by_daughter_mass_squared=normalized)
    match_phase=np.angle(np.linalg.slogdet(light[1,:,:3])[0]/np.linalg.slogdet(y[1,:,:3])[0])
    after_u,after_s=f["basis"](light[1])
    basis_err=np.max(abs(abs(after_u.conj().T@down)-np.eye(3)))
    np.testing.assert_allclose(after_s/sd,factors,rtol=2e-12)
    sb=IDENTITY.copy()
    source_res=[]
    for start,_,dt,g in reversed(lower_segments):
        light,sb=advance(light,sb,-dt,g[::-1],f)
        for j,m in enumerate(poles):
            if start==m:
                old=f["physics"]["Y0"][:,j]
                if np.linalg.norm(old):
                    source_res.append(np.linalg.norm(light[3,:,j]-old)/np.linalg.norm(old))
                light[3,:,j]=0.
    q,lep,gamma,L,ups=f["extraction"](light)
    overlaps=np.diag(L.conj().T@pt["low_E"])
    L=L@(np.diag(overlaps/abs(overlaps)))
    H=pt["low_E"].conj().T@np.linalg.solve(sv,np.linalg.solve(sb,L))
    obs=nu["observables"](H.T@pt["C"]@H,H.T@pt["b"],nu["FROZEN"][0])
    solar=JOINT["read"](JOINT["CFG"]["reactor_and_solar"])
    atm=JOINT["read"](JOINT["CFG"]["atmospheric_angle"])["atmospheric_NOvA_upper_octant"]
    boxes=[solar["DayaBay_sin2_2theta13"],solar["JUNO_sin2_theta12_NO"],[atm["central"],atm["low"],atm["high"]]]
    obs["angle_marginal_widths"]=[JOINT["one_width"](v,*a) for v,a in zip(obs["angle_readings"],boxes)]
    obs["dm31_and_32_marginal_widths"]=[JOINT["one_width"](obs[key],*nu["FREQ"][ref])
                                         for key,ref in (("dm31_eV2","juno_NO_dm31"),("dm32_eV2","accelerator_NO_dm32"))]
    moduli=dict(zip(("Vus","Vub","Vcb"), np.sin(q[3:])*np.array([np.cos(q[4]),1,np.cos(q[4])])))
    weak={k[1:]:(val-DIRECT[k][0])/math.hypot(*DIRECT[k][1:]) for k,val in moduli.items()}
    ks,_=MATCH["thresholds"](dict(moduli_MZ=moduli,weak_offsets=weak,gamma=gamma),MATCH["source_thresholds"]())
    arm=np.array([q[2]/q[1]/ks[2],1,q[0]/q[1]/ks[0],ups[0]/q[1]/(ks[1]*ks[0])])
    sign=np.array([-1,1,-1])
    strength=lambda x: (p*sign)@x/(p[0]*sd[0])
    old_lep=JOINT["LEPTONS"]
    cp=CP["test"](moduli,gamma)
    # Short-distance sensitivity with published coefficient covariance: replace only the
    # explicit one-loop mass functions, preserving the central matching conventions.
    top_ref,charm_ref=163.48,1.27
    low_default=JOINT["REFERENCE"]["FLOW"]["low"]["yukawa"]["up"]
    short_factor=CP["loop"](top_ref*ups[2]/low_default[2][0],charm_ref*ups[1]/low_default[1][0],80.379)/CP["loop"](top_ref,charm_ref,80.379)
    eps_updated=np.dot(cp["tt_ut_parts"],short_factor)
    eps_sigma=math.hypot(cp["correlated_error_at_fixed_CKM"],CP["EPS_ERROR"])
    margin=eps_sigma-abs(eps_updated-CP["EPS_EXP"])
    cp["explicit_charm_top_function_update"]={"factors":short_factor,"epsilon":eps_updated,
             "pull_conditional_on_published_QCD":(eps_updated-CP["EPS_EXP"])/eps_sigma,
             "remaining_additive_shift_for_one_width_symmetric":margin}
    ax=moduli["Vus"]/math.sqrt(1-moduli["Vus"]**2-moduli["Vub"]**2)
    low_phase=sum(np.angle(np.linalg.slogdet(light[j,:,:3])[0]/np.linalg.slogdet(low[j,:,:3])[0]) for j in (0,1))
    return dict(c=c,factors=factors,threshold_GeV=mass,CKM_jump_magnitude=basis_err, down_above=sd,down_below=after_s,
                parent_return_normalized=strength(sd),daughter_return_normalized=strength(after_s),
                matching_arg_det_change=match_phase,low_arg_det_change=low_phase,original_high_full_arg_det=pt["strong"],
                high_yukawas=pt["rep"]["Y_down_at_10PeV"],high_x=summary["x"],high_c=summary["c"],
                tan_beta=summary["positive_soft_tan_beta"], neutral_only_moment_GeV2=pt["row"]["moment_GeV2"],
                MZ_down=q[:3],MZ_up=ups,down_mass_widths=(q[:3]-JOINT["REFERENCE"]["DOWN"][:,0])/JOINT["REFERENCE"]["DOWN"][:,1],
                MZ_leptons=lep,lepton_fractional_changes=lep/old_lep[:,0]-1,
                lepton_marginal_widths=(lep-old_lep[:,0])/old_lep[:,1],
                MZ_moduli=moduli,moduli_widths=weak,gamma_deg=gamma,axial= ax,
                kaon_marginal_widths={key:(ax-mu)/err for key,(mu,err) in JOINT["KAON"].items()},
                strong_low_factors_ds_ud_bs=ks,hadron_shape_B_1_D_U=arm,three_mass_profile=project_shared(arm[:3]),
                four_mass_profile=project_shared(arm),absolute_table_IV=absolute_table(arm,ups[0]/ups[1],project_shared(arm)),
                CP=cp,neutrinos=obs,heavy_flavor_tagged_Higgs_strength=tagged,
                lower_neutral_source_relative_residual=source_res,max_rel_lepton_transport=np.linalg.norm(H-IDENTITY,2),
                finite_vev_d_and_b=[ew_projection(sd[k],mass,factors[k]) for k in (0,2)])


def run():
    rows={}
    for count,rules in [(0,["unit"]),(1,["d_only"]),(2,CFG["physical_projection_rules"])]:
        summary,pt,nu=charged_endpoint(count,CFG["step"])
        for rule in rules:
            rows[f"{count}_{rule}"]=physical(summary,pt,nu,rule,CFG["step"])
    for name,inputs in CFG["nearby_measured_inputs_in_marginal_widths"].items():
        def solve_nearby(step):
            summary,pt,nu=charged_endpoint(2,step,inputs["weak"],inputs["up"],inputs["g3"],inputs["top"],inputs["charm"])
            return physical(summary,pt,nu,"d_and_b_bisector",step)
        rows[name]=solve_nearby(CFG["step"])
        if name=="absolute_and_oscillations":
            fine_joint=solve_nearby(CFG["step"]/2)
            joint=rows[name]
            assert joint["absolute_table_IV"]["radius_at_Kuc_one"]<1.03 and all(joint["neutrinos"]["angles_in_box"])
            assert all(joint["neutrinos"]["dm31_and_32_in_box_at_central_solar"])
            assert max(abs(joint["down_mass_widths"]))<2 and joint["neutrinos"]["SK_T2K_FC68_90"]==[False,True]
            assert abs(joint["CP"]["explicit_charm_top_function_update"]["pull_conditional_on_published_QCD"])<1
            assert abs(joint["CP"]["sin2beta"]-.710)<.011 and max(abs(v) for v in joint["kaon_marginal_widths"].values())<2
    summary,pt,nu=charged_endpoint(2,CFG["step"]/2)
    fine=physical(summary,pt,nu,"d_and_b_bisector",CFG["step"]/2)
    base=rows["2_d_and_b_bisector"]
    diffs={key:abs(base[key]-fine[key]) for key in ("gamma_deg","low_arg_det_change","tan_beta")}
    diffs["shape_radius"]=abs(base["four_mass_profile"]["radius"]-fine["four_mass_profile"]["radius"])
    diffs["dm31"]=abs(base["neutrinos"]["dm31_eV2"]-fine["neutrinos"]["dm31_eV2"])
    diffs["joint_absolute_radius"]=abs(joint["absolute_table_IV"]["radius_at_Kuc_one"]-fine_joint["absolute_table_IV"]["radius_at_Kuc_one"])
    diffs["joint_delta_dm31"]=abs(joint["neutrinos"]["dm31_eV2"]-fine_joint["neutrinos"]["dm31_eV2"])
    assert diffs["joint_absolute_radius"]<1e-5 and diffs["joint_delta_dm31"]<2e-9
    assert diffs["shape_radius"]<1e-5 and diffs["dm31"]<1e-9 and diffs["gamma_deg"]<1e-5
    for row in rows.values():
        assert abs(row["parent_return_normalized"])<2e-8 and row["CKM_jump_magnitude"]<1e-10
        assert abs(row["low_arg_det_change"])<1e-10 and abs(row["matching_arg_det_change"])<1e-10
    def tidy(x):
        if isinstance(x,dict): return {k:tidy(v) for k,v in x.items()}
        if isinstance(x,(list,tuple,np.ndarray)): return [tidy(v) for v in x]
        if isinstance(x,(np.floating,float)): return float(f"{x:.11g}")
        if isinstance(x,np.bool_): return bool(x)
        return x
    return tidy(dict(weak_offsets=WEAK,rows=rows,half_step=diffs))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    if args.check:
        assert result==read("results.json")
        print("1739 verified charged thresholds, CKM, leptons and lattice profiles")
    else:
        (HERE/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf8")

