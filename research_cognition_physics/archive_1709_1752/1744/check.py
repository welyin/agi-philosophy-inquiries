"""Separate electroweak and color five-bar thresholds with two active neutral mass sources."""
import argparse
import json
import math
from pathlib import Path
import runpy
import numpy as np


HERE=Path(__file__).resolve().parent
def read(path):
    return json.loads((HERE/path).read_bytes())


CFG=read("inputs.json")
I=np.eye(3,dtype=complex)


def svd(y):
    l,s,rh=np.linalg.svd(y[:,:3])
    return l[:,::-1],s[::-1],rh.conj().T[:,::-1]


def mother(step=.125):
    model=runpy.run_path(str(HERE/CFG["parent"]))
    Mc,Mw=CFG["pole_mass_GeV"]
    gauge=read(CFG["new_weak_mass_from"])["complete_multiplets_on_previous_gauge_witness"]
    assert abs(gauge["weak_dirac_mass_repair_GeV"]-Mw)<1e-8
    # The frozen 1739 loader installs the two colored Dirac fields and their
    # own soft partners. Add the weak fields through both independent loaders.
    original=model["module"]
    def augment(flavor):
        flow=flavor["FLOW"]
        soft=math.hypot(flow["MS"],Mw)
        assert abs(soft-gauge["weak_scalar_mass_repair_GeV"])<1e-7
        flow["physics"]["thresholds"].extend([(Mw,np.array(CFG["two_weak_dirac_beta"])),
                                             (soft,np.array(CFG["two_weak_scalar_pair_beta"]))])
        flow["segments"].cache_clear()
        return flavor
    def module(which):
        obj=original(which)
        if which=="charged_neutral":
            orig=obj["module"]
            def sub(which):
                mod=orig(which)
                if which=="flavour":
                    return augment(mod)
                if which=="neutral":
                    load=mod["load"]
                    def with_weak(path):
                        return augment(load(path))
                    mod["load"]=with_weak
                    mod["endpoint"].__globals__["load"]=with_weak
                return mod
            obj["module"]=sub
            obj["shoot"].__globals__["module"]=sub
        return obj
    model["charged_endpoint"].__globals__["module"]=module
    data=model["CFG"]["nearby_measured_inputs_in_marginal_widths"][CFG["parent_point"]]
    row,pt,nu=model["charged_endpoint"](2,step,*(data[k] for k in ("weak","up","g3","top","charm")))
    flow=pt["flow"]
    poles=flow["fields"]["neutral_poles_GeV"]
    assert poles[2]==Mc
    all_segments=flow["segments"](step)
    required=[Mc,Mw,math.hypot(Mc,flow["MS"]),math.hypot(Mw,flow["MS"])]
    assert all(any(start==pole for start,_,_,_ in all_segments) for pole in required)
    segments=[part for part in all_segments if part[0]<Mw]
    y=model["initial_state"](pt)
    S=I.copy()
    for start,susy,dt,g in segments:
        assert not susy
        for j,m in enumerate(poles):
            if start==m:
                y[3,:,j]=flow["physics"]["Y0"][:,j]
        y,S=model["advance"](y,S,dt,g,flow)
    # Recover the constant right field identification using the entire high
    # complex D and CKM phase convention, as in 1712/tagged_matrices.
    high=flow["shoot"](pt["row"]["gamma_MZ_deg"],flow["angles"]("direct"),pt["low_E"],step)[0]
    Uh,_,_=svd(high[0])
    Lh,_,Rh=svd(high[1])
    Lt,_,Rt=svd(pt["D"])
    phases=(Uh.conj().T@Lh)/Lt
    np.testing.assert_allclose(abs(phases),1.,atol=3e-8)
    phases/=abs(phases)
    left=Uh*(phases[:,0]/phases[0,0])
    right=(Rh*phases[0,:].conj())@Rt.conj().T
    np.testing.assert_allclose(left.conj().T@high[1,:,:3]@right,pt["D"],atol=3e-10)
    Dl,_,Dr=svd(left.conj().T@y[1,:,:3]@right)
    cW=abs(np.vdot(left.conj().T@svd(y[0])[0][:,0],Dl[:,0]))
    fW=math.sqrt((1+cW)/2)
    F=((Dr*[cW,1.,fW])@Dr.conj().T).T
    assert np.linalg.matrix_rank(y[3][:,[0,2]])==2
    return dict(model=model,pt=pt,nu=nu,charged_summary=row,segments=segments,yW=y,S=S,
                F=F,cW=cW,fW=fW,N=y[3][:,[0,2]],Mweak=Mw,Mcolor=Mc,high=high,
                right_quark_to_GJ=right,left_quark_to_GJ=left,tagged_light_columns=Dr.conj())


def daughter(b, F_override=None, carry_weak_gain=False):
    """Apply a (possibly signed) five-bar embedding at Mw and the colored gate at Mc."""
    model,pt,nu=b["model"],b["pt"],b["nu"]
    flow=pt["flow"]
    F=b["F"] if F_override is None else np.asarray(F_override,dtype=complex)
    y=b["yW"].copy()
    y[2]=F.conj().T@y[2]
    y[3]=F.conj().T@y[3]
    poles=flow["fields"]["neutral_poles_GeV"]
    carrier=np.linalg.inv(b["S"]@pt["low_E"])
    Cparent=carrier.T@pt["C"]@carrier
    bparent=carrier.T@pt["b"]
    rho=nu["FROZEN"][0]
    mass=Cparent+rho*np.outer(bparent,bparent)
    mass_after=F.T@mass@F
    spec0=np.linalg.svd(mass,compute_uv=False)
    spec1=np.linalg.svd(mass_after,compute_uv=False)
    Q=np.linalg.qr(b["N"])[0][:,:2]
    proj=Q@Q.conj().T
    plane_error=np.linalg.norm(mass-proj.conj()@mass@proj,2)/np.linalg.norm(mass,2)
    gains=np.sqrt(np.maximum(0,np.linalg.eigvalsh(Q.conj().T@F@F.conj().T@Q)))
    zero=np.linalg.svd(mass)[2].conj().T[:,-1]
    Tdown=I.copy()
    res={}
    for start,susy,dt,g in reversed(b["segments"]):
        assert not susy
        y,Tdown=model["advance"](y,Tdown,-dt,g[::-1],flow)
        if start==b["Mcolor"]:
            left,val,right=svd(y[1])
            probability=abs(svd(y[0])[0][:,0].conj()@left)**2
            cc=np.sqrt(probability[0])
            fc=math.sqrt((1+cc)/2)
            factors=np.array([b["cW"],1.,b["fW"]] if carry_weak_gain else [cc,1.,fc])
            signed=probability*np.array([-1.,1.,-1.])
            parent_return=(signed@val)/(probability[0]*val[0])
            below_return=(signed@(val*factors))/(probability[0]*val[0])
            y[1,:,:3]=y[1,:,:3]@((right*factors)@right.conj().T)
            np.testing.assert_allclose(svd(y[1])[1],val*factors,rtol=2e-11)
            active_above_color=list(np.nonzero(np.linalg.norm(y[3],axis=0))[0])
        for j,pole in enumerate(poles):
            if start==pole:
                predicted=F.conj().T@flow["physics"]["Y0"][:,j]
                if np.linalg.norm(predicted):
                    res[str(j)]=np.linalg.norm(y[3,:,j]-predicted)/np.linalg.norm(predicted)
                y[3,:,j]=0.
        if start==b["Mcolor"]:
            active_below_color=list(np.nonzero(np.linalg.norm(y[3],axis=0))[0])
    assert active_above_color==[0,2] and active_below_color==[0]
    q,lep,gamma,L,ups=flow["extraction"](y)
    phases=np.diag(L.conj().T@pt["low_E"])
    L=L@np.diag(phases/abs(phases))
    U=np.linalg.solve(Tdown,L)
    H=carrier@F@U
    C=H.T@pt["C"]@H
    arm=H.T@pt["b"]
    np.testing.assert_allclose(C,U.T@F.T@Cparent@F@U,atol=3e-13)
    obs=nu["observables"](C,arm,rho)
    sol=nu["SOL"]
    cal=nu["observables"](pt["C"],pt["b"],rho)["scale_eV"]
    archive=read("../1739/results.json")["rows"][CFG["parent_point"]]
    for tag,sc in [("mother_solar_gain",cal),("archived_1739_gain",archive["neutrinos"]["scale_eV"])]:
        fac=(sc/obs["scale_eV"])**2
        obs["dm21_with_"+tag]=fac*sol
        obs["dm31_with_"+tag]=fac*obs["dm31_eV2"]
    zero_low=np.linalg.svd(C+rho*np.outer(arm,arm))[2].conj().T[:,-1]
    if np.linalg.cond(F)<1e8:
        pred=np.linalg.solve(U,np.linalg.solve(F,zero))
        assert abs(np.vdot(zero_low,pred/np.linalg.norm(pred)))>1-2e-10
    moduli=dict(zip(("Vus","Vub","Vcb"),np.sin(q[3:])*[np.cos(q[4]),1,np.cos(q[4])]))
    widths={k[1:]:(v-model["DIRECT"][k][0])/math.hypot(*model["DIRECT"][k][1:]) for k,v in moduli.items()}
    mk=model["MATCH"]; ks,_=mk["thresholds"](dict(moduli_MZ=moduli,weak_offsets=widths,gamma=gamma),mk["source_thresholds"]())
    shape=np.array([q[2]/q[1]/ks[2],1.,q[0]/q[1]/ks[0],ups[0]/q[1]/(ks[0]*ks[1])])
    profile=model["project_shared"](shape)
    absolute=model["absolute_table"](shape,ups[0]/ups[1],profile)
    old=model["JOINT"]["LEPTONS"]
    rel=lep/old[:,0]
    cp=model["CP"]["test"](moduli,gamma)
    central=model["JOINT"]["REFERENCE"]["FLOW"]["low"]["yukawa"]["up"]
    loop=model["CP"]["loop"]
    correction=loop(163.48*ups[2]/central[2][0],1.27*ups[1]/central[1][0],80.379)/loop(163.48,1.27,80.379)
    epsilon=np.dot(cp["tt_ut_parts"],correction)
    pull=(epsilon-model["CP"]["EPS_EXP"])/math.hypot(cp["correlated_error_at_fixed_CKM"],model["CP"]["EPS_ERROR"])
    return dict(weak_scale_c=b["cW"],color_scale_c=cc,c_difference_local_minus_shared=cc-b["cW"],
                named_quark_weak_halfangle=b["fW"],color_halfangle=fc,
                signed_return_before_and_after_actual_color_gate=[parent_return,below_return],
                lepton_filter_eigenvalues=np.linalg.eigvalsh((F+F.conj().T)/2),
                MZ_leptons=lep,MZ_relative_e_mu_after_tau_calibration=rel[:2]/rel[2]-1,
                MZ_absolute_marginal_shifts=(lep-old[:,0])/old[:,1], threshold_source_plane_gains=gains,
                neutral_mass_ratio_above_weak=spec0[0]/spec0[1],
                neutral_mass_ratio_below_weak=spec1[0]/spec1[1],
                effective_mass_outside_two_Dirac_plane=plane_error, matching_residual_neutral_columns=res,
                active_columns_on_color_sides=[active_above_color,active_below_color],
                MZ_down=q[:3],MZ_up=ups,MZ_gamma=gamma, MZ_moduli=moduli,epsilon_K=epsilon,epsilon_K_pull=pull,
                S_Bd=cp["sin2beta"],absolute_four_mass_ratios=absolute["absolute_at_Kuc_one"],
                absolute_four_ratio_radius=absolute["radius_at_Kuc_one"],profiled_four_ratio_radius=profile["radius"],
                MZ_zero_column_probabilities=abs(zero_low)**2,neutrinos=obs)


def spectrum_bounds(b,control,common):
    model,nu,flow=b["model"],b["nu"],b["pt"]["flow"]
    poles=flow["fields"]["neutral_poles_GeV"]
    total=math.log(b["Mweak"]/flow["MZ"])/(16*math.pi**2)
    first=math.log(b["Mweak"]/b["Mcolor"])/(16*math.pi**2)
    second=math.log(b["Mcolor"]/poles[0])/(16*math.pi**2)
    gmax=max(np.max(g[:,:2]) for _,_,_,g in b["segments"])
    assert gmax<.6
    ce=b["cW"]
    lepton_upper=np.linalg.norm(b["yW"][2],2)/ce
    both_upper=np.linalg.norm(b["N"],2)/ce
    one_upper=np.linalg.norm(b["N"][:,0])/ce
    # A fixed column obeys a left-linear beta. Drop nonnegative traces and
    # its positive neutral Gram term when bounding its norm on the downward leg.
    assert lepton_upper*math.exp((4.5*.6+1.5*.3**2)*total)<.02
    assert both_upper*math.exp((2.7*.6+1.5*.02**2)*first)<.30
    assert one_upper*math.exp((2.7*.6+1.5*.02**2)*(first+second))<.25
    L=1.5*.02**2*total+.5*(.30**2*first+.25**2*second)
    gain_ratio=common["threshold_source_plane_gains"][1]/common["threshold_source_plane_gains"][0]
    r0=control["neutral_mass_ratio_above_weak"]
    r1=common["neutral_mass_ratio_below_weak"]
    direct_reference=control["neutrinos"]["dm31_eV2"]
    solar,sol_err,_=nu["FREQ"]["juno_NO_dm21"]
    atm,atm_err,_=nu["FREQ"]["juno_NO_dm31"]
    accelerator,acc_err,_=nu["FREQ"]["accelerator_NO_dm32"]
    out={}
    for k in (1,2):
        sol_low,sol_high=solar-k*sol_err,solar+k*sol_err
        at_min,ac_min=atm-k*atm_err,accelerator-k*acc_err
        # Each physical Takagi singular value changes under congruence by at
        # most the square of a singular value of the transport. Use the actual
        # parent 8.374 TeV mass ratio for the sharpest single-leg estimate.
        plane_only=sol_high*r0**2*math.exp(4*L)*(1+2e-10)
        tagged=sol_high*r1**2*math.exp(4*L)*(1+2e-10)
        req=(at_min/plane_only)**.25
        rlow=max(at_min/sol_high,1+ac_min/sol_high)
        actual=common["neutrinos"]["dm31_eV2"]/solar
        gap=math.sqrt(rlow/actual)-1
        out[str(k)]=dict(solar_eV2=[sol_low,sol_high],juno_atmospheric_lower_eV2=at_min,
                    accelerator_dm32_lower_eV2=ac_min,minimum_mass_squared_ratio=rlow,
                    max_dm31_for_isometric_plane_eV2=plane_only,
                    any_shape_threshold_needed_from_dm31=req,
                    max_dm31_from_named_plane_singular_values_eV2=plane_only*gain_ratio**4,
                    max_dm31_from_actual_named_source_at_weak_eV2=tagged,
                    max_dm31_from_numerical_tagged_low_shape_eV2=sol_high*actual,
                    additional_relative_mass_ratio_increase_needed_for_two_clocks=gap,
                    split_only_dm31_at_solar_center_eV2=direct_reference)
    return dict(source_neutral_spectral_norm=np.linalg.norm(b["N"],2),
                source_lower_column_norm=np.linalg.norm(b["N"][:,0]),max_squared_weak_coupling=gmax,
                weak_to_color_and_color_to_lower_loop_time=[first,second],
                log_condition_upper_lepton_transport=L,named_source_singular_gain_ratio= gain_ratio,
                named_actual_Takagi_ratio_gain=r1/r0,intervals=out)


def solve(step):
    b=mother(step)
    control=daughter(b,I)
    common=daughter(b,b["F"])
    carried=daughter(b,b["F"],carry_weak_gain=True)
    prior=read("../1743/results.json")["cases"]
    compare={}
    for label,row in (("split",control),("common",common)):
        p=prior["complete_multiplets_"+("split_couplings" if label=="split" else "common_5bar")]
        compare[label]=dict(dm31_change_after_separating_GeV_poles=row["neutrinos"]["dm31_eV2"]-p["neutrinos"]["dm31_eV2"],
                            gamma_MZ_change_deg=row["MZ_gamma"]-p["MZ_gamma"])
    bounds=spectrum_bounds(b,control,common)
    report=bounds["intervals"]["2"]
    assert report["max_dm31_from_actual_named_source_at_weak_eV2"]<report["juno_atmospheric_lower_eV2"]
    assert report["max_dm31_from_named_plane_singular_values_eV2"]>report["juno_atmospheric_lower_eV2"]
    return dict(weak_and_color_thresholds_GeV=[b["Mweak"],b["Mcolor"]],soft_thresholds_GeV=[
                math.hypot(b["Mweak"],b["pt"]["flow"]["MS"]),math.hypot(b["Mcolor"],b["pt"]["flow"]["MS"])],
                mother_tan_beta=b["charged_summary"]["positive_soft_tan_beta"],
                cases={"same_spectrum_separate_couplings":control,"same_spectrum_shared_5bar":common},
                carrying_identical_weak_gain_to_color=dict(
                    change_four_mass_radius=carried["absolute_four_ratio_radius"]-common["absolute_four_ratio_radius"],
                    change_gamma_deg=carried["MZ_gamma"]-common["MZ_gamma"],
                    change_dm31_eV2=carried["neutrinos"]["dm31_eV2"]-common["neutrinos"]["dm31_eV2"]),
                compared_with_1743=compare,spectrum_necessary_bounds=bounds)


def tidy(x):
    if isinstance(x,dict): return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)): return [tidy(v) for v in x]
    if isinstance(x,np.integer): return int(x)
    if isinstance(x,(np.floating,float)): return float(f"{x:.12g}")
    if isinstance(x,np.bool_): return bool(x)
    return x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    coarse,fine=(solve(step) for step in CFG["steps"])
    case="same_spectrum_shared_5bar"
    a,b=(r["cases"][case] for r in (coarse,fine))
    err={"mass_radius":abs(a["absolute_four_ratio_radius"]-b["absolute_four_ratio_radius"]),
         "dm31_eV2":abs(a["neutrinos"]["dm31_eV2"]-b["neutrinos"]["dm31_eV2"]),
         "e_over_tau":abs(a["MZ_relative_e_mu_after_tau_calibration"][0]-b["MZ_relative_e_mu_after_tau_calibration"][0])}
    assert err["mass_radius"]<2e-5 and err["dm31_eV2"]<3e-9 and err["e_over_tau"]<2e-8
    result=tidy(dict(**fine,half_step_check=err))
    if args.check:
        assert result==read("results.json")
        print("1744 separate lepton/color thresholds, neutral shape and common quark data verified")
    else:
        (HERE/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf8")
