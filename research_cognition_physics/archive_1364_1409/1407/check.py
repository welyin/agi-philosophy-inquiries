"""Mixed/diagonal contact matching for one real heavy vector, with working boxes."""
import argparse
import cmath
import json
from math import asin, cos, hypot, isclose, log, pi, sqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


cfg = read("inputs.json")
frozen = read(cfg["inputs1404"])
measure = read(cfg["inputs1406"])
tree = read(cfg["tree_inputs1381"])
const = frozen["const_GeV"]
gf, mw = const["GF"], const["MW"]
v = (sqrt(2)*gf)**(-.5)
ca, ca_err = frozen["short_distance"]["CA_at5GeV"]
loop, loop_err = cfg["Rloop_SM"]
w, base = cfg["trident_axial_weight"], 1+4*cfg["weak_sine2"]


def norm(ca_value):
    return 2*sqrt(2)*gf*mw**2*ca_value/pi


def tri(u):
    return (w*(1+u)**2+(base+u)**2)/(w+base**2)


def umax(cap):
    assert cap >= 1
    linear, quad = 2*(w+base), 1+w
    c = (cap-1)*(w+base**2)
    return 2*c/(sqrt(linear**2+4*quad*c)+linear)


def rg(alpha0, mass):
    def move(a, a_scale, b_scale, b0):
        return 1/(1/a+b0*log(b_scale/a_scale)/(2*pi))
    scales = cfg["running"]
    mb, mz, mt = scales["mb"], scales["MZ"], scales["mt"]
    ab = move(alpha0,mz,mb,23/3)
    at = move(alpha0,mz,mt,23/3)
    high = move(at,mt,mass,7.)
    return (high/at)**(2/7)*(at/ab)**(6/23)


def omega(running, weak, ref_loop):
    return running * norm(weak)**2 / (4*pi**2*ref_loop)


def working_box(k):
    meas, err = frozen["osc_per_ps"]["s"]
    sm, ckm_lo, ckm_hi, had_lo, had_hi = cfg["DM_SM_ps_inv"]
    rlo = (meas-k*err)/(sm+k*(ckm_hi+had_hi))
    rhi = (meas+k*err)/(sm-k*(ckm_lo+had_lo))
    # Bound the SM phase by the unitarity triangle without using the mixing observable.
    form, width_form = tree["K_form"]
    Kmu, width_k = tree["K_product"]
    vus = (Kmu+k*width_k)/(form-k*width_form)
    vts, _, vts_up = cfg["HPQCD_inputs"]["tree_Vts"]
    vcs = sqrt(1 - vus**2 - (vts+k*vts_up)**2)
    ub_cb, *errs = tree["ub_over_cb"]
    beta_cap = 2*asin((ub_cb+k*sum(errs))*vus/vcs)
    phase, stat, syst = cfg["phi_measurement_rad"]
    phase_bound = abs(phase)+k*(stat+syst)+beta_cap+cfg["extra_control_phase_cap_rad"]
    hmax = max(hypot(r-1,sqrt(2*r*(1-cos(phase_bound)))) for r in (rlo,rhi))
    tri_c, tri_err = cfg["trident_CCFR"]
    umax_here = umax(tri_c+k*tri_err)
    assert isclose(tri(umax_here),tri_c+k*tri_err,rel_tol=1e-12)
    om_min = omega(cfg["LL_running_budget"][0], ca-k*ca_err, loop+k*loop_err)
    pmax = sqrt(hmax*umax_here/om_min)
    real_pmax = sqrt((rhi-1)*umax_here/om_min)
    # Independent true SM branching, without treating the measured mixing frequency as its SM value.
    short = frozen["short_distance"]
    S, deltaS = short["S0_published"]
    eta_b, delta_eta = short["etaB"]
    bag, dbag = frozen["bag"]["s"]
    tau, dtau = frozen["tau_ps"]["s_heavy"]
    Kmix = read("../1404/results.json")["K_SM"]
    upper_K = Kmix*((ca+k*ca_err)/ca)**2*S/(S-k*deltaS)*eta_b/(eta_b-k*delta_eta)
    dimuon_upper = upper_K*(sm+k*(ckm_hi+had_hi))*(tau+k*dtau)/(bag-k*dbag)
    dimuon_upper *= sqrt(1-(2*const["m_mu"]/const["M_Bs"])**2)*1e9*(1+cfg["dim_mu_nonfactor_budget"])
    cms = measure["CMS"]
    obs_lower = (cms["branching_in_1e-9"]-k*sum(cms["lower_widths"]))/dimuon_upper
    # All allowed H/L weights in the reference acceptance produce at most the heavy yield.
    reflected_distance = sqrt(1 + obs_lower)
    real_reflected_distance = 1+sqrt(obs_lower)
    phase_free_pmax = sqrt((1+rhi)*umax_here/om_min)
    K = measure["sm_coefficient_ratio"]["range"][1]
    rnubound = (2+(1+K*pmax)**2)/3
    belle = measure["BelleII"]
    knu_upper = belle["reference_SD_in_1e-6"]*((vts+k*vts_up)/cfg["Knu_original_CKM"])**2
    knu_upper *= 1+k*cfg["Knu_other_relative_per_width"]
    belle_low = (belle["mu"]-k*belle["width"])*belle["reference_SD_in_1e-6"]
    return {"SM_dimuon_upper_in_1e9":dimuon_upper,"abs_P_squared_lower":obs_lower,
            "SM_Knu_upper_in_1e6":knu_upper,"NP_Knu_upper_in_1e6":rnubound*knu_upper,
            "Belle_lower_in_1e6":belle_low,"Belle_separated":rnubound*knu_upper<belle_low,
            "frequency_ratio":[rlo,rhi],"abs_SM_reference_phase_max":beta_cap,
            "abs_additional_mixing_phase_max":phase_bound,"abs_h_max":hmax,"trident_u_max":umax_here,
            "min_M_over_gl_GeV":v/sqrt(umax_here),"min_omega":om_min,
            "abs_P_minus_1_max":pmax,"phase_free_max":phase_free_pmax,"real_max":real_pmax,
            "min_distance_if_ReP_nonpositive":reflected_distance,
            "min_distance_real_reflected":real_reflected_distance,
            "max_h_for_excluding_left_half_plane":om_min*reflected_distance**2/umax_here,
            "max_h_for_excluding_real_reflection":om_min*real_reflected_distance**2/umax_here,
            "independent_Belle_Rnu_upper":rnubound,"half_plane_separated":pmax<reflected_distance,
            "real_branch_separated_without_phase":real_pmax<real_reflected_distance}


def witness():
    q = cfg["witness"]
    M,u,P,bag,fb = q["M_GeV"],q["u"],q["P"],q["Bhat_s"],q["fBs_GeV"]
    eta = rg(q["alpha_s_MZ"],M)
    h = omega(eta,ca,loop)*(P-1)**2/u
    measured = frozen["osc_per_ps"]["s"][0]
    sm_mix = measured/(1+h)
    fhat = fb*sqrt(bag)
    ref_f, ef = cfg["HPQCD_inputs"]["fBs_sqrtB_GeV"]
    assert abs(fhat-ref_f) <= 2*ef
    short = frozen["short_distance"]
    coeff = gf**2*mw**2*const["M_Bs"]*short["S0_published"][0]*short["etaB"][0]/(6*pi**2)
    lam = sqrt(sm_mix*cfg["hbar_GeV_ps"]/(coeff*fhat**2))
    def ckm(s23):
        a, phase, ratio = q["s12"], cmath.exp(1j*q["quark_phase"]), q["ub_over_cb"]
        c12, c23 = sqrt(1-a*a), sqrt(1-s23*s23)
        s13 = ratio*s23/sqrt(1+(ratio*s23)**2)
        c13 = sqrt(1-s13*s13)
        return {"Vts":-c12*s23-a*c23*s13*phase,"Vtb":c23*c13,"Vcb":s23*c13,
                "Vtd":a*s23-c12*c23*s13*phase,"Vus":a*c13,"Vub":s13/phase,
                "Vcs":c12*c23-a*s23*s13*phase}
    lower,upper = .03,.05
    for _ in range(65):
        mid = (lower+upper)/2
        mix = ckm(mid)
        if abs(mix["Vts"]*mix["Vtb"]) < lam:
            lower=mid
        else:
            upper=mid
    mix = ckm((lower+upper)/2)
    vts = abs(mix["Vts"])
    assert isclose(abs(mix["Vts"]*mix["Vtb"]),lam,rel_tol=1e-12)
    assert isclose(abs(mix["Vub"]/mix["Vcb"]),q["ub_over_cb"])
    reference_phase = -2*cmath.phase(-mix["Vts"]*mix["Vtb"]/(mix["Vcs"]*mix["Vcb"]))
    ph,es,ey=cfg["phi_measurement_rad"]
    assert abs(reference_phase-ph)<2*(es+ey)
    for name in ("Vtb","Vtd"):
        center,down,up=cfg["HPQCD_inputs"]["tree_"+name]
        assert center-2*down<=abs(mix[name])<=center+2*up
    fk,efk = tree["K_form"]
    observed_k,ek = tree["K_product"]
    assert (observed_k-2*ek)/(fk+2*efk)<=abs(mix["Vus"])<=(observed_k+2*ek)/(fk-2*efk)
    ref,lo,hi = cfg["HPQCD_inputs"]["tree_Vts"]
    assert ref-2*lo <= vts <= ref+2*hi
    bs = read("../1404/results.json")["K_SM"]*sm_mix*frozen["tau_ps"]["s_heavy"][0]/bag
    bs *= sqrt(1-(2*const["m_mu"]/const["M_Bs"])**2)*1e9
    b_ref = measure["sm_Bs"]
    assert abs(bs-b_ref["branching_in_1e-9"])<=2*b_ref["width"]
    assert abs(bag-frozen["bag"]["s"][0])<2*frozen["bag"]["s"][1]
    obs = bs*P**2
    for d in (measure["CMS"],measure["LHCb"]):
        center = d["branching_in_1e-9"]
        assert center-2*sum(d["lower_widths"]) <= obs <= center+2*sum(d["upper_widths"])
    lep = M*sqrt(u)/v
    quark = -norm(ca)*(P-1)*sqrt(2)*gf*M**2*lam/(pi*lep)
    D = eta/(4*sqrt(2)*gf*M**2*loop)*(quark/lam)**2
    assert isclose(D,h)
    result_nu = (2+(1-measure["sm_coefficient_ratio"]["center"]*(P-1))**2)/3
    belle = measure["BelleII"]
    scale = belle["reference_SD_in_1e-6"]
    bare_knu = scale*(lam/cfg["Knu_original_CKM"])**2
    knu = bare_knu*result_nu
    assert scale*(belle["mu"]-2*belle["width"]) <= knu <= scale*(belle["mu"]+2*belle["width"])
    assert tri(u)<=cfg["trident_CCFR"][0]+2*cfg["trident_CCFR"][1]
    return {"u":u,"P":P,"eta":eta,"h":h,"trident_rate_ratio":tri(u),"true_SM_DM_per_ps":sm_mix,
            "fBs_sqrtB_GeV":fhat,"Vts_abs":vts,"true_SM_Bs_in_1e9":bs,"observed_Bs_in_1e9":obs,
            "abs_CKM_witness":{n:abs(x) for n,x in mix.items()},"phi_SM_rad":reference_phase,
            "abs_lambda_t":lam,"true_SM_Knu_in_1e6":bare_knu,"actual_Knu_in_1e6":knu,
            "Rnu":result_nu,"gl":lep,"abs_gq":abs(quark),"arg_1_plus_h":0.}


def calculate():
    alpha_lo, alpha_hi = cfg["strong_alpha_MZ_range"]
    mlo,mhi = cfg["running"]["M_in_GeV"]
    etas = [rg(alpha_hi,mhi),rg(alpha_lo,mlo)]
    assert cfg["LL_running_budget"][0] < etas[0] < etas[1] < cfg["LL_running_budget"][1]
    boxes = {str(k):working_box(k) for k in (1,2)}
    assert boxes["1"]["half_plane_separated"] and not boxes["2"]["half_plane_separated"]
    # Dimensionless couplings at equal masses. L is the mixed coefficient.
    quarks, leptons = (1,1j), (1,0)
    D = sum(x*x/2 for x in quarks)
    T = sum(y*y for y in leptons)
    L = sum(x*y for x,y in zip(quarks,leptons))
    assert D == 0 and L == T == 1
    return {"LL_QCD_extrema":etas,"boxes":boxes,"common_two_width_real_witness":witness(),
            "two_vector_example":{"D_real":D.real,"trident":T,"semileptonic":L.real}}


def agrees(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(agrees(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(agrees(x,y) for x,y in zip(a,b))
    return isclose(a,b,rel_tol=1e-10,abs_tol=1e-20)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--check",action="store_true")
    modes.add_argument("--write-new",action="store_true")
    args = parser.parse_args()
    result=calculate()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as dest:
            json.dump(result,dest,indent=2,allow_nan=False)
            dest.write("\n")
    else:
        assert agrees(result,read("results.json"))
        print("1407 single-vector phases, trident and common real-branch witness verified")
