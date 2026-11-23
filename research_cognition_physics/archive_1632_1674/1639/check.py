"""Cover the entire complex strange-bottom pole at each fixed spectrum with enclosing disks."""
import argparse
import cmath
import itertools
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
def read(path):
    return json.loads((HERE / path).read_bytes())


def bounds(row, n):
    return row[0] - n*sum(row[1:]), row[0] + n*sum(row[1:])


def phase_support(lo, hi, direction):
    """Enclose cos(phi-direction) and its endpoints for all angles in [lo,hi]."""
    vals = [math.cos(lo-direction), math.cos(hi-direction)]
    for k in range(-3, 4):
        if lo <= direction+k*math.pi <= hi:
            vals.append((-1)**k)
    return np.array([min(vals), max(vals)])


def domain(cfg, high, n, k_range):
    spec, pole, mix, had, clk, weak, cal, rare, ff = [
        read(cfg[k]) for k in ("weak_fixture", "weak_poles", "weak_operators", "lattice",
                               "mixing_clock", "absolute_weak", "rare_calibration",
                               "absolute_rates", "rare_form_factors")]
    w = dict(spec["witness"])
    if high:
        w.update(spec["high_mass"])
    old = runpy.run_path(str((HERE/"../1635/check.py").resolve()))
    u = old["ckm"](*w["ordinary_rotation_s12_s23_s13_delta"])
    ms, mb = w["light_masses_GeV"][1:]
    M = w["M_D_GeV"]
    rb = math.sqrt(pole["right_doublet_probability"])
    ell_b = mb*rb/M
    sb2 = 2*(pole["left_difference"]/pole["k"]-ell_b**2)
    h = math.sqrt(1-sb2-ell_b**2)
    lt, lc = [u[j, 1].conjugate()*u[j, 2]*h for j in (2, 1)]
    # For all physical r_s, |ell_s| <= ms/M; polar sqrt is Lipschitz on this cone.
    ls = ms/M
    min_sq = 1-sb2-ell_b**2-ls**2
    dh = (2*ls*ell_b + ls**2)/(math.sqrt(min_sq) + h)
    dl = 2*dh
    ax = lc.conjugate()**2/abs(lc)**2
    dax = 4*dl/(abs(lc)-dl)
    d = weak["const_GeV"]
    gf, mw, mbs = d["GF"], d["MW"], d["M_Bs"]
    norm = gf**2*mw**2/(4*math.pi**2)
    pre = 2*mbs/(3*had["hbar_GeV_ps"])
    eta, s0 = weak["short_distance"]["etaB"], weak["short_distance"]["S0_published"]
    bs = had["HPQCD_inputs"]["fBs_sqrtB_GeV"]
    t_lo, t_hi = np.prod(np.array([bounds(eta,n), bounds(s0,n)]), axis=0)
    zrange = clk["eta_Z_RGI_range"]
    a4,a5 = mix["LR_B4_B5_weights"]
    b4,b5 = mix["LR_F2B4_B5_GeV2"]
    urange = a4*np.array(bounds(b4,n))+a5*np.array(bounds(b5,n))
    refs = ff["HPQCD_tableI_eq15"]
    gnu = refs["Xt_ref"]/(math.pi*refs["sine_squared_ref"]*refs["alpha_inverse_ref"])
    L, dL = -gnu*lt, gnu*dl
    k = sum(k_range)/2
    dk = (k_range[1]-k_range[0])/2
    e = 2*ms*mb/M**2
    alpha = clk["alpha_s_ref"]
    mt, mz = mix["mt_GeV"], mix["MZ_GeV"]
    run_alpha = runpy.run_path(str((HERE/"../1636/check.py").resolve()))["alpha_step"]
    at = run_alpha(alpha,mz,mt,5)
    am = run_alpha(at,mt,M,6)
    qcd = (at/alpha)**(6/23)*(am/at)**(6/21)
    light = (1-ms**2/M**2)*(1-mb**2/M**2)
    xt = weak["short_distance"]["xt"][0]
    h1 = -(xt-7)/(4*(xt-1))-(xt**2-2*xt+4)/(2*(xt-1)**2)*math.log(xt)
    lr = -norm*xt*(math.log(M/mw)+h1)*light
    hard = mix["RR_hard_Qd_color_weak_factor"]*(2*(M/mix["vev_GeV"])**2*light)**2/((4*math.pi)**2*M**2)
    zz = gf/math.sqrt(2)*(1+e**2) + qcd*hard
    coefs = []
    coefficient_errs = []
    for f,t,z,bag in itertools.product(bounds(bs,n), (t_lo,t_hi), zrange, urange):
        c = pre*np.array([f**2*norm*t*lt**2, bag*lr*lt, f**2*z*zz])
        dlparts = pre*np.array([f**2*norm*t*(2*abs(lt)+dl)*dl, abs(bag*lr)*dl, 0])
        coefs.append(ax*c)
        coefficient_errs.append(abs(c)*dax + dlparts)
    coef = np.array(coefs)
    coef_err = np.max(coefficient_errs,axis=0)
    ome = bounds(clk["omega_ps_inverse"],n)
    phases = [bounds(clk["phase_inputs_rad"][j],n) for j in
              ("LHCb_KK_2015_18","CMS_KK_13TeV_final")]
    cap = clk["extra_decay_phase_cap_rad"]
    phi_range = [max(i[0] for i in phases)-cap, min(i[1] for i in phases)+cap]
    directions = np.linspace(0,np.pi,cfg["projection_directions"],endpoint=False)
    project = np.exp(-1j*directions)
    clock_range = []
    for alpha in directions:
        support = phase_support(*phi_range,alpha)
        values = np.outer(ome,support)
        clock_range.append([np.min(values),np.max(values)])
    clock_range = np.array(clock_range)
    cms_record, r3 = cal["CMS"], weak["LHCb_run3_2024_25"]["s_in_1e9"]
    cms = (cms_record["branching_in_1e-9"]-n*sum(cms_record["lower_widths"]),
           cms_record["branching_in_1e-9"]+n*sum(cms_record["upper_widths"]))
    lhcb = (r3[0]-n*(r3[1]+r3[3]), r3[0]+n*(r3[2]+r3[4]))
    fbs = ff["fBs_GeV"]
    tau = weak["tau_ps"]["s_heavy"]
    scales = 1e9*rare["Hmu_heavy"]*(np.array(bounds(fbs,n))/fbs[0])**2
    scales *= np.array(bounds(tau,n))/tau[0]
    scales *= 1 + np.array([-1,1])*ff["mu_rate_transfer_width"]
    muons = [(cal["CMS"]["accepted_light_to_heavy_weight_range"][0],cms),
             (ff["run3_working_light_response"],lhcb)]
    # The positive semidefinite spin-orthogonal scalar may only raise a rate.
    upper_p = min(math.sqrt(r[1]/(scales[0]*weight)) for weight,r in muons)
    fmax = (upper_p+k_range[1]*(abs(L)+dL))/(1-e) + 1e-12
    b_lo, b_hi = bounds([cal["BelleII"]["mu"],cal["BelleII"]["width"]],n)
    b_norm = rare["Hnu"]/(1e-6*cal["BelleII"]["reference_SD_in_1e-6"])
    b_norm_bounds = b_norm*(1+n*refs["hadronic_rate_width"]*np.array([-1,1]))
    phi_mid = sum(phi_range)/2
    axis_mid = ax*cmath.exp(-1j*phi_mid)
    axis_radius = dax+2*math.sin((phi_range[1]-phi_range[0])/4)
    pad = cfg["pad_per_ps"]

    def rule(F, half):
        rad = math.sqrt(2)*half
        absF = abs(F)
        if absF-rad > fmax+1e-12:
            return "outside"
        N = L+(1+e)*F
        rn = dL+(1+e)*rad
        if ((abs(N)+rn)**2*b_norm_bounds[1] < b_lo-1e-9 or
                max(0,abs(N)-rn)**2*b_norm_bounds[0] > b_hi+1e-9):
            return "Belle"
        P = k*L-(1-e)*F
        rp = k_range[1]*dL+dk*abs(L)+(1-e)*rad
        p2err = (2*abs(P)+rp)*rp
        for weight, window in muons:
            low = ((1+weight)*abs(P)**2 + (1-weight)*(P*P*axis_mid).real)/2
            err = p2err+(1-weight)*abs(P)**2*axis_radius/2
            scalar = ((mbs/w["h_mass_GeV"])**2*(absF+rad))**2
            if (scales[1]*((abs(P)+rp)**2+scalar) < window[0]-1e-8 or
                    scales[0]*(low-err) > window[1]+1e-8):
                return "dimuon"
        poly = coef[:,0]+coef[:,1]*F+coef[:,2]*F**2
        curve_rad = (np.max(np.abs(coef[:,1]))+np.max(np.abs(coef[:,2]))*(2*absF+rad))*rad
        error = (curve_rad + float(coef_err@np.array([1.,absF+rad,(absF+rad)**2]))
                 + pad + cfg.get("extra_mix_radius_ps_inverse",0.))
        proj = np.real(poly[:,None]*project)
        if np.any(np.max(proj,axis=0)+error < clock_range[:,0]-pad):
            return "mix"
        if np.any(np.min(proj,axis=0)-error > clock_range[:,1]+pad):
            return "mix"
        return None

    return rule, fmax, {"lambda_product_disk":dl,"control_axis_disk":dax,
        "mix_phase_range_rad":phi_range,"omega_range_ps_inverse":list(ome),
        "F_max":fmax, "center_lambda_t":[lt.real,lt.imag],"center_lambda_c":[lc.real,lc.imag],
        "frequency_roundoff_pad":pad,"k_range":list(k_range),"Hmu_1e9_bounds":scales.tolist(),
        "Belle_relative_hadronic_halfwidth":n*refs["hadronic_rate_width"]}


def cover(cfg, high, n, k_range):
    rule, edge, info = domain(cfg,high,n,k_range)
    by_reason = dict(outside=0,Belle=0,dimuon=0,mix=0,unresolved=0)
    remaining = []
    stack = [(0j,edge,0)]
    max_depth = 0
    while stack:
        F,half,depth = stack.pop()
        max_depth = max(max_depth,depth)
        passed = rule(F,half)
        if passed:
            by_reason[passed] += 1
        elif depth == cfg["max_depth"]:
            if len(remaining)<3:
                remaining.append([F.real,F.imag,half])
            by_reason["unresolved"] += 1
        else:
            h = half/2
            stack.extend((F+h*(sx+1j*sy),h,depth+1) for sx,sy in itertools.product((-1,1),repeat=2))
    return dict(**info,depth=max_depth,leaf_counts=by_reason,unresolved_examples=remaining)


def normalize(value):
    if isinstance(value, dict):
        return {k:normalize(v) for k,v in value.items()}
    if isinstance(value, list):
        return [normalize(v) for v in value]
    if isinstance(value, float):
        return float(f"{value:.11g}")
    return value


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    cfg=read("inputs.json")
    kr=read(cfg["rare_calibration"])["sm_coefficient_ratio"]["range"]
    out={}
    for name,high in [("common_witness",False),("high_mass_witness",True)]:
        out[name]={str(n):cover(cfg,high,n,kr) for n in (1,2)}
        for key in ("1","2"):
            assert out[name][key]["leaf_counts"]["unresolved"] == 0
        residual=cfg["residual_test"]
        delta=residual["two_width_extra_ps_inverse"][int(high)]
        extra_cfg=dict(cfg,extra_mix_radius_ps_inverse=delta,max_depth=residual["max_depth"])
        edges=np.linspace(*kr,residual["k_slices"]+1)
        counts=[]
        for a,b in zip(edges[:-1],edges[1:]):
            data=cover(extra_cfg,high,2,[a,b])
            assert data["leaf_counts"]["unresolved"] == 0
            counts.append(sum(data["leaf_counts"].values()))
        out[name]["two_width_added_complex_mix_disk"]={"excluded_radius_ps_inverse":delta,
            "k_edges":edges.tolist(),"enclosing_leaves_per_slice":counts}
    fixed=runpy.run_path(str(HERE/"frequency_certificate.py"))["calculate"]()
    for key in out:
        out[key]["one_width_frequency_lower_at_central_k_ps_inverse"]=fixed[key]["1"]["clock_lower_ps_inverse"]
    result=normalize(out)
    for key in out:
        print(key,[(n,out[key][n]["depth"],out[key][n]["leaf_counts"]) for n in ("1","2")])
    if args.check:
        assert result==read("results.json")
    else:
        (HERE/"results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")

if __name__=="__main__":
    main()
