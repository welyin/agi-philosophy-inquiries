"""Restore the physical right columns across the split down-exchange thresholds."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())
LOOP = 16*math.pi**2
COLOR = np.array([4/15, 0, 2/3])


def read(key):
    return json.loads((HERE / CFG[key]).read_bytes())


def module(key):
    return runpy.run_path(str(HERE / CFG[key]))


def gram(y):
    return y @ y.conj().T


def project_down(yd, g, ratio):
    exchange = module("tree_exchange")
    left, sv, right = np.linalg.svd(yd)
    link = (right.conj().T * [-1., 1., -1.]) @ left.conj().T
    k = yd @ link / sv[0]
    a, w = exchange["dilation"](yd @ link, sv[0])
    heavy = np.column_stack((-math.sqrt(ratio**2-1)*w @ link.conj().T, a))
    _, masses, vh = np.linalg.svd(heavy, full_matrices=True)
    np.testing.assert_allclose(masses, [ratio]*3+[1.]*2, rtol=3e-13)
    high_basis, mid_basis, zero_basis = (vh.conj().T[:,sl] for sl in
                                          (slice(0,3), slice(3,5), slice(5,8)))
    full = np.column_stack((np.zeros((3,3)), g*w.conj().T))
    high, mid, light = (full @ v for v in (high_basis, mid_basis, zero_basis))
    k2 = k @ k
    for block, wanted in [(gram(high), g*g*k2/ratio**2),
                         (gram(mid), g*g*(np.eye(3)-k2)),
                         (gram(light), g*g*(1-ratio**-2)*k2),
                         (gram(full), g*g*np.eye(3))]:
        np.testing.assert_allclose(block, wanted, atol=3e-16)
    # The arbitrary right frames in the null space do not set the three physical
    # down labels. The caller retained the measured light columns and their frame.
    return k, mid


def low_flow(flow, rep, step, poles, neutral_y, mass, ratio, g, down_input,
             add_color=True, restore_columns=True):
    mz, soft = flow["MZ"], flow["MS"]
    eig = np.linalg.svd
    thresholds = list(flow["physics"]["thresholds"])
    if add_color:
        thresholds.extend([(mass,COLOR)]*2)
    y = np.zeros((4,3,5), complex)
    y[0,:,:3] = np.diag(flow["central"]("up"))
    y[1,:,:3] = flow["low_ckm"](rep["gamma_MZ_deg"],flow["angles"]("direct")) @ np.diag(down_input)
    E0 = np.array(rep["charged_left_MZ"]["real"])+1j*np.array(rep["charged_left_MZ"]["imag"])
    y[2,:,:3] = E0 @ np.diag(flow["central"]("charged_leptons"))
    lam = flow["physics"]["quartic"]["lambda_paper"][0]/4
    gauge = np.array([flow["physics"]["weak"][v][0]**2 for v in ("g1_GUT","g2","g3")])
    events = sorted(set([mz,mass,soft]+list(poles)+[m for m,_ in thresholds if mz<m<soft]))
    for start, end in zip(events, events[1:]):
        for i, p in enumerate(poles):
            if start==p:
                y[3,:,i] = neutral_y[:,i]
        if start==mass:
            yd = y[1,:,:3]
            d,vals,_ = eig(yd)
            up = eig(y[0,:,:3])[0][:,-1]
            weights = abs(up.conj()@d)**2
            match = {"mass_desc":vals, "prob_desc":weights, "return_over_yd":
                     np.dot(weights, vals*[-1,1,-1])/(weights[2]*vals[2])}
            if restore_columns:
                _, mid = project_down(yd,g,ratio)
                y[1,:,3:] = mid
        beta_g = np.array([4.1,-19/6,-7.])
        for m,c in thresholds:
            if m<=start:
                beta_g+=c
        n = max(1, math.ceil(math.log(end/start)/step))
        dt = math.log(end/start)/n/LOOP
        def rhs(state, la, gg):
            traces = np.array([gram(f) for f in state])
            t = np.array([3,3,1,1])@np.trace(traces,axis1=1,axis2=2).real
            t4 = np.array([3,3,1,1])@np.trace(traces@traces,axis1=1,axis2=2).real
            gy, gw = .6*gg[0],gg[1]
            dl = 24*la**2+(4*t-9*gw-3*gy)*la-2*t4+3/8*(2*gw**2+(gy+gw)**2)
            return flow["beta"](state,gg,False),dl
        for _ in range(n):
            midpoint, last = (1/(1/gauge-c*dt*beta_g) for c in (1,2))
            a, av = rhs(y,lam,gauge)
            b, bv = rhs(y+dt*a/2,lam+dt*av/2,midpoint)
            c, cv = rhs(y+dt*b/2,lam+dt*bv/2,midpoint)
            d, dv = rhs(y+dt*c,lam+dt*cv,last)
            y+=dt*(a+2*b+2*c+d)/6
            lam+=dt*(av+2*bv+2*cv+dv)/6
            gauge=last
    return match, {"lambda":lam,"yt":eig(y[0],compute_uv=False)[0],
                   "gauge":np.sqrt(gauge),"down_Gram_eig":np.linalg.eigvalsh(gram(y[1])),
                   "Higgs_down_trace":3*np.linalg.norm(y[1],"fro")**2}


def run():
    point = read("gauge_point")
    choice = point["central_conditional_combined_point"]
    poles = np.array(read("original_spectra")["shared_witness"]["neutral_poles_GeV"])
    flow, rep = module("charged_frame")["charged_background"](.10)
    down = flow["central"]("down")
    g, ratio = choice["g"],choice["R"]
    mass = choice["five_color_poles_GeV"][0]
    goal_b = g*math.sqrt(1-ratio**-2)
    target = point["central_combined_Higgs_GeV2"]
    neutral_y = flow["physics"]["Y0"]
    weights = np.sum(abs(neutral_y)**2,axis=0)*poles**2/(8*np.pi**2)
    q = read("higgs_curvature")["q0_GeV"]
    def curvature(s):
        return -s**3*np.dot(weights, np.log((s*poles/q)**2)-1)
    old_s = target["neutral_rescaling"]
    def scale_for(neutral_target):
        lo,hi=old_s-.01,old_s+.01
        assert curvature(lo)>neutral_target>curvature(hi)
        for _ in range(55):
            mid=(lo+hi)/2
            if curvature(mid)>neutral_target:
                lo=mid
            else:
                hi=mid
        return (lo+hi)/2
    s=scale_for(target["neutral"])
    def fitted(step):
        values = down.copy()
        def scan(vals, color, columns):
            return low_flow(flow,rep,step,s*poles,neutral_y*math.sqrt(s),mass,ratio,g,vals,color,columns)
        first = scan(values,True,False)[0]
        for _ in range(8):
            data,_=scan(values,True,False)
            vals,p=data["mass_desc"],data["prob_desc"]
            n_d=(p[1]*vals[1]-p[0]*vals[0])/p[2]
            if abs(vals[0]/goal_b-1)<3e-13 and abs(data["return_over_yd"])<3e-13:
                break
            values[[0,2]] *= [n_d/vals[2],goal_b/vals[0]]
        else:
            raise RuntimeError("rank saturation / signed return not closed")
        _,soft=scan(values,True,True)
        _,gauge_only=scan(values,True,False)
        old_soft=low_flow(flow,rep,step,poles,neutral_y,mass,ratio,g,down,False,False)[1]
        return first, values, data, soft, gauge_only, old_soft
    for _ in range(5):
        first, fitted_low, data, soft, gauge_only, old_soft = fitted(CFG["step_dlogmu"])
        k2=(data["mass_desc"]/goal_b)**2
        color=-3*mass**2*g**2/(8*np.pi**2)*(1.5+2*sum(k2)*math.log(ratio))
        residual=curvature(s)+color-target["target_negative_curvature"]
        if abs(residual)<1e-9:
            break
        s=scale_for(target["target_negative_curvature"]-color)
    else:
        raise RuntimeError("charged and neutral curvature shooting did not close")
    *_, fine, __, ___ = fitted(CFG["step_dlogmu"]/2)
    assert abs(soft["lambda"]-fine["lambda"])<1e-9
    assert abs(soft["yt"]-fine["yt"])<1e-9
    np.testing.assert_allclose(old_soft["yt"],rep["top_soft"],rtol=2e-8)
    np.testing.assert_allclose(old_soft["lambda"],rep["lambda_30TeV_full_flavour"],rtol=2e-8)
    mod=read("charged_measurements")["yukawa"]["down"]
    widths=(fitted_low-np.array(mod)[:,0])/np.array(mod)[:,1]
    assert max(abs(widths))<1
    end = flow["END"]
    old_g2=flow["segments"](.1)[-1][-1][-1]
    masslist=np.array(choice["five_color_poles_GeV"])
    loss=COLOR*np.log(end/masslist).sum()/(2*np.pi)
    g_end=np.sqrt(4*np.pi/(4*np.pi/old_g2-loss))
    old_high = np.array(rep["Y_down_at_10PeV"])
    soft_condition=math.sqrt(soft["down_Gram_eig"][-1]/soft["down_Gram_eig"][0])
    required_spread=math.log((old_high[-1]/old_high[0])/soft_condition)
    # Gauge-only diagnostic after the soft scale; it is separate from Yukawa running.
    return {"original_and_updated_neutral_scale": [old_s,s],
            "1712_actual_cluster_Yukawa_squared":[np.sum(abs(neutral_y[:,:2])**2),np.sum(abs(neutral_y[:,2:])**2)],
            "old_s_neutral_curvature_shift_GeV2":curvature(old_s)-target["neutral"],
            "refitted_color_curvature_GeV2":color,
            "refitted_neutral_curvature_GeV2":curvature(s),
            "new_combined_curvature_residual_GeV2":residual,
            "frozen_light_retune": {"before_matching_b_over_fixed_target_minus_1":first["mass_desc"][0]/goal_b-1,
                "before_return_over_yd":first["return_over_yd"], "down_yukawas_MZ_after":fitted_low,
                "down_yukawas_MZ_before":down,"after_widths_each":widths,"return_over_yd":data["return_over_yd"],
                "after_matching_b_fractional_residual":data["mass_desc"][0]/goal_b-1},
            "common_tree_projection_desc_b_s_d":{
                "Y_light":data["mass_desc"],"Y_mid":g*np.sqrt(1-k2/ratio**2),"Y_above":np.full(3,g),
                "intermediate_to_light_trace_ratio":sum(1-k2/ratio**2)/sum(k2*(1-ratio**-2)),
                "two_low_pole_fraction_of_full_trace":sum(1-k2)/3,
                "three_high_pole_fraction_of_full_trace":sum(k2)/(3*ratio**2),
                "commutator_determinant_mid_over_light":-1/(ratio**2-1)**3},
            "original_1712_10PeV_effective_down_singular_asc":old_high,
            "common_factor_min_worst_fractional_error_for_tree_iso_vs_old_singular":(old_high[-1]-old_high[0])/(old_high[-1]+old_high[0]),
            "unprojected_high_texture_necessary_log_distortion":{
                "soft_active_condition_number":soft_condition,"old_10PeV_condition_number":old_high[-1]/old_high[0],
                "integrated_generator_spreads_plus_threshold_terms_min":required_spread,
                "mean_dimensionless_loop_generator_spread_if_zero_matching_cost":required_spread*LOOP/math.log(end/flow["MS"])},
            "gauge_10PeV_old_new":{"old":np.sqrt(old_g2),"actual":g_end,"inverse_alpha_loss":loss},
            "one_light_Higgs_below_30TeV":{"old_1712":old_soft,"new_gauge_neutral_but_only_light_columns":gauge_only,
                "actual_active_columns":soft, "half_step_abs_delta_lambda":abs(soft["lambda"]-fine["lambda"]),
                "half_step_abs_delta_top":abs(soft["yt"]-fine["yt"])}}


def tidy(x):
    if isinstance(x,dict): return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)): return [tidy(a) for a in x]
    if isinstance(x,(float,np.floating)): return float(f"{x:.11g}")
    return x


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    check=ap.parse_args().check
    results=tidy(run())
    out=HERE/"results.json"
    if check:
        assert results==json.loads(out.read_bytes())
        print("1721 matched the rank-five split spectrum, original neutral arms, and charged/soft transport")
    else:
        with out.open("w",encoding="utf8") as target:
            json.dump(results,target,indent=2,allow_nan=False)
            target.write("\n")
