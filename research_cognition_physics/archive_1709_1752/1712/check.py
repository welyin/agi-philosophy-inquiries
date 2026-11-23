"""Shoot y_d from reciprocal zero return on the same neutral-cluster boundary."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


CFG = read("inputs.json")
MODEL = runpy.run_path(str(HERE / CFG["generator"]))
FL = MODEL["fl"]
PORT = MODEL["portal"]
SPECTRUM = read(CFG["measured_masses"])["yukawa"]
DIRECT = read(CFG["charged_moduli"])["direct_semileptonic_moduli"]


def select_upper_yukawas(rep, step):
    y = np.zeros((4, 3, 4), complex)
    y[0, :, :3] = np.diag(FL["central"]("up"))
    y[1, :, :3] = FL["low_ckm"](rep["gamma_MZ_deg"], FL["angles"]("direct")) @ np.diag(FL["central"]("down"))
    y[2, :, :3] = MODEL["sim"]["zer"]["unpack"](rep["charged_left_MZ"]) @ np.diag(FL["central"]("charged_leptons"))
    poles = MODEL["spectra"]["neutral_poles_GeV"]
    for start, susy, dt, gauge in FL["segments"](step):
        if start >= poles[2]:
            up, _ = FL["basis"](y[0])
            dw, down = FL["basis"](y[1])
            weights = abs(up[:, 0].conj() @ dw)**2
            return weights, down
        assert not susy
        for j, pole in enumerate(poles):
            if start == pole:
                y[3, :, j] = FL["physics"]["Y0"][:, j]
        for j in range(0, len(gauge)-1, 2):
            beta = FL["beta"]
            k1 = beta(y, gauge[j], False)
            k2 = beta(y+dt*k1/2, gauge[j+1], False)
            k3 = beta(y+dt*k2/2, gauge[j+1], False)
            k4 = beta(y+dt*k3, gauge[j+2], False)
            y += dt*(k1+2*k2+2*k3+k4)/6
    raise ValueError("upper cluster absent")


def run_point(sign, us=0., strange=0., bottom=0., k=1, ub=-1., cb=0., step=.20):
    widths = {"Vus": us, "Vub": ub, "Vcb": cb, "ys": strange, "yb": bottom}
    port_values = FL["cfg"]["direct_semileptonic_moduli"]
    port_values["Vcb"][0] = DIRECT["Vcb"][0] + cb * math.hypot(*DIRECT["Vcb"][1:])
    for i, z in [(1, strange), (2, bottom)]:
        central, sigma = SPECTRUM["down"][i]
        FL["low"]["yukawa"]["down"][i][0] = central + z*sigma
    PORT.update({"yd_widths": 0., "Vub_widths": ub, "phase_k": k, "third_port": CFG["third_port"]})
    # Remove the R1702 mass-calibrated root as a numerical seed on this branch.
    graph_context = MODEL["state"]
    old_read = graph_context["read"]
    def uncalibrated_start(name):
        if name == graph_context["cfg"]["endpoints"]:
            return {"charged_endpoints": {"direct": {str(k): {"gamma_MZ_deg": 60., "large_modulus": 1.}}}}
        return old_read(name)
    graph_context["read"] = uncalibrated_start
    case_context = MODEL["case"].__globals__
    previous_upper = case_context["upper_cabibbo"]
    boundary = {}
    def read_same_scale(rep, grid):
        weights, masses = select_upper_yukawas(rep, grid)
        boundary.update(weights=weights, masses=masses, gamma=rep["gamma_MZ_deg"])
        high_sines = np.sin(rep["CKM_angles_at_10PeV_rad"])
        ub2 = high_sines[1]**2
        us2 = high_sines[0]**2 * (1-ub2)
        boundary["uv_weights"] = np.array([1-us2-ub2, us2, ub2])
        boundary["uv_masses"] = np.array(rep["Y_down_at_10PeV"])
        return math.sqrt(weights[1])
    case_context["upper_cabibbo"] = read_same_scale
    # Start with the mass formula evaluated at the measured light charged vertex.
    p_s = (DIRECT["Vus"][0] + us * math.hypot(*DIRECT["Vus"][1:]))**2
    p_b = (DIRECT["Vub"][0] + ub * math.hypot(*DIRECT["Vub"][1:]))**2
    s0, b0 = FL["central"]("down")[1:]
    d = (p_s*s0+sign*p_b*b0)/(1-p_s-p_b)
    try:
        for iteration in range(1, 9):
            FL["low"]["yukawa"]["down"][0] = [float(d), 0.]
            output = MODEL["case"](CFG["neutral_rule"], us, step=step)
            weights, masses = boundary["weights"], boundary["masses"]
            target = (weights[1]*masses[1]+sign*weights[2]*masses[2])/weights[0]
            residual = target / masses[0] - 1
            if abs(residual) < 3e-12:
                break
            d *= 1 + residual
        else:
            raise RuntimeError(("mass shooting not closed", widths, residual))
    finally:
        case_context["upper_cabibbo"] = previous_upper
        graph_context["read"] = old_read
    def relative_return(w, v):
        return (w @ (v * [-1, 1, sign])) / (w[0]*v[0])
    assert abs(relative_return(weights, masses)) < 5e-12
    return {"sign": sign, "k": k, "input_widths": widths, "y_d_MZ": d, "iterations": iteration,
            "upper_probabilities": weights.tolist(), "upper_yukawas": masses.tolist(),
            "upper_relative_return": relative_return(weights, masses),
            "relative_return_at_10PeV": relative_return(boundary["uv_weights"], boundary["uv_masses"]),
            "quark_neutrino": output}


def all_phases_positive_sign_box():
    """One-loop comparison on the SM portion from MZ to the upper neutral mass."""
    tau = .028
    scalar_trace = 3*(1.2**2+2*.006**2)+9*.025**2+3*.02**2+.3**2
    log_ratio = 1.5*(1.2**2+.025**2)*tau
    # Bootstrap for sorted Yukawas. Below the upper cluster at most one neutrino
    # column is active. Its insertion norm is bounded by the positive Higgs moment.
    mu = MODEL["spectra"]["neutral_poles_GeV"]
    actual_tau = math.log(mu[2]/FL["MZ"])/(16*math.pi**2)
    assert actual_tau < tau
    norm_insert = math.sqrt(8*math.pi**2*MODEL["spectra"]["neutral_only_Higgs_mass_squared_GeV2"])/mu[0]
    def high_bound(start, own):
        return start*math.exp(tau*(1.5*own**2+scalar_trace))
    tested = [(SPECTRUM["up"][-1][0]+SPECTRUM["up"][-1][1],1.2),
              (sum(SPECTRUM["up"][1]),.006), (sum(SPECTRUM["down"][-1]),.025),
              (sum(SPECTRUM["charged_leptons"][-1]),.02), (norm_insert,.3)]
    assert all(high_bound(start, own)<own for start,own in tested)
    def endpoints(rec):
        return rec[0]-math.hypot(*rec[1:]), rec[0]+math.hypot(*rec[1:])
    d_hi = endpoints(SPECTRUM["down"][0])[1]
    s_lo, s_hi = endpoints(SPECTRUM["down"][1])
    # For d above this independent test cap the final lower bound already holds;
    # below it, use the imposed cap only to control spectral denominators.
    d_cap = .06*s_lo
    down_ratio = max(d_cap/s_lo, s_hi/endpoints(SPECTRUM["down"][2])[0])
    up_ratio = max(endpoints(SPECTRUM["up"][i])[1]/endpoints(SPECTRUM["up"][i+1])[0] for i in (0,1))
    ratio = max(down_ratio, up_ratio)*math.exp(log_ratio)
    sylvester = (1+ratio**2)/(1-ratio**2)
    assert sylvester < 1.01
    coefficient = 1.5*1.01
    w0 = math.hypot(endpoints(DIRECT["Vub"])[1], endpoints(DIRECT["Vcb"])[1])
    w = w0*math.exp(tau*coefficient*(1.2**2+.006**2+2*.025**2))
    assert w < .046
    w = .046
    drift = tau*coefficient*(.025**2+.006**2+1.2**2*w*w/2+(1.2**2+.006**2)*w*w)
    ratio_drift = tau*1.5*(.025**2+.006**2+1.2**2*w*w)
    min_p = (endpoints(DIRECT["Vus"])[0]-drift)**2
    lower = s_lo*min_p/(1-min_p)*math.exp(-ratio_drift)
    assert d_cap > lower > d_hi
    return {"actual_log_interval_over_loop": actual_tau, "trace_cap": scalar_trace,
            "largest_active_neutral_insertion": norm_insert, "max_spectral_gap_multiplier": sylvester,
            "heavy_misalignment_upper_bound": w, "max_drift_Vus": drift,
            "max_drift_log_yd_over_ys": ratio_drift, "positive_branch_necessary_yd_MZ": lower,
            "initial_d_dichotomy_cap": d_cap, "one_width_yd_upper": d_hi}


def assess(row):
    obs = row["quark_neutrino"]
    widths = row["input_widths"]
    cb = DIRECT["Vcb"][0]+widths["Vcb"]*math.hypot(*DIRECT["Vcb"][1:])
    ub = DIRECT["Vub"][0]+widths["Vub"]*math.hypot(*DIRECT["Vub"][1:])
    ratio = ub/cb
    tree = read(CFG["tree_alternatives"])
    row["heldout_yd_widths"] = (row["y_d_MZ"]-SPECTRUM["down"][0][0])/SPECTRUM["down"][0][1]
    row["abs_Vub_over_Vcb"] = ratio
    row["other_tree_ratio_widths"] = {
        key: (ratio-v)/math.hypot(*(up if ratio>=v else down))
        for key, (v, down, up) in tree["ratio_ports"].items()}
    row["other_Bs_phi_widths"] = {
        key: (obs["quarks"]["phi_s_rad"]-raw[0])/math.hypot(*raw[1:])
        for key, raw in tree["phi_s_rad"].items()}
    return row


def main():
    cases = {label: assess(run_point(**params)) for label,params in CFG["points"].items()}
    primary = cases["interior_all_one_width"]
    obs = primary["quark_neutrino"]
    assert abs(primary["heldout_yd_widths"])<1 and max(abs(v) for v in primary["input_widths"].values())<1
    assert all(value for name, value in obs["tests"].items() if name!="all_frequencies_fixed_solar")
    solar = CFG["solar_for_interior_one_width_point_eV2"]
    assert obs["allowed_solar_eV2"][0]<solar<obs["allowed_solar_eV2"][1]
    r2 = obs["mass_ratio_after_transport"]**2
    fixed = cases["central_solar_within_two_widths"]
    assert all(fixed["quark_neutrino"]["tests"].values()) and abs(fixed["heldout_yd_widths"])<1
    fine = run_point(**CFG["points"]["interior_all_one_width"], step=.10)
    errors = {"relative_generated_yd": abs(fine["y_d_MZ"]/primary["y_d_MZ"]-1),
              "gamma_deg": abs(fine["quark_neutrino"]["quarks"]["gamma_deg"]-obs["quarks"]["gamma_deg"]),
              "relative_mass_ratio": abs(fine["quark_neutrino"]["mass_ratio_after_transport"]/math.sqrt(r2)-1)}
    assert max(errors.values())<3e-8
    return {"cases": cases, "all_phase_one_width_box": all_phases_positive_sign_box(),
            "interior_solar_d31_d32_eV2": [solar, r2*solar, (r2-1)*solar], "half_step_differences": errors}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    results = MODEL["sim"]["tidy"](main())
    target = HERE/"results.json"
    if args.check:
        assert json.loads(target.read_bytes())==results
        print("1712 reciprocal mass shooting and full one-width positive-sign exclusion verified")
    else:
        target.write_text(json.dumps(results, indent=2, allow_nan=False)+"\n", encoding="utf-8")
