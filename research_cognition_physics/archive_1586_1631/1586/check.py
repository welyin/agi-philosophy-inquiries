"""BHL inverse-Yukawa integral and finite-window rate comparison."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
K = 8 * math.pi**2
A = 4.5


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def integration(mu, g, b, c, mz, n, matter_cap=0.0):
    t = np.linspace(0, math.log(mu/mz), n + 1)
    step = float(t[1])
    log_kernel = matter_cap * t / K
    for g_i, b_i, c_i in zip(g, b, c):
        log_kernel += (c_i / b_i) * np.log1p(-b_i * g_i**2 * t / K)
    kernel = np.exp(log_kernel)
    assert np.all(np.diff(kernel) < 0)
    # This empirical check is complemented by the analytical C>D bound in proof.md.
    # Left rectangle sums are upper bounds. Expose a conservative roundoff guard.
    upper = np.cumsum(kernel[:-1]) * step + 1e-9
    lower = float(step * np.sum(kernel[1:]) - 1e-9)
    assert lower < upper[-1]
    return t, kernel, upper, lower


def calculate():
    spec = read("inputs.json")
    g_info = read(spec["gauge_source"])
    yuk = read(spec["yukawa_source"])["yukawa"]
    lamb = read(spec["quartic_source"])
    g_cent, g_error = map(np.array, zip(*(g_info[key] for key in ("g1_GUT", "g2", "g3"))))
    mz = g_info["mu_ref_GeV"]
    n = spec["n_rectangles"]
    widths = spec["box_widths"]
    ymax = yuk["up"][-1][0] + widths*yuk["up"][-1][1]
    umin = 1/ymax**2
    gmin = g_cent - widths*g_error
    c = spec["gauge_c_top"]
    b = spec["gauge_b_SM"]
    masses = spec["windows_GeV"]
    tmax = math.log(masses[-1]/mz)
    samples = {}
    for mu in masses:
        t, kernel, high, low = integration(mu, g_cent, b, c, mz, n)
        middle = (high[-1] + low)/2
        j = A*middle/K
        top = yuk["up"][-1][0]
        samples[str(mu)] = dict(inverse_IR_for_pole=j, y_IR_for_pole=1/math.sqrt(j),
                                inverse_UV_for_real_top=(1/top**2-j)/kernel[-1],
                                spread_in_j=A*(high[-1]-low)/K)
    _, _, high, _ = integration(masses[-1], gmin, b, c, mz, n)
    pole_bound = A*high[-1]/K
    cap = spec["bootstrap_caps"]
    dcap = 6*cap["other_up"]**2+9*cap["down"]**2+3*cap["lepton"]**2
    t, kernel, high, _ = integration(masses[-1], gmin, b, c, mz, n, dcap)
    lower_u = (umin - A*high/K) / kernel[:-1]
    assert np.min(lower_u)>1
    # In each sector gauge and opposite-sector terms only decrease growth of a singular value.
    max_rate = 3 + dcap + 1.5*max(cap[k]**2 for k in ("other_up","down","lepton"))
    growth = math.exp(max_rate*tmax/(2*K))
    small_final = {}
    for name, capkey in (("up","other_up"),("down","down"),("charged_leptons","lepton")):
        small = yuk[name][:-1] if name=="up" else yuk[name]
        small_final[name] = max(x+widths*e for x,e in small)*growth
        assert small_final[name]<cap[capkey]
    # Include the entire mesh interval which straddles 1 TeV.
    window = t[1:] >= math.log(masses[0]/mz)
    cap_top_pole = A*high[-1]/K
    assert umin > cap_top_pole
    general = dict(weak_top_max=ymax, inverse_weak_min=umin,
                   top_only_pole_y_lower=1/math.sqrt(pole_bound),
                   full_sm_pole_y_lower_under_bootstrap_caps=1/math.sqrt(cap_top_pole),
                   top_only_needed_integral_positive_D=K*math.log(umin/pole_bound),
                   needed_sup_extra_a_top_only=A*(umin/pole_bound-1),
                   needed_integral_extra_D_after_ordinary_cap=K*math.log(umin/pole_bound)-dcap*tmax,
                   needed_sup_extra_a_after_ordinary_cap=A*(umin/cap_top_pole-1),
                   extra_a_constant_at_central_max_scale=A*(1/(yuk["up"][-1][0]**2*
                          samples[str(masses[-1])]["inverse_IR_for_pole"])-1),
                   finite_boundary_inverse_lower=float(np.min(lower_u[window])),
                   all_scales_inverse_lower=float(np.min(lower_u)),
                   matter_pointwise_cap=dcap,matter_integral_cap=dcap*tmax,
                   small_flavour_end_bounds=small_final,
                   scalar_lambda_MZ_box=[(lamb["lambda_paper"][0]+z*lamb["lambda_paper"][1])/
                                           lamb["lambda_paper_factor"] for z in (-widths,widths)])
    alt = spec["two_doublet"]
    _, _, high, _ = integration(masses[-1], gmin, spec["gauge_b_two_doublets"], c, mz, n)
    i_alt = float(high[-1])*math.exp(alt["other_yukawa_signed_response_integral_upper"]/K)
    ycomp = math.sqrt(K/(A*i_alt))
    proj = ymax/ycomp
    gtb = alt["Vtb_floor"]*math.sqrt(ycomp**2-ymax**2)
    mt, mh = alt["mt_floor_GeV"], alt["charged_mass_test_GeV"]
    gamma = (mt/(32*math.pi)*(gtb*alt["vertex_transport_amplitude_floor"])**2*
             (1-(mh/mt)**2)**2*alt["finite_width_transfer_floor"])
    bt = gamma/(alt["other_top_width_cap_GeV"]+gamma)
    product = alt["ATLAS_observed_upper_product_envelope"]
    competing = dict(y_comp_lower=ycomp, vev_share_upper=proj, vev_squared_share_upper=proj**2,
                     charged_yukawa_floor_at_MZ=gtb,width_t_to_Hb_floor_GeV=gamma,
                     branching_top_floor=bt,
                     conditional_H_to_cb_upper=product/(bt*alt["effective_signal_transfer_floor"]))
    return {"central_sample_not_interval_proof":samples,"whole_window":general,
            "second_doublet":competing}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    check = parser.parse_args().check
    result = calculate()
    destination = HERE / "results.json"
    if check:
        actual = read("results.json")
        def compare(left, right):
            if isinstance(left, dict):
                return set(left) == set(right) and all(compare(left[k],right[k]) for k in left)
            if isinstance(left, list):
                return len(left)==len(right) and all(compare(a,b) for a,b in zip(left,right))
            return math.isclose(left,right,rel_tol=5e-10,abs_tol=1e-13)
        if not compare(actual, result):
            raise SystemExit("results differ")
        print("1586: continuous bounds, ordinary families and charged top branch verified")
    else:
        destination.write_text(json.dumps(result, indent=2)+"\n",encoding="utf-8")
