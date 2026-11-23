"""Independent one-loop quadratures and a same-mass Higgs portal test."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def nodes(count):
    point, weight = np.polynomial.legendre.leggauss(count)
    return (1+point)/2, weight/2


def kernels(fermion, light_sq, heavy_sq, q_sq, muon, count):
    # Flatten the narrow t=1 endpoint with t=1-x^4.
    x, wx = nodes(count)
    t = 1-x**4
    wt = 4*x**3*wx
    mass_gap = heavy_sq-light_sq
    low = x**4*fermion**2 + t*light_sq
    high = x**4*fermion**2 + t*heavy_sq
    triangle0 = np.sum(wt*np.log1p(t*mass_gap/low))/mass_gap
    dipoles = []
    for charge in (0, 1):
        factor = t*x**4*(charge-x**4)/(low*high)
        dipoles.append(float(2*muon**2*np.dot(wt, factor)/triangle0))
    # Attach the Higgs to each scalar eigenstate. External muon mass set to zero.
    s, ws = nodes(count)
    momentum = q_sq * t[:, None]**2 * (s*(1-s))[None, :]
    triangle_pair = [float(np.dot(wt*t, (1/(a[:, None]-momentum)) @ ws))
                     for a in (low, high)]
    return float(triangle0), dipoles, triangle_pair


def scalar_yukawa(diagonal, mixing, fermion, muon, higgs, vev, count):
    low, high = diagonal-mixing, diagonal+mixing
    norm, dipoles, triangles = kernels(fermion, low, high, higgs**2, muon, count)
    lo, hi = triangles
    k0, slope = (lo+hi)/(2*norm), vev**2*(hi-lo)/(2*mixing*norm)
    yuk_product = muon*16*math.pi**2/(fermion*mixing*norm)
    # Independently use the charged Dirac fermion/neutral complex scalar F1C and F2C.
    external, internal = 0., 0.
    for scalar_sq, sign in [(low, -1), (high, 1)]:
        ratio = fermion**2/scalar_sq
        f1 = 2*(2+3*ratio-6*ratio**2+ratio**3+6*ratio*math.log(ratio))/(1-ratio)**4
        f2 = -1.5*(3-4*ratio+ratio**2+2*math.log(ratio))/(1-ratio)**3
        external += muon**2*yuk_product*f1/(192*math.pi**2*scalar_sq)
        internal += sign*muon*fermion*yuk_product*f2/(48*math.pi**2*scalar_sq)
    assert math.isclose(internal, dipoles[0], rel_tol=2e-11)
    witness = []
    for lam in (.5, 2.5):
        bare = diagonal - lam*vev**2/2
        global_floor = bare - (mixing/vev)**2/(2*lam)
        kappa = k0+lam*slope
        rate = kappa**2
        assert global_floor > 0
        witness.append({"lambda": lam, "bare_scalar_squared": bare,
                        "global_extra_scalar_squared_floor": global_floor,
                        "kappa_on_shell": kappa, "signal_strength_with_transfer": [.75*rate, 1.25*rate],
                        "minimum_rate_transfer_for_experimental_lower_edge": .6/rate})
    assert witness[0]["signal_strength_with_transfer"][0] > .6
    assert witness[0]["signal_strength_with_transfer"][1] < 2.2
    assert witness[1]["signal_strength_with_transfer"][1] < .6
    return {"scalar_masses": [math.sqrt(low), math.sqrt(high)], "fermion_mass": fermion,
            "matched_yukawa_product": yuk_product, "internal_flip_Q0": internal,
            "external_flip_equal_yukawa_strengths_Q0": external,
            "soft_or_on_shell_intercept": k0, "slope_per_portal_coupling": slope,
            "two_full_portal_choices": witness}


def calculate():
    data = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    muon = data["muon_mass_GeV"]
    bench = []
    for sample in data["benchmarks"]:
        fermion, a, b = sample["masses_GeV"]
        dipoles = kernels(fermion, a*a, b*b, 0, muon, 120)[1]
        bench.append([round(v*1e11,7) for v in dipoles])
    arg = (300.**2,.7*300.**2,5000.,muon)
    on = scalar_yukawa(*arg,125.2,data["v_GeV"],96)
    check = scalar_yukawa(*arg,125.2,data["v_GeV"],160)
    for key in ("slope_per_portal_coupling", "soft_or_on_shell_intercept", "internal_flip_Q0"):
        assert math.isclose(on[key],check[key],rel_tol=1e-11)
    soft = kernels(5000.,27000.,153000.,0,muon,96)
    soft_kappa = sum(soft[2])/(2*soft[0])
    assert 1 < soft_kappa < on["soft_or_on_shell_intercept"]
    exp, wp20 = data["experiment_world"], data["WP20_theory"]
    historical_lower = exp[0]-wp20[0]-2*math.hypot(exp[1],wp20[1])
    same_a = (on["internal_flip_Q0"]+on["external_flip_equal_yukawa_strengths_Q0"])*1e11
    assert -88 < same_a < 164 and same_a < historical_lower
    return {"authors_benchmarks_Q0_Q1_at_1e_minus_11": bench, "higgs_dimuon": on,
            "zero_portal_soft_kappa": soft_kappa, "higgs_signal_strength_work_interval": [.6,2.2],
            "full_one_loop_equal_yukawas_Q0_in_1e_minus_11": same_a,
            "WP20_with_latest_experiment_two_width_lower": historical_lower}


def equivalent(left,right):
    if isinstance(left,dict):
        return left.keys()==right.keys() and all(equivalent(left[k],right[k]) for k in left)
    if isinstance(left,list):
        return len(left)==len(right) and all(equivalent(a,b) for a,b in zip(left,right))
    return math.isclose(left,right,rel_tol=1e-10,abs_tol=1e-23)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=calculate()
    dest=HERE/"independent_results.json"
    if args.write_new:
        with dest.open("x",encoding="utf-8") as handle:
            json.dump(result,handle,indent=2)
            handle.write("\n")
    else:
        assert equivalent(result,json.loads(dest.read_text(encoding="utf-8")))
    print(json.dumps(result))
