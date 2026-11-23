"""Round 1532. Evaluate the electron-dipole kernels and project the shared direction."""
import argparse
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent


@lru_cache(maxsize=None)
def quadrature(order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return (nodes + 1) / 2, weights / 2


def photon_kernel(mass, boson_mass, order=140):
    # Symmetry x <-> 1-x; x=t^4/2 removes the integrable endpoint logarithm.
    t, weights = quadrature(order)
    x = t ** 4 / 2
    ratio = x * (1 - x) * (boson_mass / mass) ** 2
    delta = ratio - 1
    near = abs(delta) < 1e-7
    kernel = np.ones_like(delta)
    kernel[near] -= delta[near] / 2
    kernel[~near] = np.log(ratio[~near]) / delta[~near]
    return float(np.dot(weights, 4 * t ** 3 * kernel))


def loop(name, benchmark, order=140):
    top = name == "top"
    color, charge, isospin = (3, 2 / 3, .5) if top else (1, -1, -.5)
    mass = benchmark["top_MS_GeV" if top else "tau_GeV"]
    h, z, s = (benchmark[k] for k in ("higgs_GeV", "z_GeV", "sin2_weak"))
    f_h = photon_kernel(mass, h, order)
    f_z = photon_kernel(mass, z, order)
    weak_kernel_abs = (h*h*f_h - z*z*f_z) / (h*h-z*z)
    v_over_q_product = abs((-.5 + 2*s) * (isospin - 2*charge*s)
                           / (-charge * 4*s*(1-s)))
    assert 0 < weak_kernel_abs < f_h
    prefactor = (4 * math.sqrt(2) * benchmark["fermi_GeV_minus2"]
                 * benchmark["electron_GeV"] * benchmark["hbar_c_GeV_cm"]
                 / (benchmark["alpha_inverse"] * (4*math.pi)**3))
    scale = prefactor * color * charge**2
    return {"photon_e_cm": scale*f_h, "Z_absolute_e_cm": scale*v_over_q_product*weak_kernel_abs,
            "photon_kernel": f_h, "Z_absolute_kernel": weak_kernel_abs}


def calculate(data):
    ref, boxes = data["benchmark"], data["diagnostic_box"]
    mol = data["HfF"]
    cap = mol["upper_abs_e_cm_90"]
    loops = {name: loop(name, ref) for name in ("top", "tau")}
    for name, value in loops.items():
        for k, v in loop(name, ref, 92).items():
            assert math.isclose(v, value[k], rel_tol=1e-10)
    def edges(item):
        return (item["photon_e_cm"]-item["Z_absolute_e_cm"],
                item["photon_e_cm"]+item["Z_absolute_e_cm"])
    top, tau = edges(loops["top"]), edges(loops["tau"])
    low, high = [sum(pair) for pair in zip(top, tau)]
    ke_low = boxes["abs_kappa_e"][0]
    k_tau_low = boxes["kappa_tau"][0]
    sl, sh = boxes["sin2_weak"]
    min_ref = {**ref, "top_MS_GeV": boxes["top_MS_GeV"][0],
               "tau_GeV": boxes["tau_GeV"][0], "higgs_GeV": boxes["higgs_GeV"][1]}
    # Uniform absolute bound for |ve vf /(Qe Qf)|. Kernel_Z < kernel_photon.
    max_ev = max(abs(-.5+2*sl), abs(-.5+2*sh))
    guaranteed = 0
    for name, q, t3 in [("top", 2/3, .5), ("tau", -1, -.5)]:
        vf = max(abs(t3-2*q*sl), abs(t3-2*q*sh))
        frac = max_ev * vf / (4*abs(q)*sl*(1-sl))
        guaranteed += loop(name, min_ref)["photon_e_cm"]*(1-frac)
    guaranteed *= ke_low*(1-boxes["relative_higher_order_budget"])
    nominal_angle = math.degrees(math.atan(cap/low))
    blind = [-tau[1]/top[0], -tau[0]/top[1]]
    total = {"benchmark_kernels": loops, "shared_unit_response_abs_e_cm": [low, high],
             "HfF_nominal_abs_eta_upper": cap/low, "nominal_alpha_degrees_upper": nominal_angle,
             "two_molecule_projected_abs_angle_if_only_de_CS":
                 math.degrees(math.atan(mol["joint_ThO_bound_on_de"]/low)),
             "stress_box_abs_response_lower_e_cm": guaranteed,
             "stress_box_alpha_degrees_upper": math.degrees(math.atan(cap/(guaranteed*k_tau_low))),
             "independent_top_per_tau_blind_line_slope_envelope": blind,
             "CMS_combined_68_interval_degrees": [-9,23], "CMS_Run3_68_interval_degrees": [6,69],
             "cancellation_at_registered_angles": {}}
    for angle in mol["Tau_angle_tests_degrees"]:
        eta = math.tan(math.radians(angle))
        total["cancellation_at_registered_angles"][str(angle)] = {
            "linearized_eta_for_unit_even_tau": eta,
            "min_extra_dipole_e_cm_at_central": low*eta-cap,
            "min_fraction_cancelled_central": 1-cap/(low*eta),
            "abs_electron_scalar_needed_without_other_source_at_central": cap/(low*eta)}
    for label, obs in [("CMS_combined_68_interval_degrees", data["CMS"]["Run2_and_Run3"]),
                       ("CMS_Run3_68_interval_degrees", data["CMS"]["Run3_only"])]:
        assert total[label] == [obs["central"]-obs["minus"], obs["central"]+obs["plus"]]
    assert 0 < cap/low < .001 and blind[0] < blind[1] < 0
    assert nominal_angle < total["stress_box_alpha_degrees_upper"] < .1
    return total


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = calculate(json.loads((BASE/"inputs.json").read_text(encoding="utf-8-sig")))
    path = BASE/"results.json"
    if args.write_new:
        with path.open("x", encoding="utf-8") as output:
            json.dump(result, output, ensure_ascii=False, indent=2, allow_nan=False)
            output.write("\n")
    else:
        recorded = json.loads(path.read_text(encoding="utf-8-sig"))
        def same(a, b):
            if isinstance(a, dict):
                return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
            if isinstance(a, list):
                return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
            return math.isclose(a, b, rel_tol=1e-10, abs_tol=0)
        assert same(result, recorded)
        print("1532: quadrature, photon/Z bounds, experimental readouts and cancellation verified")


if __name__ == "__main__":
    main()
