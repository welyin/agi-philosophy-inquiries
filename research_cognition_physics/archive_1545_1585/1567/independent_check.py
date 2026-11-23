"""Independent parity and dispersive-phase checks; no fit to rare branching fractions."""
import argparse
import cmath
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"


def read(relative):
    return json.loads((HERE / relative).read_text(encoding="utf-8-sig"))


def box(data, widths):
    center, *errors = data
    radius = widths * sum(errors)
    return center - radius, center + radius


def trig_range(lo, hi, cosine=False):
    func = math.cos if cosine else math.sin
    shift = 0 if cosine else math.pi / 2
    candidates = [func(lo), func(hi)]
    for k in range(math.ceil((lo - shift) / math.pi), math.floor((hi - shift) / math.pi) + 1):
        candidates.append(func(shift + k * math.pi))
    return min(candidates), max(candidates)


def mixing_bound(widths, mirror, weak_source, charged, mixing):
    usnum = box(weak_source["K_product"], widths)
    usden = box(weak_source["K_form"], widths)
    us = usnum[0] / usden[1], usnum[1] / usden[0]
    cb = box(charged["BelleII_Vcb_BGL"], widths)
    rho = box(weak_source["ub_over_cb"], widths)
    ub_hi = cb[1] * rho[1]
    c13_min = math.sqrt(1 - ub_hi ** 2)
    s12_min, s12_max = us[0], us[1] / c13_min
    c12_min = math.sqrt(1 - s12_max ** 2)
    s23_max = cb[1] / c13_min
    c23_min = math.sqrt(1 - s23_max ** 2)
    cs_min = c12_min * c23_min - s12_max * s23_max * ub_hi
    rmin, rmax = us[0] * rho[0], us[1] * rho[1] / cs_min
    phase_transfer = (math.asin(s23_max * ub_hi / (s12_min * c23_min)) +
                      math.asin(s12_max * s23_max * ub_hi / (c12_min * c23_min)))
    gamma = weak_source["gamma_conditional_deg"][0]
    sign = -1 if mirror else 1
    tau_min = math.radians(gamma - widths * 5.1 - 1.0) - phase_transfer
    tau_max = math.radians(gamma + widths * 5.2 + 1.0) + phase_transfer
    si = trig_range(tau_min, tau_max)
    co = trig_range(tau_min, tau_max, cosine=True)
    den = [1 + sign * r * c for r, c in itertools.product((rmin, rmax), co)]
    a_low = 2 * math.atan(rmin * si[0] / max(den))
    a_high = 2 * math.atan(rmax * si[1] / min(den))
    top = (a_low, a_high) if mirror else (-a_high, -a_low)
    # Explicit matching budget for the baseline mixing phase (not from data).
    baseline_phase = (top[0] - .010, top[1] + .010)
    phi_data = box(mixing["phase_inputs_rad"]["LHCb_KK_2015_18"], widths)
    phi_m = (phi_data[0] - .020, phi_data[1] + .020)
    # H direction -2*tau modulo 2*pi on both experimental DK branches.
    sep = trig_range(2 * tau_min + baseline_phase[0], 2 * tau_max + baseline_phase[1])
    assert sep[0] > 0
    numer = trig_range(phi_m[0] - baseline_phase[1], phi_m[1] - baseline_phase[0])
    fraction = max(abs(numer[0]), abs(numer[1])) / sep[0]
    freq_upper = box(mixing["omega_ps_inverse"], widths)[1]
    # Additional 0.1% conversion budget, e.g. any difference Delta m / (2 |M12|).
    freq_cap = fraction * freq_upper / .999
    return {"abs_lambda_u_over_c": [rmin, rmax], "phase_translation_deg": math.degrees(phase_transfer),
            "tau_mod_pi_deg": [math.degrees(tau_min), math.degrees(tau_max)],
            "baseline_phase_rad": baseline_phase, "physical_phase_rad": phi_m,
            "min_transverse_sine": sep[0], "max_fraction": fraction,
            "max_twice_scalar_mixing_ps_inverse": freq_cap}


def unit_tests():
    algebra_error = 0.0
    # Nontrivial phases, complex CKM factor, all four allowed spectral sign patterns.
    for eb, es, angle in itertools.product((-1, 1), (-1, 1), (0., .21, 1.3, math.pi / 2)):
        ed = -es
        c = cmath.exp(1j * angle)
        for eq in (es, ed):
            x = (.8 * eb + .06 * eq) * (.03 + .02j)
            right, left = x * eb * c, x * eq * c.conjugate()
            for same, crossed in ((1., -.27), (-2., .12), (0., 1.)):
                direct = (same * (right**2 + left**2) + 2 * crossed * right * left) / x**2
                reduced = 2 * (same * (c**2).real + crossed * eb * eq * abs(c)**2)
                algebra_error = max(algebra_error, abs(direct - reduced))
        algebra_error = max(algebra_error, abs((eb*c + es*c.conjugate()) -
                                                (eb*c - ed*c.conjugate())))
    assert algebra_error < 1e-13
    # Independent real coefficients recoverable by projected complex-plane area.
    ptop, ph, ht, hh = -.034, -2.3, 1.14, -.3
    total = ht * cmath.exp(1j * ptop) + hh * cmath.exp(1j * ph)
    recovered = abs(math.sin(cmath.phase(total) - ptop) / math.sin(ph - ptop))
    assert abs(recovered - abs(hh / total)) < 1e-14
    return {"all_parity_and_sign_cases": True, "signed_two_ray_reconstruction": True}


def compute():
    weak = read("../1381/inputs.json")
    charged = read("../1412/inputs.json")
    mixing = read("../1468/inputs.json")
    branches = {}
    for width in (1, 2):
        for mirror in (False, True):
            label = "gamma_plus_pi" if mirror else "gamma"
            branches[f"{width}_width_{label}"] = mixing_bound(width, mirror, weak, charged, mixing)
    return {"tests": unit_tests(), "inputs": {"decay_phase_budget_rad": .020,
            "baseline_phase_budget_rad": .010, "mass_clock_relative_budget": .001,
            "gamma_measurement_deg": [68.7, -5.1, 5.2], "gamma_transport_deg": 1.,
            "phase_data": mixing["phase_inputs_rad"]["LHCb_KK_2015_18"],
            "clock_ps_inverse": mixing["omega_ps_inverse"]}, "bounds": branches}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    command = parser.add_mutually_exclusive_group(required=True)
    command.add_argument("--check", action="store_true")
    command.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    data = json.loads(json.dumps(compute()))
    if args.write_new:
        with OUT.open("x", encoding="utf-8") as dst:
            json.dump(data, dst, ensure_ascii=False, indent=2)
            dst.write("\n")
    else:
        if data != read("independent_results.json"):
            raise SystemExit("independent response changed")
        print("1567 independent mixing and parity identities passed")
