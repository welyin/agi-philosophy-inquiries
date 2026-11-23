"""Neutral kaon cross-check independent of the charged-B implementation."""
import argparse
import cmath
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "independent_results.json"


def phase_space(kaon, pion):
    return (kaon * (1.0 - (pion / kaon) ** 2)) ** 3


def pion_momentum(kaon, pion, invisible):
    return math.sqrt((kaon**2 - (pion + invisible)**2) *
                     (kaon**2 - (pion - invisible)**2)) / (2 * kaon)


def phase_envelope(charge_asymmetry):
    cosine = math.sqrt(1 - charge_asymmetry**2)
    return charge_asymmetry**2 / (2 * (1 + cosine)), (1 + cosine) / 2


def row_phase_bounds(box):
    us_min, us_max = box["Vus"]
    cb_max = box["Vcb"][1]
    s13 = cb_max * box["ratio_ub_cb"][1]
    c13 = math.sqrt(1 - s13**2)
    s12, s23 = us_max / c13, cb_max / c13
    c12, c23 = math.sqrt(1 - s12**2), math.sqrt(1 - s23**2)
    lu_min = us_min * c12 * c13
    lc_min = (us_min * c23 - s23 * s13) * (c12 * c23 - s12 * s23 * s13)
    lt_max = (s12 * s23 + s13) * (s23 + s12 * s13)
    return {"u_upper": 0.0, "c_upper": s12*s23*s13/(lu_min*lc_min),
            "t_lower": box["absolute_J_lower"]/(us_max*lt_max)}


def calculate():
    input_box = json.loads((HERE.parent / "1598/results.json").read_text(encoding="utf-8"))["two_width_box"]
    # KLOE independent lifetimes in ns; masses in MeV, broad transfer boxes.
    tl, tc = 50.92, 12.347
    h_min = ((tl - 2*(0.17+0.25)) / (tc + 2*0.030) *
             phase_space(497.60, 135.00) / phase_space(493.71, 139.55) * 0.95**2)
    h_max = ((tl + 2*(0.17+0.25)) / (tc - 2*0.030) *
             phase_space(497.63, 134.96) / phase_space(493.65, 139.60) * 1.05**2)
    d_min = (3322 - 2*(58+47))*1e-6
    lower, upper = phase_envelope(d_min)
    conjugate_ratio = math.sqrt((1+d_min)/(1-d_min))
    xi_box = math.atan((6.98+2*(0.62+1.44))/(2.99-2*(0.32+0.59))*1e-4)
    xi, eps = 0.001, 0.003
    phase = row_phase_bounds(input_box)
    upper_u = h_max * (phase["u_upper"] + xi + eps)**2
    upper_c = h_max * (phase["c_upper"] + xi + eps)**2
    lower_t = h_min * (phase["t_lower"] - xi - eps)**2 / (1 + eps**2)
    # Errors are fractions of the original charged or flavor-tagged neutral
    # amplitude, before projecting the two neutral sources onto K_L.
    reference_c = phase["c_upper"] + xi + eps
    reference_t = (phase["t_lower"] - xi - eps) / math.sqrt(1 + eps**2)

    def robust_bounds(fraction):
        top_min = h_min * (max(0, reference_t-fraction)/(1+fraction))**2
        charm_max = h_max * ((reference_c+fraction)/(1-fraction))**2
        return top_min, charm_max

    left, right = 0.0, reference_t
    for _ in range(60):
        middle = (left + right) / 2
        top_min, charm_max = robust_bounds(middle)
        if top_min > charm_max:
            left = middle
        else:
            right = middle
    # Validate the two amplitude forms for a complex epsilon and arbitrary phases.
    e = 0.0022 * cmath.exp(0.76j)
    norm = math.sqrt(2*(1 + abs(e)**2))
    p, q = (1+e)/norm, (1-e)/norm
    for phi in (-1.7, -0.2, 0.0, 0.81):
        direct = abs(p*cmath.exp(1j*phi) - q*cmath.exp(-1j*phi))**2 / 2
        reduced = abs(1j*math.sin(phi) + e*math.cos(phi))**2/(1+abs(e)**2)
        assert abs(direct-reduced) < 5e-16
    assert 0 < xi_box < xi
    assert 0 < lower_t < h_max*upper and upper_c < lower_t
    return {
        "isospin_and_clock_H": [h_min, h_max],
        "KTeV_charge_asymmetry_lower": d_min,
        "neutral_over_charged_unrestricted_weak_phase": [h_min*lower, h_max*upper],
        "minimum_conjugate_amplitude_ratio_for_exact_blindness": conjugate_ratio,
        "lattice_SM_two_width_max_abs_A0_phase_radians": xi_box,
        "registered_phase_allowances_radians_and_epsilon_modulus": [xi, eps],
        "lambda_row_relative_to_lambda_u_sines": phase,
        "neutral_over_charged_phase_referenced": {
            "u_upper": upper_u, "c_upper": upper_c, "t_lower": lower_t,
            "t_lower_over_c_upper": lower_t/upper_c},
        "four_percent_amplitude_error_bounds_t_lower_and_c_upper": list(robust_bounds(0.04)),
        "separation_certificate_error_fraction": left,
        "KOTO_zero_to_1MeV_pion_momentum_shift_MeV":
            pion_momentum(497.611,134.977,0)-pion_momentum(497.611,134.977,1)}


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    options = parser.parse_args()
    result = calculate()
    if options.write:
        with OUTPUT.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2, ensure_ascii=False)
            stream.write("\n")
    elif options.check:
        expected = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if result != expected:
            raise SystemExit("independent results mismatch")
        print("1599 independent check passed")
    else:
        print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
