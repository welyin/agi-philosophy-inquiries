"""Independent moment bounds for real neutral vectors with complex flavour currents."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RESULT = HERE / "independent_results.json"


def trident(u):
    axial = 1.13
    vector = 1 + 4 * 0.231
    return (axial * (1 + u)**2 + (vector + u)**2) / (axial + vector**2)


def inverse_rate(rate):
    # Strict monotonicity on the physical ray gives a unique boundary.
    lo, hi = 0.0, 1.0
    for _ in range(65):
        mid = (lo + hi) / 2
        if trident(mid) > rate:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def alpha_c10(ca):
    # Negative conventional SM C10.  This is alpha times C10, not C10 itself.
    return -2 * math.sqrt(2) * 1.1663787e-5 * 80.379**2 * ca / math.pi


def omega(width):
    return 0.70 * (alpha_c10(.4694 - width * .0036) / (2 * math.pi))**2 / (.00131 + width * .00001)


def leading_log(mass):
    """Integrate the one-loop inverse coupling in each active-flavour interval."""
    a_mz = .118
    inv_b = 1/a_mz + (23/3) * math.log(4.2/91.1876) / (2*math.pi)
    inv_t = 1/a_mz + (23/3) * math.log(163.07/91.1876) / (2*math.pi)
    inv_m = inv_t + 7 * math.log(mass/163.07) / (2*math.pi)
    return (inv_b/inv_t)**(6/23) * (inv_t/inv_m)**(2/7)


def independent_box(width):
    # All experimental entries are primitive magnitudes from the registered sources.
    ratio_hi = (17.757 + width * .021) / (17.59 - width * (1.22 + .78))
    ratio_lo = (17.757 - width * .021) / (17.59 + width * (.33 + .78))
    vus = (.2157 + width * .0006) / (.9704 - width * .0032)
    vts = .04169 + width * .00039
    vcs = math.sqrt(1 - vus**2 - vts**2)
    beta_twice = 2 * math.asin(vus / vcs * (.083 + width * .008))
    phase = .039 + width * (.022 + .006) + .020 + beta_twice
    radial_max = max(abs(x * cmath.exp(1j * phase) - 1) for x in (ratio_lo, ratio_hi))
    u = inverse_rate(.82 + width * .28)
    distance = math.sqrt(u * radial_max / omega(width))
    # Substitute the independent true SM frequency for |lambda_t|^2*f_B^2
    # in the partial width, keeping the hadronic and electroweak factors explicit.
    mu = .1056583755
    sm_upper = (3 / (8*math.pi)) * (alpha_c10(.4694 + width*.0036) * mu / 80.379)**2
    sm_upper *= (1.615 + width*.009) * (17.59 + width*(.33+.78)) * 1e9
    sm_upper *= 1.02 * math.sqrt(1-(2*mu/5.36688)**2)
    sm_upper /= (2.313 - width*.008) * (.5521 - width*.00062) * (1.232 - width*.053)
    min_norm_sq = (3.83 - width * (.36 + .16 + .13)) / sm_upper
    min_distance_left = math.sqrt(1 + min_norm_sq)
    real_distance = math.sqrt(u * (ratio_hi - 1) / omega(width))
    k_upper = 4.97*(vts/.04185)**2*(1 + width*.08)
    nu_upper = k_upper * (2 + (1+.70*distance)**2)/3
    return {
        "oscillation_modulus_range": [ratio_lo, ratio_hi], "max_beta_twice": beta_twice,
        "max_total_phase_with_control": phase, "max_abs_h": radial_max, "u_upper": u,
        "omega_lower": omega(width), "max_abs_P_minus_1": distance,
        "SM_dimuon_upper_in_1e9": sm_upper, "CMS_lower_abs_P_squared": min_norm_sq,
        "left_halfplane_distance_lower": min_distance_left,
        "BK_nunu_upper_in_1e6": nu_upper, "Belle_lower_in_1e6": 4.97*(4.6-width*1.3),
        "left_halfplane_ruled_out_with_clock": distance < min_distance_left,
        "max_abs_P_minus_1_if_real_without_clock": real_distance,
        "real_reflected_branch_ruled_out_without_clock": real_distance < 1 + math.sqrt(min_norm_sq)
    }


def multi_vector(b, p, lepton, weights):
    """Construct two actual flavour couplings. Loop weights multiply the quadratic moment."""
    l = np.array(lepton, dtype=float)
    eta = np.array(weights)
    effective_l = l / np.sqrt(eta)
    norm = np.linalg.norm(effective_l)
    parallel = effective_l / norm
    perpendicular = np.array([-parallel[1], parallel[0]])
    x = p / norm * parallel + cmath.sqrt(b - (p / norm)**2) * perpendicular
    q = x / np.sqrt(eta)
    assert abs(np.dot(q * eta, q) - b) < 2e-14
    assert abs(q @ l - p) < 2e-14
    budget = abs(p)**2 / norm**2 + abs(b - p**2 / norm**2)
    assert abs(np.vdot(x, x).real - budget) < 2e-14
    return {"loop_weights": list(weights), "lepton": list(lepton),
            "flavour_complex": [[a.real, a.imag] for a in q],
            "sharp_weighted_positive_budget": float(budget),
            "sum_lepton_squares": float(l @ l)}


def run():
    boxes = {str(w): independent_box(w) for w in (1, 2)}
    # Start from Lagrangian matching, with v*g_q/(M*lambda_t) and v*g_l/M.
    ca = alpha_c10(.4694)
    eta, loop, u = leading_log(1000.), .00131, .272
    q = -2 * ca / (math.pi * math.sqrt(u))
    h = eta * q**2 / (4 * loop)
    z = math.pi * q * math.sqrt(u) / ca
    assert abs(z + 2) < 1e-14 and abs(trident(u) - 1.38) < .02
    mix_required = 17.757 / (1 + h)
    assert 17.59 - 2 * (1.22 + .78) <= mix_required <= 17.59 + 2 * (.33 + .78)
    # Different raw amplitudes, identical |1+P-1|=1 and the same trident and mixing.
    # Both take the same physical positive lepton norm and two inequivalent thresholds.
    l = (math.sqrt(.03), math.sqrt(.03))
    weights = (leading_log(1000.), leading_log(2000.))
    same_branch = multi_vector(.00005, 0., l, weights)
    other_branch = multi_vector(.00005, -2 * ca / math.pi, l, weights)
    assert boxes["1"]["left_halfplane_ruled_out_with_clock"]
    assert not boxes["2"]["real_reflected_branch_ruled_out_without_clock"]
    gf, mw, m_bs, f_bs, bag = 1.1663787e-5, 80.379, 5.36688, .231, 1.127
    sm_kernel = gf**2 * mw**2 * m_bs * 2.313 * .5521 / (6*math.pi**2)
    ckm_sq = mix_required * 6.582119569e-13 / (sm_kernel * f_bs**2 * bag)
    mu = .1056583755
    # Decay amplitude before lifetime conversion, evaluated with the very same CKM.
    width = gf**2 / (16*math.pi**3) * alpha_c10(.4694)**2 * ckm_sq
    width *= f_bs**2 * m_bs * mu**2 * math.sqrt(1 - (2*mu/m_bs)**2)
    dimuon = width * 1.615 / 6.582119569e-13 * 1e9
    b_to_k = 4.97 * ckm_sq / .04185**2
    nu_ratio = (2 + (1 + 2*.64)**2)/3
    b_to_k_np = b_to_k * nu_ratio
    assert 3.38 < dimuon < 3.94
    assert 4.97*(4.6-2*1.3) < b_to_k_np < 4.97*(4.6+2*1.3)
    return {"CCFR_raw_yield_ratio": 37 / 45.3,
            "CCFR_error_independent_propagation": math.hypot(12.4/45.3, 37*2.3/45.3**2),
            "boxes": boxes,
            "one_vector_at_1_TeV": {"u": u, "P": 1 + z,
                 "h": h, "trident": trident(u), "SM_frequency_required": mix_required,
                 "true_abs_CKM_product": math.sqrt(ckm_sq), "SM_Bs_in_1e9": dimuon,
                 "rescaled_SM_BK_nunu_in_1e6": b_to_k, "BK_nunu_in_1e6": b_to_k_np},
            "two_equal_rate_branches": {"masses_GeV": [1000, 2000],
                 "h": .00005/(4*loop), "trident": trident(.06),
                 "P_plus_one": same_branch, "P_minus_one": other_branch}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    flags = parser.add_mutually_exclusive_group()
    flags.add_argument("--write", action="store_true")
    flags.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        RESULT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    else:
        expected = json.loads(RESULT.read_text(encoding="utf-8"))
        def compare(actual, saved):
            if isinstance(actual, dict):
                assert actual.keys() == saved.keys()
                for key in actual:
                    compare(actual[key], saved[key])
            elif isinstance(actual, list):
                assert len(actual) == len(saved)
                for a, b in zip(actual, saved):
                    compare(a, b)
            elif isinstance(actual, bool):
                assert actual is saved
            else:
                assert math.isclose(actual, saved, rel_tol=2e-13, abs_tol=2e-15)
        compare(result, expected)
        print("1407 matching, phase boxes and sharp multi-vector cost verified")
