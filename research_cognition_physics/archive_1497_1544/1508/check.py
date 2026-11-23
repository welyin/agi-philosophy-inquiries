"""Check anomaly sums, alloy response and the normalized central colour weight."""

import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def charges(q, ell, include_nu):
    # Every row is a left-handed Weyl fermion, including conjugate right fields.
    rows = [(6, F(1, 6), q), (3, F(-2, 3), -q), (3, F(1, 3), -q),
            (2, F(-1, 2), ell), (1, F(1), -ell)]
    if include_nu:
        rows.append((1, F(0), -ell))
    return {
        "SU3_squared_X": 2*q-q-q,
        "SU2_squared_X": (3*q+ell)/2,
        "Y_squared_X": sum(n*y*y*x for n, y, x in rows),
        "Y_X_squared": sum(n*y*x*x for n, y, x in rows),
        "X_cubed": sum(n*x**3 for n, y, x in rows),
        "gravity_X": sum(n*x for n, y, x in rows),
        "sum_X_squared": sum(n*x*x for n, y, x in rows),
        "sum_Y_X": sum(n*y*x for n, y, x in rows)
    }


def certificate():
    a = read(HERE / "inputs.json")
    static_limit = a["published_gBL_95"]
    n = a["neutrons_per_atomic_mass_unit"]
    def alloy(w):
        assert math.isclose(sum(w.values()), 1.0)
        return sum(n[name]*weight for name, weight in w.items())
    nt = alloy(a["titanium_mass_fractions"])
    np_ = alloy(a["platinum_mass_fractions"])
    assert abs(nt-.54043) < 1e-5 and abs(np_-.59668) < 1e-5 and np_ > nt
    e_squared = 4*math.pi/a["inverse_alpha_at_zero_rounded"]
    response = a["eta_per_epsilon_squared_massless"] / e_squared
    reconstructed = math.sqrt(a["physical_vector_positive_eta_95"]/response)
    assert 1.06e-25 < reconstructed < 1.08e-25 and reconstructed < static_limit

    full = charges(F(1, 3), F(-1), True)
    minimal = charges(F(1, 3), F(-1), False)
    assert all(full[k] == 0 for k in list(full)[:6])
    assert minimal["X_cubed"] == minimal["gravity_X"] == -1
    for q, ell in [(F(3, 7), F(8, 9)), (F(-1), F(2))]:
        for right in [True, False]:
            c = charges(q, ell, right)
            assert c["Y_squared_X"] == -c["SU2_squared_X"]
            assert c["Y_X_squared"] == 0
            assert c["X_cubed"] == (0 if right else ell**3)
    bxx = 3*F(2, 3)*full["sum_X_squared"]
    byx = 3*F(2, 3)*full["sum_Y_X"]
    assert bxx == F(32, 3) and byx == F(16, 3)

    gs, half = read(HERE / a["g3_source"])["g3"]
    gsl = gs - a["g3_range_half_widths"]*half
    hard_limit = (2/3)*(static_limit/gsl)**2
    transfer = a["comparison_transport_g_hard_over_g_static"]
    t0 = read(HERE / a["reference_t_source"])["U3_counterexample"]["relative_abelian_squared_coupling"]
    assert math.isclose(3*(8+t0*t0)/(8+t0)**2, .37, rel_tol=1e-12)
    g_example = a["example_g_nonzero"]
    # Subtracted exact formulae avoid losing 50-power shifts in a float near 1.
    dx = (9/4)*hard_limit*(16+hard_limit)/(8+hard_limit)**2
    dy = hard_limit*(48-21*hard_limit)/(8*(8+hard_limit)**2)
    return {
        "anomalies_per_generation_q_one_third_ell_minus_one":
            {name: {k: str(val) for k, val in values.items()} for name, values in
             [("with_right_nu", full), ("without_right_nu", minimal)]},
        "bXX_three_families": float(bxx), "bYX_three_families": float(byx),
        "N_over_Ar_Ti_alloy": nt, "N_over_Ar_Pt_alloy": np_,
        "difference_Pt_minus_Ti": np_-nt, "eta_over_g_squared_massless": response,
        "positive_eta95_reconstructed_g": reconstructed, "g95_rounded_used": static_limit,
        "g3_two_width_lower": gsl, "t_static_equal_hard_upper": hard_limit,
        "t_with_factor_two_upper": hard_limit*transfer**2,
        "ideal_x_drop_upper": dx, "ideal_y_drop_upper": dy,
        "relative_aligned_one_loop_shift_in_inverse_g2_upper":
            float(bxx)*static_limit**2*a["perturbative_aligned_log_interval_budget"]/(8*math.pi**2),
        "round1507_reference_t": t0,
        "minimum_g_hard_for_that_example": gsl*math.sqrt(1.5*t0),
        "minimum_hard_to_static_ratio_for_that_example": gsl*math.sqrt(1.5*t0)/static_limit,
        "allowed_nonzero_example": {"g": g_example, "eta": response*g_example**2,
                                    "t_at_central_g3": (2/3)*(g_example/gs)**2}
    }


def compare(actual, saved):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for key in actual:
            compare(actual[key], saved[key])
    elif isinstance(actual, float):
        assert math.isclose(actual, saved, rel_tol=3e-13, abs_tol=0), (actual, saved)
    else:
        assert actual == saved


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = certificate()
    path = HERE / "results.json"
    if args.write:
        with path.open("x", encoding="utf-8") as dst:
            json.dump(result, dst, indent=2, ensure_ascii=False, allow_nan=False)
            dst.write("\n")
    else:
        compare(result, read(path))
        print("Round 1508 checks passed.")
