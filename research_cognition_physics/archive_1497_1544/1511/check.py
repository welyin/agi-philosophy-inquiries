"""Compact axion anomaly tests for the flavor texture and a competing matter realization."""
import argparse
import json
from fractions import Fraction as F
from math import pi
from pathlib import Path

HERE = Path(__file__).resolve().parent


def anomaly(field):
    color = field["color_dimension"]
    weak = field["weak_dimension"]
    assert color in (1, 3) and weak in (1, 2)
    weight = field.get("multiplicity", 1) * (field.get("PQ") or field.get("PQ_difference", 0))
    return (weight * weak * (color == 3), weight * color * (weak == 2),
            weight * color * weak * field["y"]**2)


def descends(field, q):
    triality = int(field["color_dimension"] == 3)
    parity = int(field["weak_dimension"] == 2)
    return (2 * triality + 3 * parity + field["y"]) % q == 0


def fits(coeffs, q):
    k, w, b = coeffs
    return b % q == 0 and (24*k + 18*w + b) % (q*q) == 0


def calculate():
    inp = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    groups = inp["groups"]
    old = json.loads((HERE/inp["prior_QCD"]).read_text(encoding="utf-8"))
    original = anomaly(inp["new_PQ_charged_Weyl_fields_of_old_texture"][0])
    fields = inp["vectorlike_competitor"]
    other = tuple(map(sum, zip(*(anomaly(f) for f in fields))))
    assert original == (-2,-3,-6) and other == (1,0,30)
    assert all(fits(original, q) for q in groups)
    assert [q for q in groups if fits(other,q)] == [1,3]
    assert [q for q in groups if all(descends(f,q) for f in fields)] == [1,3]
    for k in range(-12,13):
        for w in range(-12,13):
            for q in groups:
                assert fits((k,w,30*k-18*w),q) == (54*k % (q*q) == 0)
    # For |k|=1 the unsigned slope is minimized by the nearest point of
    # the q=6 lattice to the interval of calibrated QCD offsets.
    c, dc = old["photon_lattice"]
    widths = old["window_widths"]
    cl, ch = c - widths*dc, c + widths*dc
    z, stat, syst = old["chisus_4throot_MeV"]
    zl, zh = (z-widths*(stat+syst))/1000, (z+widths*(stat+syst))/1000
    norm = 1/(old["alpha_inverse"]*2*pi)
    ratio = F(5,3)
    candidate = [norm*(cl-float(ratio))/zh**2, norm*(ch-float(ratio))/zl**2]
    left, right = F(2,3), F(8,3)
    assert float(left) <= cl <= ch <= float(right)
    q6_floor = norm*min(cl-float(left),float(right)-ch)/zh**2
    mirrored = F(20,9)
    assert fits((1,0,40),2)
    alias_offset = (ratio + mirrored)/2
    assert cl < float(alias_offset) < ch and abs(ratio-alias_offset) == abs(mirrored-alias_offset)
    assert candidate[1] < q6_floor
    return {
        "texture_K3_K2_K1": original, "competing_K3_K2_K1": other,
        "q_compatible_at_given_ratio": {"odd_k":[1,3],"even_k":[1,2,3,6]},
        "all_q_forbid_odd_k_at_weak_to_color_ratio":"3/2",
        "common_offset_for_q2_unsigned_alias":str(alias_offset),
        "q2_alias_E_over_N":str(mirrored),
        "unsigned_slope_GeV_minus2": {"target_5_over_3":candidate,"q6_unit_winding_floor":q6_floor,
                                     "separation":q6_floor-candidate[1]},
        "signed_bare_ratio_distance_to_nearest_q2_unit_winding":str(F(1,9)),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    options = parser.parse_args()
    # Normalize tuples as JSON arrays for the read-only reproduction.
    result = json.loads(json.dumps(calculate()))
    output = HERE/"results.json"
    if options.write_new:
        with output.open("x", encoding="utf-8") as file:
            json.dump(result, file, indent=2)
            file.write("\n")
    else:
        assert result == json.loads(output.read_text(encoding="utf-8"))
    print("Verified anomaly coefficients, quotient tests and calibrated slope gap.")
