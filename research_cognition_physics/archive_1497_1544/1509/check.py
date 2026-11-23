"""Finite-window subjet count predictions and sharp sample-bias thresholds."""

import argparse
import itertools
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def mean_box(pair, width, padding=0.0):
    center, error = pair
    return [center - width*error - padding, center + width*error + padding]


def mixture_scale(r, fraction_high, fraction_low):
    return (1 + fraction_high*(r-1))/(1 + fraction_low*(r-1))


def forward(low, r, fraction_high, fraction_low):
    return 1 + (low-1)*mixture_scale(r, fraction_high, fraction_low)


def invert(high, low, fraction_high, fraction_low):
    delta = (high-low)/(fraction_high-fraction_low)
    quark = low - fraction_low*delta
    assert quark > 1
    return {"gluon_minus_quark": delta, "gluon": quark+delta, "quark": quark,
            "extra_subjet_ratio": 1+delta/(quark-1)}


def bounds(values):
    values = list(values)
    return [min(values), max(values)]


def predictions(r, high, low, fh, fl):
    scales = [mixture_scale(r, f, g) for f, g in itertools.product(fh, fl)]
    predicted = bounds(1+(m-1)*s for m, s in itertools.product(low, scales))
    if predicted[0] > high[1]:
        # With independent bounded additive biases in both reads, their
        # maximal combined displacement is (1+s)*b for a fixed slope s.
        budget = (predicted[0]-high[1])/(1+min(scales))
    elif predicted[1] < high[0]:
        budget = (high[0]-predicted[1])/(1+max(scales))
    else:
        budget = 0
    return {"high_from_low": predicted, "minimum_per_sample_bias": budget}


def certificate():
    a = load(ROOT / "inputs.json")
    h, l = a["raw_H_M"], a["raw_L_M"]
    fh, fl = a["fraction_H"], a["fraction_L"]
    x = a["ideal_SU3"]
    center = invert(h[0], l[0], fh[1], fl[1])
    assert math.isclose(center["gluon_minus_quark"], .2/.26, rel_tol=1e-12)
    assert abs(center["extra_subjet_ratio"] - a["published_uncorrected_inversion"][0]) < .02
    cases = {}
    for w in a["working_widths"]:
        high, low = mean_box(h, w), mean_box(l, w)
        ranges = [high, low, [fh[0], fh[2]], [fl[0], fl[2]]]
        inv = [invert(*t) for t in itertools.product(*ranges)]
        raw_ratio = bounds(b["extra_subjet_ratio"] for b in inv)
        rec = a["corrected_species_ratio"]
        corrected = [rec["central"]-w*(rec["stat"]+rec["syst_minus"]),
                     rec["central"]+w*(rec["stat"]+rec["syst_plus"])]
        cases[str(w)] = {
            "gluon_minus_quark_detector_level": bounds(b["gluon_minus_quark"] for b in inv),
            "extra_subjet_ratio_detector_level": raw_ratio,
            "transfer_factor_to_raw_for_SU3": [v/x for v in raw_ratio],
            "forward_equal": predictions(1, *ranges),
            "forward_SU3": predictions(x, *ranges),
            "corrected_ratio_arithmetic_box": corrected,
            "transfer_factor_to_corrected_for_SU3": [v/x for v in corrected]
        }
    padding = a["rounding_padding_per_printed_mean_optional"]
    pad = predictions(1, mean_box(h, 2, padding), mean_box(l, 2, padding), fh, fl)
    # Exact equal-response witness at the two-width boundary.
    lo_high, hi_low = mean_box(h, 2)[0], mean_box(l, 2)[1]
    common = (lo_high+hi_low)/2
    b = cases["2"]["forward_equal"]["minimum_per_sample_bias"]
    assert math.isclose(common+b, lo_high) and math.isclose(common-b, hi_low)
    assert cases["2"]["forward_SU3"]["minimum_per_sample_bias"] == 0
    assert cases["1"]["forward_SU3"]["minimum_per_sample_bias"] > 0
    return {
        "central_raw_inversion": center,
        "central_high_from_low_if_equal": l[0],
        "central_high_from_low_if_casimir_transfers_unchanged": forward(l[0], x, fh[1], fl[1]),
        "two_width_equal_bias_witness": {"common": common, "bias_high": b, "bias_low": -b},
        "extra_printed_rounding_two_width_equal_bias": pad["minimum_per_sample_bias"],
        "two_width_SU3_witness": {"low": mean_box(l, 2)[0], "fh": fh[0], "fl": fl[2],
                                 "high": forward(mean_box(l, 2)[0], x, fh[0], fl[2])},
        "windows": cases
    }


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for name in a:
            compare(a[name], b[name])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for left, right in zip(a, b):
            compare(left, right)
    else:
        assert math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-14), (a, b)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = ROOT / "results.json"
    actual = certificate()
    if args.write:
        with path.open("x", encoding="utf-8") as f:
            json.dump(actual, f, ensure_ascii=False, indent=2, allow_nan=False)
            f.write("\n")
    else:
        compare(actual, load(path))
        print("Round 1509 checks passed.")
