"""Independent finite-polynomial bounds and a double-Higgs contact holdout."""

from argparse import ArgumentParser
from pathlib import Path
import json
from math import pi, sqrt
import numpy as np

HERE = Path(__file__).resolve().parent


def measure(coeff, x):
    n = np.arange(len(coeff))
    d = 2 * n + 1
    masses = np.asarray(coeff) * x ** n
    return sum(masses), sum(d * masses), sum(d * (d - 1) * masses)


def compute():
    # Re-read historical scientific inputs independently of the author's result.
    old = json.loads((HERE.parent / "1435" / "inputs.json").read_text(encoding="utf-8"))
    masses = json.loads((HERE.parent / "1386" / "inputs.json").read_text(encoding="utf-8"))["yukawa"]
    v = old["v_GeV"]
    output = {}
    for sector, channel in (("down", "bb_ZZ"), ("charged_leptons", "tautau_ZZ")):
        obs, _, upwards = old["atlas_table7"][channel]["observed"]
        ref, ref_err = old["atlas_table7"][channel]["SM"]
        k = sqrt((obs + 2 * upwards) / ((ref - 2 * ref_err) * old["rate_transfer_window"][0]))
        central, error = masses[sector][2]
        y = central - 2 * error
        caps_result = {}
        for cap in (1., 4 * pi):
            # Find root using finite tails of length 120 and interval bisection.
            ns = np.arange(2, 121, dtype=float)
            left, right = 0., 0.5
            threshold = (3 - k) * y / 2
            for _ in range(64):
                mid = (left + right) / 2
                if cap * sum((ns - 1) * mid ** ns) < threshold:
                    left = mid
                else:
                    right = mid
            x = right
            c = np.full(121, -cap)
            c[0] = 0.
            c[1] = (y - measure(c, x)[0] + c[1] * x) / x
            m, g, _ = measure(c, x)
            assert max(abs(m - y), abs(g - k * y)) < 1e-13 and max(abs(c)) <= cap
            caps_result[str(round(cap, 6))] = round(v / sqrt(2 * x), 8)
        test_x = 0.15
        coeff = [0., 2 * y / test_x, -y / test_x**2]
        m, g, curvature = measure(coeff, test_x)
        assert max(abs(m - y), abs(g - y), abs(curvature + 8*y)) < 1e-13 and max(map(abs,coeff)) < 1
        q = test_x/(1-test_x)
        # If |G| <= k y and c0=0, the signed double-h contact obeys this limit.
        t_upper = 7*k - 15 + 8*q**3/y
        assert t_upper < 0
        output[channel] = {
            "k_upper": round(k, 10), "no_direct_scale_upper_GeV": caps_result,
            "test_x": test_x, "test_scale_GeV": round(v/sqrt(2*test_x), 8),
            "single_h_SM_witness_c1_c2": [round(a, 10) for a in coeff[1:]],
            "witness_signed_T_over_y": -8,
            "T_over_y_upper_if_only_actual_single_h_bound": round(t_upper, 8),
            "necessary_x_for_SM_one_and_two_h": round(y**(1/3)/(1+y**(1/3)), 9),
        }
    # All three linear functionals can be tested on complex coefficients without fitting them.
    rng = np.random.default_rng(1438)
    for x in (0.11, 0.36, 0.59):
        c = rng.uniform(0, 1, 48) * np.exp(2j * pi * rng.uniform(0, 1, 48))
        y, g, t = measure(c, x)
        n = np.arange(3, len(c))
        remainder = sum((n-1)*(n-2)/2*c[3:]*x**n)
        assert abs((t-7*g+15*y)/8 - remainder-c[0]) < 1e-12
        assert abs(remainder) <= (x/(1-x))**3
    return output


if __name__ == "__main__":
    cli = ArgumentParser()
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = cli.parse_args()
    result = compute()
    target = HERE / "independent_results.json"
    if args.write:
        with target.open("x", encoding="utf-8") as file:
            json.dump(result, file, indent=2)
            file.write("\n")
    else:
        assert result == json.loads(target.read_text(encoding="utf-8"))
    print(json.dumps(result))
