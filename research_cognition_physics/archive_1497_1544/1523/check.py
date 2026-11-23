"""One-loop top-singular-value comparison on the first shared MSSM window."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOOP = 8 * math.pi**2
B = np.array([33/5, 1., -3.])
C = np.array([13/15, 3., 16/3])


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf8"))


def matching_gauge(source, table, widths, low, sign):
    names = ("g1_GUT", "g2", "g3")
    start = np.array([source[key][0] + sign * widths * source[key][1] for key in names])
    sm_b = np.array([4.1, -19/6, -7])
    casimir = np.array([0., 2., 3.])
    old = 1 / np.sqrt(start**-2 - sm_b * math.log(low/source["mu_ref_GeV"])/LOOP
                      - casimir/(48*math.pi**2))
    published = [table[key][0] + sign*widths*table[key][1] for key in ("g1", "g2", "g3")]
    return old, np.array(published)


def gauge_factor(time, gauge):
    denominators = 1 - np.asarray(time)[..., None] * B * gauge**2 / LOOP
    assert np.all(denominators > 0)
    return np.exp(np.log(denominators) @ (C / B))


def integral(time, gauge, count):
    # E decreases strictly. Endpoint Riemann sums bracket its integral.
    x = np.linspace(0., time, count+1)
    values = gauge_factor(x, gauge)
    centre = (time/count) * (float(np.sum(values[1:-1])) + (1+values[-1])/2)
    error = (time/count) * (1-values[-1])/2
    return centre, error


def calculate_window(gauge, yt, low, high, sine):
    end = math.log(high/low)
    factor = gauge_factor(end, gauge)
    total, err = integral(end, gauge, 32768)
    def cap(at_end, integ):
        return 1/math.sqrt(6*integ/LOOP + (0 if math.isinf(at_end) else factor/at_end**2))

    h0 = yt/sine
    out = {}
    for name, upper in (("one_loop_pole", math.inf), ("yt_squared_4pi", math.sqrt(4*math.pi))):
        limits = sorted(cap(upper, total + q*err) for q in (-1, 1))
        cut = cap(upper, total)
        min_sine = yt/cut
        def residual(t):
            s, _ = integral(t, gauge, 512)
            e = gauge_factor(t, gauge)
            return 1/h0**2 - 6*s/LOOP - (0 if math.isinf(upper) else e/upper**2)
        if residual(end) < 0:
            lo, hi = 0., end
            for _ in range(34):
                mid = (lo+hi)/2
                if residual(mid) > 0:
                    lo = mid
                else:
                    hi = mid
            tc = (lo+hi)/2
            # Verify a conservative time bracket against monotone quadrature.
            times = [tc-.002, tc+.002]
            for pos, t in enumerate(times):
                s, er = integral(t, gauge, 32768)
                r = 1/h0**2-6*s/LOOP-(0 if math.isinf(upper) else gauge_factor(t,gauge)/upper**2)
                assert (r > 6*er/LOOP) if pos == 0 else (r < -6*er/LOOP)
            hit = [low*math.exp(t) for t in times]
        else:
            hit = None
        out[name] = dict(maximum_initial_top_for_endpoint=cut, quadrature_bracket=limits,
                         necessary_tan_beta_at_zero_threshold=min_sine/math.sqrt(1-min_sine**2),
                         necessary_delta_at_tan_beta_one=yt/(sine*cut)-1,
                         comparison_crossing_GeV_bracket=hit)
    return out


def compute():
    cfg = read("inputs.json")
    old, table = read(cfg["gauge_base"]), cfg["sm_dr_at_one_TeV"]
    low = cfg["low_scale_GeV"]
    high = min(m for m in read(cfg["thresholds"])["witness"]["even_doublets_GeV"] if m > low)
    sine = math.sin(math.atan(cfg["tan_beta_at_low_scale"]))
    result = {}
    for label, gauge_sign, yt_sign in (("central", 0, 0), ("slow_growth", 1, -1), ("fast_growth", -1, 1)):
        n = cfg["experimental_half_widths"]
        gauges = matching_gauge(old, table, n, low, gauge_sign)
        yt = table["yt"][0] + yt_sign*n*table["yt"][1]
        for tag, gauge in zip(("1518", "table4"), gauges):
            value = calculate_window(gauge, yt, low, high, sine)
            assert 0 < value["one_loop_pole"]["necessary_delta_at_tan_beta_one"]
            assert value["yt_squared_4pi"]["comparison_crossing_GeV_bracket"][1] < high
            result[tag+"_"+label] = {"gauge_at_low":gauge.tolist(), "matched_top":yt/sine, **value}
    return {"first_heavy_doublet_GeV":high, "windows":result}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for key in a:
            compare(a[key],b[key])
    elif a is None:
        assert b is None
    else:
        np.testing.assert_allclose(a,b,rtol=1e-8,atol=1e-12)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--check",action="store_true")
    flags.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    data=compute()
    if args.write_new:
        with (HERE / "results.json").open("x", encoding="utf8") as out:
            json.dump(data,out,indent=2)
            out.write("\n")
    else:
        compare(data,read("results.json"))
    print("Verified top transport bounds at the first physical heavy threshold.")
