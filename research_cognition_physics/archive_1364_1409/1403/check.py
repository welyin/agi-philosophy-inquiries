"""Gauge-selected scales and determinant ratios in the no-threshold one-loop SM."""
import argparse
import importlib.util
import json
from math import ceil, exp, expm1, log, pi
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
K = 8 * pi**2
B = np.array([41 / 10, -19 / 6, -7])


def read(name):
    return json.loads((ROOT / name).read_text(encoding="utf-8-sig"))


CFG = read("inputs.json")
spec = importlib.util.spec_from_file_location("r1397", ROOT / CFG["quark_and_CKM"])
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)


def gauge(g0, t):
    return 1 / (1 / g0**2 - B * t / K)


def gauge_root(g2, g3):
    return K * (1 / g2**2 - 1 / g3**2) / (B[1] - B[2])


def log_gauge_integral(g1, g3, t):
    def integrated(g, b):
        return -K / b * np.log1p(-b * g**2 * t / K)
    return (6 * integrated(g1, B[0]) - 24 * integrated(g3, B[2])) / (2 * K)


def bracketing_root(f, lo, hi):
    assert f(lo) * f(hi) <= 0
    for _ in range(55):
        mid = (lo + hi) / 2
        if f(mid) * f(lo) <= 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def central(g0, us, ds, es, V, t, steps):
    state = np.array([np.diag(us), V * ds, np.diag(es)], dtype=complex)
    identity = np.eye(3)

    def derivative(s, y):
        g1, g2, g3 = gauge(g0, s)
        Hu, Hd, He = [m @ m.conj().T for m in y]
        T = np.trace(3 * (Hu + Hd) + He).real
        Su = 17 * g1 / 20 + 9 * g2 / 4 + 8 * g3
        Sd = g1 / 4 + 9 * g2 / 4 + 8 * g3
        Se = 9 * (g1 + g2) / 4
        return np.array([((T - Su) * identity + 1.5 * (Hu - Hd)) @ y[0],
                         ((T - Sd) * identity + 1.5 * (Hd - Hu)) @ y[1],
                         ((T - Se) * identity + 1.5 * He) @ y[2]]) / (2 * K)

    step = t / steps
    for i in range(steps):
        s = i * step
        a = derivative(s, state)
        b = derivative(s + step / 2, state + step * a / 2)
        c = derivative(s + step / 2, state + step * b / 2)
        d = derivative(s + step, state + step * c)
        state += step * (a + 2*b + 2*c + d) / 6
    return abs(np.linalg.det(state[1]) / np.linalg.det(state[2]))


def endpoints(data, k):
    central_values, widths = np.array(data).T
    return central_values - k * widths, central_values + k * widths


def certificate(gauges, us, ds, es, k, T, ref):
    gl, gh = endpoints(gauges, k)
    ul, uh = endpoints(us, k)
    dl, dh = endpoints(ds, k)
    el, eh = endpoints(es, k)
    tlo = gauge_root(gh[1], gl[2])
    thi = gauge_root(gl[1], gh[2])
    assert 0 < tlo < thi < T
    cap, ecap = CFG["traces"]["Qcap"], CFG["traces"]["Ecap"]
    q = float(uh @ uh + dh @ dh)
    E0 = eh @ eh
    Emax = E0 * exp((6 * cap + 5 * ecap) * T / (2 * K))
    assert Emax < ecap and q < cap
    # Stop a hypothetical first escape at either cap. The lepton cap cannot be
    # reached, and the analytic Riccati comparison checks the quark cap on each slab.
    ns = ceil(T * 4)
    h = T / ns
    max_q = q
    for i in range(ns):
        _, weak, color = gauge(gl, (i + 1) * h)
        damp = 16 * color + 4.5 * weak - 2 * ecap
        arg = damp * h / (2 * K)
        inv_next = 1 / q + expm1(arg) * (1 / q - 9 / damp)
        assert inv_next > 0
        q = 1 / inv_next
        max_q = max(q, max_q)
    assert max_q < cap
    # Q+E bounds the absolute trace contribution. Both envelopes strictly decrease.
    tmax_g = gauge(gh, T)
    tmin_g = gauge(gl, T)
    biggest_slope = 6 * tmax_g[0] - 24 * tmin_g[2] + 1.5 * (cap + ecap)
    assert biggest_slope < 0
    rlo, rhi = np.prod(dl) / np.prod(eh), np.prod(dh) / np.prod(el)
    def lower(t):
        return log(rlo) + log_gauge_integral(gl[0], gh[2], t) - 1.5*(cap+ecap)*t/(2*K)
    def upper(t):
        return log(rhi) + log_gauge_integral(gh[0], gl[2], t) + 1.5*(cap+ecap)*t/(2*K)
    vol = [exp(lower(thi)), exp(upper(tlo))]
    gains = [r**(1/3) for r in vol]
    ans = {"log_gauge_crossing": [tlo, thi], "gauge_crossing_GeV": [ref*exp(tlo), ref*exp(thi)],
           "R_at_crossing": vol, "required_kd_over_ke": gains, "E_cap_certificate": Emax,
           "max_piecewise_Q_bound": max_q, "max_16pi2_logR_slope": biggest_slope,
           "min_additive_logR_shift": -log(vol[1])}
    trial_shift = CFG["illustrative_inverse_gauge_threshold"]
    early = tlo - K * trial_shift / (B[1] - B[2])
    late = thi + K * trial_shift / (B[1] - B[2])
    assert 0 < early < late < T
    ans["min_log_volume_shift_with_gauge_budget"] = -upper(early)
    if rlo > 1:
        # These scalar envelopes bracket every actual R=1 scale.
        t1 = bracketing_root(lower, 0, T)
        t2 = bracketing_root(upper, 0, T)
        assert t2 < tlo
        ans["all_equal_volume_crossings_GeV"] = [ref*exp(t1),ref*exp(t2)]
        ans["min_gauge_inverse_square_mismatch_at_equal_volume"] = (tlo-t2)*(B[1]-B[2])/K
    return ans


def calculate():
    lowg = read(CFG["MZ"])
    all_g = [[lowg["g1_GUT"], lowg["g2"], lowg["g3"]], CFG["high_gauge"]]
    leps = [read(CFG["MZ_charged_leptons"])["yukawa"]["charged_leptons"],
            read(CFG["high_charged_leptons"])["scales"]["1e16"][3:]]
    output = {}
    for index, (label, us, ds, angles, g2, radian) in enumerate(prior.data()):
        ref, T = CFG["reference_GeV"][index], CFG["end_log"][index]
        gs, es = all_g[index], leps[index]
        assert gs[1] == g2
        g0, u0, d0, e0, a0 = (np.array(v)[:, 0] for v in (gs, us, ds, es, angles))
        V = prior.mix(np.sin(a0[:3]) if radian else a0[:3], a0[3])
        t = gauge_root(g0[1], g0[2])
        R = central(g0, u0, d0, e0, V, t, 256)
        fine = central(g0, u0, d0, e0, V, t, 512)
        assert abs(fine-R) < 2e-9
        gg = np.sqrt(gauge(g0, t))
        assert abs(gg[1]-gg[2]) < 1e-12 and gg[0] > gg[1]
        boxes = {str(k): certificate(gs,us,ds,es,k,T,ref) for k in CFG["widths"]}
        assert all(b["R_at_crossing"][0] < fine < b["R_at_crossing"][1] < 1 for b in boxes.values())
        output[label] = {"log_crossing": t, "crossing_GeV": ref*exp(t),
                         "g1_g2_g3_at_crossing": gg.tolist(), "R": fine,
                         "kd_over_ke": fine**(1/3), "RK_refinement": abs(fine-R), "boxes": boxes}
    return output


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    group = args.add_mutually_exclusive_group(required=True)
    group.add_argument("--check", action="store_true")
    group.add_argument("--write-new", action="store_true")
    opts = args.parse_args()
    result = calculate()
    path = ROOT / "results.json"
    if opts.write_new:
        with path.open("x", encoding="utf-8") as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        old = read("results.json")
        def agrees(a,b):
            if isinstance(a,dict):
                return a.keys()==b.keys() and all(agrees(a[k],b[k]) for k in a)
            return np.allclose(a,b,rtol=5e-9,atol=1e-12)
        assert agrees(result,old)
        print("1403 scale, full-window determinant flow and two matching prescriptions verified")
