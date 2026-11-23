"""SM proxy inverse boundary, measured top holdout and phase weights."""
import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
LOOP = 16 * math.pi**2
BETA = (41 / 6, -19 / 6, -7)


def gauge(t, low):
    return tuple(g*g / (1 - 2*b*g*g*t/LOOP) for g, b in zip(low, BETA))


def critical_y_squared(x, y):
    return math.sqrt(2*y*y + (x+y)**2) / 4


def transport(t, state, low):
    x, y, z = gauge(t, low)
    q, lam = state
    dq = q * (9*q - 17*x/6 - 9*y/2 - 16*z)
    proxy = critical_y_squared(x, y)
    dl = 24*lam**2 + (12*q - 9*y - 3*x)*lam + 6*(proxy**2 - q**2)
    return dq / LOOP, dl / LOOP


def propagate(start, end, state, low, count):
    step = (end - start) / count
    q, lam = state
    for idx in range(count):
        t = start + idx*step
        a, b = transport(t, (q, lam), low)
        c, d = transport(t+step/2, (q+step*a/2, lam+step*b/2), low)
        e, f = transport(t+step/2, (q+step*c/2, lam+step*d/2), low)
        g, h = transport(t+step, (q+step*e, lam+step*f), low)
        q += step*(a+2*c+2*e+g)/6
        lam += step*(b+2*d+2*f+h)/6
    return q, lam


def terminal(t, low, count=256):
    x, y, _ = gauge(t, low)
    q, lam = propagate(t, 0, (critical_y_squared(x, y), 0), low, count)
    return math.sqrt(q), lam


def solve(scale, low, target, bounds, count=256):
    left, right = (math.log(m/scale) for m in bounds)
    assert terminal(left, low, count)[1] < target < terminal(right, low, count)[1]
    for _ in range(38):
        middle = (left + right)/2
        if terminal(middle, low, count)[1] < target:
            left = middle
        else:
            right = middle
    t = (left+right)/2
    yt = terminal(t, low, count)[0]
    x, y, z = gauge(t, low)
    dx, dy = (2*b*v*v/LOOP for b, v in zip(BETA, (x, y)))
    qcrit = critical_y_squared(x, y)
    dG = 3/8 * (4*y*dy + 2*(x+y)*(dx+dy))
    dqt = transport(t, (qcrit, 0), low)[0]
    beta_prime = (dG-12*qcrit*dqt)/LOOP
    assert beta_prime > 0
    return {"t": t, "field_coordinate": scale*math.exp(t), "low_y": yt,
            "high_q": qcrit, "beta_lambda_prime": beta_prime}


def analyze(cfg):
    ref, cal, wid = cfg["reference"], cfg["matching"], cfg["marginal_width"]
    low = (ref["gY"], ref["g2"], ref["g3"])
    unit, bounds = cfg["matching_scale"], cfg["window"]
    central = solve(unit, low, ref["lambda"], bounds, 512)
    coarse = solve(unit, low, ref["lambda"], bounds)
    assert abs(coarse["t"] - central["t"]) < 1e-6
    assert abs(coarse["low_y"] - central["low_y"]) < 1e-8
    two = 2*wid["alpha_s"]*cal["g3_per_alpha"]
    small = (low[0], low[1], low[2]*(1-two))
    large = (low[0], low[1], low[2]*(1+two))
    quartic_step = 2*ref["lambda"]*wid["mh"]*cal["lambda_per_relative_GeV"]
    min_t = solve(unit, large, ref["lambda"]-quartic_step, bounds)
    max_t = solve(unit, small, ref["lambda"]+quartic_step, bounds)
    y_floor = terminal(min_t["t"], small)[0]
    y_ceiling = terminal(max_t["t"], large)[0]
    x0, y0, _ = gauge(0, low)
    xmax, _, zmin = gauge(math.log(bounds[-1]/unit), small)
    h_floor = 16*zmin + 17*x0/6 - 11*y0/6 - 9*critical_y_squared(xmax, y0)
    barrier = 16*zmin - 9*y0/2 - xmax/6
    assert h_floor > 0 and barrier > 0
    high = cfg["masina_quiros"]["critical_top"]
    def y_of_mass(mass):
        return ref["top_yukawa"]*(1+cal["top_per_relative_GeV"]*(mass-ref["top_mass"]))
    compare = {}
    for pdf, (mass, err) in cfg["cms_normalized_ttjet"].items():
        diff = mass - high
        compare[pdf] = {"difference_GeV": diff, "difference_in_quoted_widths": diff/err,
                        "extra_halfwidth_for_one_width_overlap": max(0, abs(diff)-err),
                        "two_width_mass_range": [mass-2*err, mass+2*err],
                        "top_yukawa_central": y_of_mass(mass),
                        "top_yukawa_one_width": [y_of_mass(mass-err), y_of_mass(mass+err)]}
    old = cfg["iacobellis_masina_2016"]
    old_budget = sum(old["errors"])
    assert all(m-e < old["mass"]+old_budget and m+e > old["mass"]-old_budget
               for m,e in cfg["cms_normalized_ttjet"].values())
    tn = central["t"]
    step = .01
    left, right = terminal(tn-step, low), terminal(tn+step, low)
    restored = propagate(0, tn, (central["low_y"]**2, ref["lambda"]), low, 1024)
    assert abs(restored[0]-central["high_q"]) < 1e-8 and abs(restored[1]) < 1e-8
    history = json.loads((ROOT/"../1493/independent_results.json").read_text(encoding="utf-8-sig"))
    hist = solve(history["boundary_scale_GeV"],history["gauge_ordinary_Y_2_3"],
                 history["input_lambda"],bounds)
    # At a true stationary point its displacement makes no first-order
    # contribution to the change of vacuum energy.
    phi = cfg["masina_quiros"]["critical_field"]
    height = cfg["masina_quiros"]["barrier_fourth_root"]
    cutoff = bounds[1]
    shift_per_c6 = (phi/cutoff)**2 / 2
    c6_barrier = 4*(height/phi)**4 / shift_per_c6
    cutoff_t = math.log(1e16 / unit)
    x, y, _ = gauge(cutoff_t, low)
    q_endpoint = 1.1*critical_y_squared(x, y)
    unequal = propagate(cutoff_t, 0, (q_endpoint, 0), low, 512)
    finite_boundary={"endpoint_GeV":1e16, "low_lambda":unequal[1], "low_y":math.sqrt(unequal[0]),
                     "endpoint_beta":transport(cutoff_t, (q_endpoint, 0),low)[1]}
    assert finite_boundary["endpoint_beta"] < 0 and unequal[1] > 0
    return {"leading_central":central, "older_matching_173p34":hist,
            "endpoint_zero_without_stationarity":finite_boundary,
            "dimension_six_diagnostic":{"cutoff":cutoff, "d_4V_over_phi4_dc":shift_per_c6,
                                        "first_order_barrier_sized_c6":c6_barrier},
            "terminal_rows":{str(m):terminal(math.log(m/unit),low) for m in cfg["terminal_mass_scales"]},
            "double_width_proxy":{"field_interval":[min_t["field_coordinate"],max_t["field_coordinate"]],
                                  "y_outer_interval":[y_floor,y_ceiling],
                                  "lambda_change":quartic_step,
                                  "minimum_scaled_crossing_slope":h_floor,
                                  "minimum_barrier_coefficient":barrier},
            "local_shooting":{"d_lambda_ref_per_log_scale":(right[1]-left[1])/(2*step),
                              "d_y_ref_per_lambda_ref":(right[0]-left[0])/(right[1]-left[1])},
            "high_order":{"conditional_yt_200": y_of_mass(high),"CMS":compare,
                          "historical_2016_separate_error_rectangle":[old["mass"]-old_budget,old["mass"]+old_budget]}}


def tidy(x):
    if isinstance(x, dict):
        return {k:tidy(v) for k,v in x.items()}
    if isinstance(x, (list, tuple)):
        return [tidy(v) for v in x]
    return float(f"{x:.9g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    cfg = json.loads((ROOT/"inputs.json").read_text(encoding="utf-8-sig"))
    result = tidy(analyze(cfg))
    output = ROOT/"results.json"
    if args.check:
        assert result == json.loads(output.read_text(encoding="utf-8-sig"))
        print("1608: inverse double-zero, nuisance enclosure, matching and two PDF readings passed.")
    else:
        with output.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write("\n")
