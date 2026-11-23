"""Two gauge-Higgs boundary conditions in the same MSbar convention."""
import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / "inputs.json").read_text(encoding="utf-8"))
LOOP = 16 * math.pi**2
SM = (41/6, -19/6, -7)
EXTRA = (3/2, 1/6, 0)
UNIT = DATA["mu0"]


def gauge(t, start, gate):
    return tuple(1 / (1/g**2 - 2*(b*t + d*max(0., t-gate))/LOOP)
                 for g, b, d in zip(start, SM, EXTRA))


def match_t(start, gate):
    gap = (1/start[0]**2 - 3/start[1]**2) * LOOP / 2
    slope = SM[0] - 3*SM[1]
    increment = EXTRA[0] - 3*EXTRA[1]
    t0 = gap/slope
    return t0 if t0 <= gate else (gap+increment*gate)/(slope+increment)


def top_beta(q, x, y, z):
    return q*(9*q - 17*x/6 - 9*y/2 - 16*z)/LOOP


def quartics(q, gauge_sq, state, active):
    x, y, _ = gauge_sq
    l, h, a, b = state
    gh = 3*q - 3*x/4 - 9*y/4
    gauge_h = 3*(3*y*y + x*x + 2*x*y)/8
    beta = 24*l*l + 4*gh*l + gauge_h - 6*q*q
    if not active:
        return (beta/LOOP, 0., 0., 0.)
    ge = -27*x/4 - 9*y/4
    gauge_e = 3*(3*y*y + 81*x*x + 18*x*y)/8
    portal = a*a + (a+b)**2
    beta_h = 24*h*h + 4*ge*h + gauge_e + portal
    beta_a = (4*(l+h)*(3*a+b) + 4*a*a + 2*b*b + 2*(gh+ge)*a
              + 3*(3*y*y + 9*x*x - 6*x*y)/4)
    beta_b = 4*(l+h)*b + 8*a*b + 4*b*b + 2*(gh+ge)*b + 9*x*y
    return tuple(v/LOOP for v in (beta+portal, beta_h, beta_a, beta_b))


def evolve(func, t0, t1, state, n=96):
    if t0 == t1:
        return state
    step = (t1-t0)/n

    def add(values, derivative, h):
        return tuple(v+h*d for v, d in zip(values, derivative))

    for i in range(n):
        t = t0+i*step
        k1 = func(t, state)
        k2 = func(t+step/2, add(state, k1, step/2))
        k3 = func(t+step/2, add(state, k2, step/2))
        k4 = func(t+step, add(state, k3, step))
        state = tuple(v + step*(a+2*b+2*c+d)/6
                      for v, a, b, c, d in zip(state, k1, k2, k3, k4))
    return state


def forward_prediction(divisor, mass=None, shift=0., steps=96):
    ref = DATA["central"]
    start = tuple(ref[key] for key in ("gY", "g2", "g3"))
    gate = math.log(mass/UNIT) if mass is not None else math.inf
    end = match_t(start, gate)
    middle = min(end, gate)
    y2 = gauge(end, start, gate)[1]

    def top(t, q):
        return (top_beta(q[0], *gauge(t, start, gate)),)

    q_gate = evolve(top, 0., middle, (ref["yt"]**2,), steps)
    q_end = evolve(top, middle, end, q_gate, steps)
    potential = (y2/divisor + shift, y2/2, 0., -y2/2) if divisor == 6 else (y2/divisor + shift, 0., 0., 0.)

    def flow(active):
        def rhs(t, values):
            g = gauge(t, start, gate)
            q, *p = values
            return (top_beta(q, *g), *quartics(q, g, p, active))
        return rhs

    at_gate = evolve(flow(True), end, middle, (*q_end, *potential), steps)
    low = evolve(flow(False), middle, 0., at_gate, steps)
    assert abs(low[0]-ref["yt"]**2) < 1e-9
    return {"match_GeV": UNIT*math.exp(end), "tree_lambda_at_match": y2/divisor,
            "predicted_lambda200": low[1], "low_residual": ref["lambda"]-low[1],
            "quartics_at_gate": list(at_gate[1:])}


def central_point(divisor, mass=None):
    out = forward_prediction(divisor, mass)
    ref = DATA["central"]["lambda"]
    shift = 0.
    for _ in range(5):
        value = forward_prediction(divisor, mass, shift)["predicted_lambda200"]
        slope = (forward_prediction(divisor, mass, shift+.00001)["predicted_lambda200"] - value)/.00001
        shift += (ref-value)/slope
    assert abs(forward_prediction(divisor, mass, shift, 160)["predicted_lambda200"] - ref) < 2e-10
    out["required_matching_delta_lambda"] = shift
    return out


def enclosure():
    """Continuous scalar comparison, including any one complete extra doublet gate."""
    box = DATA["box"]
    gauge_box = [box[k] for k in ("gY", "g2", "g3")]
    lo = tuple(pair[0] for pair in gauge_box)
    hi = tuple(pair[1] for pair in gauge_box)
    # Lower and upper matching times; the hypercharge and weak endpoints oppose.
    earliest = match_t((hi[0], lo[1]), 0.)
    latest = match_t((lo[0], hi[1]), math.inf)
    ql, qu = [v*v for v in box["yt"]]
    ll, lu = box["lambda"]
    lower, upper = math.inf, -math.inf
    n = 1800
    dt = latest/n
    speed = .17
    wiggle = speed*dt + 1e-11
    for k in range(n):
        a, c = k*dt, (k+1)*dt
        xl = gauge(a, lo, math.inf)[0]
        xu = gauge(c, hi, 0.)[0]
        yl, zl = gauge(c, lo, math.inf)[1:]
        yu, zu = gauge(a, hi, 0.)[1:]
        # A strict speed fence is verified over the entire prospective cell.
        q0, q1 = ql-wiggle, qu+wiggle
        l0, l1 = ll-wiggle, lu+wiggle
        assert .04 < l0 < l1 < q0 < q1 < 1.1
        gauge_min = 3*(3*yl**2 + xl**2 + 2*xl*yl)/8
        gauge_max = 3*(3*yu**2 + xu**2 + 2*xu*yu)/8

        def qmin(left, right):
            return (9*left**2-right*(17*xu/6+9*yu/2+16*zu))/LOOP

        def qmax(left, right):
            return (9*right**2-left*(17*xl/6+9*yl/2+16*zl))/LOOP

        def lmin(left, right):
            # For fixed l the term 12*q*l-6*q^2 decreases with q above l.
            return (24*left**2+12*q1*left-6*q1*q1-right*(9*yu+3*xu)+gauge_min)/LOOP

        def lmax(left, right):
            return (24*right**2+12*q0*right-6*q0*q0-left*(9*yl+3*xl)+gauge_max)/LOOP

        assert min(qmin(q0, q1), lmin(l0, l1)) > -speed
        assert max(qmax(q0, q1), lmax(l0, l1)) < speed
        if c >= earliest:
            lower = min(lower, l0-yu/6)
            upper = max(upper, l1-yl/2)
        # Endpoints of four separate comparison trajectories, with look-ahead.
        ql += dt*qmin(ql-wiggle, ql+wiggle) - 1e-13
        qu += dt*qmax(qu-wiggle, qu+wiggle) + 1e-13
        ll += dt*lmin(ll-wiggle, ll+wiggle) - 1e-13
        lu += dt*lmax(lu-wiggle, lu+wiggle) + 1e-13
    return {"match_GeV": [UNIT*math.exp(earliest), UNIT*math.exp(latest)],
            "G2_required_shift_at_least": lower,
            "Sp6_required_shift_at_most": upper}


def gauge_threshold_trade():
    ref = DATA["central"]
    start = tuple(ref[k] for k in ("gY", "g2", "g3"))

    def mismatch(t):
        def rhs(s, v):
            g = gauge(s, start, math.inf)
            q, l = v
            return (top_beta(q, *g), quartics(q, g, (l, 0., 0., 0.), False)[0])
        l = evolve(rhs, 0., t, (ref["yt"]**2, ref["lambda"]), 160)[1]
        return l - gauge(t, start, math.inf)[1]/6

    left, right = (math.log(m/UNIT) for m in (12000, 25000))
    assert mismatch(left) > 0 > mismatch(right)
    for _ in range(42):
        mid = (left+right)/2
        if mismatch(mid) > 0:
            left = mid
        else:
            right = mid
    t = (left+right)/2
    x, y, _ = gauge(t, start, math.inf)
    return {"G2_quartic_only_match_GeV": UNIT*math.exp(t),
            "sin2_thetaW": x/(x+y), "inverse_gauge_residual": 1/x-3/y}


def calculate():
    windows = enclosure()
    assert windows["G2_required_shift_at_least"] > 0.
    assert windows["Sp6_required_shift_at_most"] < 0.
    sm = central_point(6)
    return {"G2_SM_below_match": sm, "Sp6_SM_below_match": central_point(2),
            "G2_active_complex_doublet": {str(m): central_point(6,m) for m in DATA["hard_doublet_examples_GeV"]},
            "SM_below_match_weak_boundary_released": gauge_threshold_trade(),
            "continuous_comparison_bounds": windows}


def stable(obj):
    if isinstance(obj, dict):
        return {k: stable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [stable(v) for v in obj]
    return float(f"{obj:.9g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = stable(calculate())
    destination = ROOT / "results.json"
    if args.check:
        if json.loads(destination.read_text(encoding="utf-8")) != result:
            raise SystemExit("1621 result mismatch")
    else:
        with destination.open("x", encoding="utf-8") as f:
            f.write(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print("1621 common gauge scale, scalar matching and continuous bounds verified")
