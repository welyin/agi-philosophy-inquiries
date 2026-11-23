"""Distance from the measured up port to a self-adjoint return window."""
import argparse
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(rel):
    return json.loads((HERE / rel).read_text(encoding="utf-8-sig"))


def product(a, b):
    return sum(x * y for x, y in zip(a, b))


def closest(port, signed, tol):
    """Return the exact optimal angle and a checked witness, modulo bisection."""
    weights = [x * x for x in port]
    initial = product(weights, signed)
    if abs(initial) <= tol:
        return 0.0, port
    target = math.copysign(tol, initial)
    shifted = [x - target for x in signed]
    if not min(shifted) < 0 < max(shifted):
        # No interior level. The spectra in this test have no endpoint exactly zero.
        return None
    scale = max(map(abs, shifted))
    g = [x / scale for x in shifted]
    lo, hi = -1 / max(g), -1 / min(g)
    for _ in range(120):
        t = (lo + hi) / 2
        den = [1 + t * x for x in g]
        deriv_negative = sum(p * x / (y * y) for p, x, y in zip(weights, g, den))
        if deriv_negative > 0:
            lo = t
        else:
            hi = t
    t = (lo + hi) / 2
    w = [v / (1 + t * x) for v, x in zip(port, g)]
    norm = math.sqrt(product(w, w))
    w = [v / norm for v in w]
    overlap = product(port, w)
    # sin^2 = sum of all squared real 2x2 minors for positive optimal phases.
    sin2 = sum((port[i] * w[j] - port[j] * w[i]) ** 2 for i,j in ((0,1),(0,2),(1,2)))
    theta = math.degrees(math.atan2(math.sqrt(sin2), overlap))
    assert abs(product([x*x for x in w], signed) - target) < 2e-16
    assert abs(overlap**2 - sum(p/(1+t*x) for p,x in zip(weights,g))) < 4e-12
    return theta, w


def solve():
    rows = load("../1556/results.json")
    high_down = load("../1397/inputs.json")["high_scale"]["down"]
    first = high_down[0]
    d_window = [first[0]-2*sum(first[1:]), first[0]+2*sum(first[1:])]
    out = {"signed_light_magnitude_interval":d_window, "seeds":{}}
    for source in ("release_down_A", "release_down_B"):
        fit = dict(zip(rows["labels"],rows["candidates"][source]["values"]))
        d, s, b = [fit[k] for k in ("yd","ys","yb")]
        a12, a13 = [fit[k] for k in ("theta12_q","theta13_q")]
        port = [math.cos(a12)*math.cos(a13),math.sin(a12)*math.cos(a13),math.sin(a13)]
        selected = {}
        strong_angle_cases = {}
        signed_branches = []
        for sign_s, sign_b, post_sign in itertools.product((-1,1),repeat=3):
            initial = [d,sign_s*s,sign_b*b]
            cof = [initial[1]*initial[2],initial[0]*initial[2],initial[0]*initial[1]]
            # An unrestricted c could put the initial strange or bottom eigenvalue
            # into the measured d interval. Then BOTH other masses exceed the
            # strange upper edge, so no relabelling saves these fixed seeds.
            strange_upper = high_down[1][0] + 2*sum(high_down[1][1:])
            for other in (1,2):
                endpoints = [(post_sign*m-initial[other])/cof[other] for m in d_window]
                for k in set(range(3)) - {other}:
                    interval = [initial[k]+c*cof[k] for c in endpoints]
                    assert min(interval)>strange_upper or max(interval)<-strange_upper
            at = []
            for mv in d_window:
                c = (post_sign*mv-d)/cof[0]
                lam = [x+c*y for x,y in zip(initial,cof)]
                assert max(abs(lam[0]),abs(lam[1])) < abs(lam[2]) and abs(lam[0]) < abs(lam[1])
                at.append((mv,c,lam,product([x*x for x in port],lam)))
            assert at[0][3]*at[1][3] > 0
            signed_branches.append((sign_s,sign_b,post_sign,at))
        for tolerance in (0,1e-8,3e-8,5e-8,6e-8,1e-7):
            candidates = []
            for sign_s,sign_b,post_sign,at in signed_branches:
                for mv,c,lam,return_at_original in at:
                    near = closest(port,lam,tolerance)
                    if near is not None:
                        theta,w = near
                        candidates.append((theta,abs(return_at_original),
                            {"min_angle_deg":theta,"probability_leak":math.sin(math.radians(theta))**2,
                             "conditional_c":c,"sign_s_b_dpost":[sign_s,sign_b,post_sign],
                             "abs_signed_dpost":mv,"initial_return":return_at_original,
                             "port_abs":port,"attaining_port_abs":w,"signed_output":lam}))
            selected[str(tolerance)] = min(candidates,key=lambda q:(q[0],q[1]))[-1]
            if tolerance == 0:
                keys = {"same_det_sign":lambda r:r["sign_s_b_dpost"][2]==1,
                        "final_det_positive":lambda r:math.prod(r["signed_output"])>0,
                        "both_det_positive":lambda r:r["sign_s_b_dpost"][2]==1 and math.prod(r["signed_output"])>0}
                for key, include in keys.items():
                    strong_angle_cases[key] = min(x[0] for x in candidates if include(x[2]))
        out["seeds"][source] = {"minimum_mispointing_for_return_budget":selected,
            "minimum_angle_deg_if_strong_angle_signs_fixed":strong_angle_cases,
            "minimum_return_at_exact_port": min(abs(x[3]) for *_,at in signed_branches for x in at)}
    return out


def compact(x):
    if isinstance(x,dict):
        return {k:compact(v) for k,v in x.items()}
    if isinstance(x,list):
        return list(map(compact,x))
    if isinstance(x,float):
        return float(format(x, ".11g"))
    return x


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument("--check",action="store_true")
    g.add_argument("--write-new",action="store_true")
    args=p.parse_args()
    result=compact(solve())
    path=HERE/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,ensure_ascii=False,allow_nan=False)
            f.write("\n")
    elif result!=load("independent_results.json"):
        raise SystemExit("Independent angle certificate differs")
    for name, row in result["seeds"].items():
        print(name,row["minimum_return_at_exact_port"],
              {k:v["min_angle_deg"] for k,v in row["minimum_mispointing_for_return_budget"].items()})
