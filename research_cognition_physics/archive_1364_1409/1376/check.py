"""Count physical crossings and compare light spectra in the same threshold EFT."""
import argparse
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
K = 8 * math.pi**2
B = np.array([-41/6, 19/6, 7.])
SHIFT = np.array([-28/3, -28/3, 0.])
TRACE = np.array([1., -3., 4/3])


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def load_flavour_source(name):
    spec = importlib.util.spec_from_file_location("old_common_flavour", HERE / name)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    return old


def inverse_input(g):
    xs = 1/np.asarray(g)**2
    xs[0] *= 5/3
    return xs


def running(xs, t, ell):
    return xs + (B*t + SHIFT*max(0., t-ell))/K


def old_root(xs):
    return K*(TRACE @ xs)/7


def weight(xs):
    return 3*xs[2]/(4*xs[1]-3*xs[2])


def roots(xs, ell, end):
    first = old_root(xs)
    # The cusp is one distinct root when threshold and old crossing coincide.
    if ell < first:
        return []
    second = (8*ell-3*first)/5
    candidates = [first] if ell == first or second > end else [first, second]
    return [t for t in candidates if 0 <= t <= end]


def high_branch(xs, t):
    ell = (5*t + 3*old_root(xs))/8
    return ell, weight(running(xs, t, ell))


def log_mass_for_weight(xs, target):
    first = old_root(xs)
    t = (4*(K*xs[1] + 3.5*first) - (3 + 3/target)*K*xs[2])
    t /= 4/3 + 7*(3 + 3/target)
    ell, check = high_branch(xs, t)
    assert math.isclose(check, target, rel_tol=2e-12)
    return ell


def flavour_bounds(old, gauge_lo, end):
    upper, lower = old.Y + old.Y_ERR, old.Y - old.Y_ERR
    u0 = np.sum(upper[:2]**2)
    e0 = np.sum(upper[2]**2)
    beta = 7*gauge_lo[2]**2/K
    z = (1 + beta*end)**(-1/7)
    den = 1 - 9*u0*7*(1-z)/(2*K*beta)
    ecap = 0.001
    emin_g2 = 1/(gauge_lo[1]**-2+B[1]*end/K)
    assert 0 < z**8 < den <= 1 and u0 < 1 and 2*ecap < 4.5*emin_g2
    emax = e0*math.exp((6+5*ecap)*end/(2*K))
    assert emax < ecap
    fac = den**(-1/6)
    rd = [lower[1, 1]/upper[1, 0]/fac, upper[1, 1]/lower[1, 0]*fac]
    re = [lower[2, 1]/upper[2, 0],
          upper[2, 1]/lower[2, 0]*math.exp(1.5*emax*end/(2*K))]
    w = [old.minimum_flavour_weight(rd[1], re[0]), old.minimum_flavour_weight(rd[0], re[1])]
    return {"comparison_D": den, "z_to_eighth": z**8, "E_upper": emax,
            "rd_box": rd, "re_box": re, "minimum_weight_box": w}


def initial_matrices(old):
    ckm = old.rotation(1, 2, .04193) @ old.rotation(0, 2, .00370, 1.139) @ old.rotation(0, 1, .2251)
    return np.array([np.diag(old.Y[0]), ckm @ np.diag(old.Y[1]), np.diag(old.Y[2])])


def transport(xs, ell, finish, steps, old):
    def rhs(time, mats):
        y, w, c = 1/running(xs, time, ell)
        up, down, lep = mats @ mats.conj().transpose(0, 2, 1)
        trace = np.trace(3*(up+down)+lep).real
        gauge = np.array([17*y/12+9*w/4+8*c, 5*y/12+9*w/4+8*c, 15*y/4+9*w/4])
        matrix = (trace-gauge)[:, None, None]*np.eye(3) + 1.5*np.array([up-down, down-up, lep])
        return (matrix @ mats)/(2*K)

    mats = initial_matrices(old)
    intervals = [0, min(ell, finish), finish]
    for start, stop in zip(intervals, intervals[1:]):
        if stop <= start:
            continue
        h = (stop-start)/steps
        for index in range(steps):
            t = start + h*index
            a = rhs(t, mats)
            b = rhs(t+h/2, mats+h*a/2)
            c = rhs(t+h/2, mats+h*b/2)
            d = rhs(t+h, mats+h*c)
            mats += h*(a+2*b+2*c+d)/6
    dval, evals = np.sort(np.linalg.svd(mats[1:], compute_uv=False), axis=1)
    rd, re = dval[1]/dval[0], evals[1]/evals[0]
    return np.array([rd, re, old.minimum_flavour_weight(rd, re)])


def compute():
    cfg = read("inputs.json")
    data = read(cfg["gauge_and_yukawa"])
    old = load_flavour_source(cfg["flavour_initial_values_code"])
    mz = data["mu_ref_GeV"]
    end = math.log(cfg["upper_scale_GeV"]/mz)
    rows = np.array([data[k] for k in ("g1_GUT", "g2", "g3")])
    center, err = rows[:, 0], rows[:, 1]
    xs = inverse_input(center)
    corners = [inverse_input(center+err*s) for s in itertools.product((-1, 1), repeat=3)]
    envelope = flavour_bounds(old, center-err, end)
    wl, wh = envelope["minimum_weight_box"]
    energy = lambda t: mz*math.exp(t)
    thresholds = {(label): [energy(fun(v)) for v in corners] for label, fun in [
        ("cusp", old_root), ("high_end", lambda v: (5*end+3*old_root(v))/8),
        ("necessary", lambda v: log_mass_for_weight(v, wl)),
        ("sufficient", lambda v: log_mass_for_weight(v, wh))]}
    cutoff = {label: [min(vals), max(vals)] for label, vals in thresholds.items()}
    assert all(min(running(v, end, (5*end+3*old_root(v))/8)) > 0
               and high_branch(v, end)[1] > wh for v in corners)
    tests = []
    above_d = []
    for mass in cfg["threshold_benchmarks_GeV"]:
        ell = math.log(mass/mz)
        times = roots(xs, ell, end)
        sample = {"threshold_GeV": mass, "roots": []}
        assert all(len(roots(v, ell, end)) == len(times) for v in corners)
        for j, t in enumerate(times):
            matched = running(xs, t, ell)
            assert abs(TRACE @ matched) < 2e-14 and min(matched)>0
            ratios = transport(xs, ell, t, 128, old)
            assert np.allclose(ratios, transport(xs, ell, t, 64, old), rtol=2e-9)
            w = weight(matched)
            ww = [weight(running(v, roots(v, ell, end)[j], ell)) for v in corners]
            row = {"mu_GeV": energy(t), "w": w, "w_box": [min(ww), max(ww)],
                   "central_rd_re_wmin": ratios.tolist(),
                   "full_box_verdict": "pass" if min(ww)>=wh else "fail" if max(ww)<wl else "unresolved"}
            sample["roots"].append(row)
            if j == 1:
                for logfactor in (-math.log(2), math.log(2)):
                    dec = ell+logfactor
                    matched_elsewhere = xs+B*dec/K+SHIFT*(dec-ell)/K+(B+SHIFT)*(t-dec)/K
                    above_d.append(max(abs(matched_elsewhere-matched)))
        tests.append(sample)
    assert max(above_d)<4e-15
    assert [len(x["roots"]) for x in tests] == [0, 2, 2, 1]
    assert tests[2]["roots"][1]["full_box_verdict"]=="pass"
    assert tests[1]["roots"][1]["full_box_verdict"]=="fail"
    # One central-data solution of equality; uniqueness is not assumed.
    left, right = math.log(2e12/mz), math.log(1e13/mz)
    def excess(ell):
        t = (8*ell-3*old_root(xs))/5
        return weight(running(xs,t,ell))-transport(xs,ell,t,48,old)[2]
    assert excess(left)<0<excess(right)
    for _ in range(35):
        mid=(left+right)/2
        if excess(mid)>0:
            right=mid
        else:
            left=mid
    ell=(left+right)/2
    t=(8*ell-3*old_root(xs))/5
    return {"global_flavour_box": envelope, "threshold_landmarks_GeV":cutoff,
            "benchmarks": tests, "central_same_scale_equality":{
                "threshold_GeV":energy(ell),"mu_GeV":energy(t),
                "w":weight(running(xs,t,ell))},
            "arbitrary_decoupling_scale_max_abs_difference":max(above_d)}


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    if isinstance(a,(int,float)):
        return math.isclose(a,b,rel_tol=2e-9,abs_tol=3e-14)
    return a==b


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    opts=parser.add_mutually_exclusive_group(required=True)
    opts.add_argument("--write-new",action="store_true")
    opts.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=compute()
    target=HERE/"results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
    else:
        assert same(result,read("results.json"))
        print("1376: thresholds, matching-scale invariance and same-scale spectra verified")
