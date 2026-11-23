"""Interval certificate on all real mass-column branches, independent of the angle solver."""

import argparse
import itertools
import json
import math
from decimal import Decimal
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "independent_results.json"
ONE = F(1)
GRID = 10**23


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8"), parse_float=Decimal)


def sqrt_pair(x):
    lower = math.isqrt(x.numerator * GRID**2 // x.denominator)
    return F(lower, GRID), F(lower + 1, GRID)


def mass_data():
    rows = read("../1402/inputs.json")["scales"]["1e16"]
    d, e = ([F(row[0]) for row in rows[a:a+3]] for a in (0, 3))
    t = [F(0), max(d[0], e[0])*min(d[1], e[1]),
         max(d[1], e[1])*min(d[2], e[2]), (F(6, 5)*max(d[2], e[2]))**2]

    def residues(s):
        return [math.prod(p-x*x for x in s) / math.prod(p-q for q in t if q != p) for p in t]

    wd, we = residues(d), residues(e)
    ratio = math.prod(d) / math.prod(e)
    assert sum(wd) == sum(we) == 1 and ratio**2 * we[0] == wd[0]
    omega = [ratio**2 * we[j] / wd[j] for j in (1, 2, 3)]
    squares = []
    for y in d:
        row = [p*w/(p-y*y)**2 for p, w in zip(t[1:], wd[1:])]
        squares.append([x/sum(row) for x in row])
    return omega, squares


def x_intervals(w, cmin):
    a, b, c = 9-4*w, 6+4*w, 1-w
    small, large = [], []
    for cp in (cmin, ONE):
        radlo, radhi = sqrt_pair((b*cp)**2-4*a*c)
        for den in (b*cp+radlo, b*cp+radhi):
            small.append(2*c/den)
            large.append(den/(2*a))
    return [(min(xs), max(xs)) for xs in (small, large)]


def intrinsic_bounds(xs, s2, phi_hi):
    cmin = 1-phi_hi**2/2
    lower = []
    for pair in xs:
        vals = [1+2*x*c for x in pair for c in (cmin, ONE)]
        assert min(vals)*max(vals) > 0
        lower.append(min(abs(v) for v in vals))
    d = F(0)
    for i, j in itertools.combinations(range(3), 2):
        dx = max(abs(x-y) for x in xs[i] for y in xs[j])
        sn = min(ONE, 2*dx*phi_hi/(lower[i]*lower[j]))
        d += s2[i]*s2[j]*sn**2
    assert 4*d < 1
    return sqrt_pair(2*d/(1+sqrt_pair(1-4*d)[0]))[1]


def weak_box(widths, transported=False):
    inp = read("../1381/inputs.json")
    def box(data):
        middle, *errors = map(F, data)
        error = widths * sum(errors)
        return middle-error, middle+error

    kl, kh = box(inp["K_product"])
    fl, fh = box(inp["K_form"])
    ul, uh = kl/fh, kh/fl
    cl, ch = box(read("../1412/inputs.json")["BelleII_Vcb_BGL"])
    rl, rh = box(inp["ub_over_cb"])
    if transported:
        ul, uh = F(998,1000)*ul, F(1002,1000)*uh
        ch *= F(115,100)
    bl, bh = cl*rl, ch*rh
    udmin = sqrt_pair(1-uh**2-bh**2)[0]
    s23max = ch / sqrt_pair(1-bh**2)[0]
    cdmin = ul*sqrt_pair(1-s23max**2)[0]-s23max*bh
    gamma = read("../1490/inputs.json")["chargedB"]
    low = F(gamma["central_degrees"])-widths*F(gamma["lower_error"])-F(gamma["additional_degrees"])
    high = F(gamma["central_degrees"])+widths*F(gamma["upper_error"])+F(gamma["additional_degrees"])
    assert 0 < low < high < 90
    z = low * F(333, 106*180)
    sinlo = z-z**3/6+z**5/120-z**7/5040
    jlo = udmin*bl*cdmin*cl*sinlo
    return uh, ch, bh, jlo, cl, rh, udmin*cdmin*rl*sinlo


def bound(degrees, omega, squares, box, scaled=False):
    phi_hi = degrees * F(22, 7*180)
    roots = [x_intervals(w, 1-phi_hi**2/2) for w in omega]
    uh, ch, bh, jlo, cl, rh, normalized_lo = box
    sensitivity = ch*(uh+bh) + uh*bh*(1+ch)
    slope_b, slope_s = uh*(1+rh)/cl, rh*(1+uh)
    if scaled:
        jlo, sensitivity = normalized_lo, slope_b+slope_s
    rows = []
    for indices in itertools.product((0, 1), repeat=3):
        selected = [ri[k] for ri, k in zip(roots, indices)]
        es, eb = [intrinsic_bounds(selected, s2, phi_hi) for s2 in squares[1:]]
        ju = ch*(uh*eb+bh*es)+uh*bh*(eb+ch*es)
        if scaled:
            ju = slope_b*eb+slope_s*es
        rows.append(("".join(map(str, indices)), ju, max(F(0), (jlo-ju)/sensitivity)))
    return rows


def certificate():
    omega, squares = mass_data()
    all_boxes = {}
    for widths in (1, 2):
        box = weak_box(widths)
        windows = {}
        for deg in map(F, ("0.5", "1", "1.5", "2", "3")):
            rows = bound(deg, omega, squares, box)
            worst = max(rows, key=lambda x: x[1])
            windows[str(float(deg))] = dict(
                max_J_upper=float(worst[1]), worst_branch=worst[0],
                positive_gap_branch_count=sum(j < box[3] for _, j, _ in rows),
                minimum_complex_distance=float(worst[2]),
                branches={i: dict(J_upper=float(j), required_distance=float(r)) for i,j,r in rows},
            )
        # Locate the largest certified radius on a 0.001-degree publication grid.
        # Every accepted radius covers the WHOLE interval starting at zero.
        lo, hi = 0, 4000
        while hi-lo > 1:
            mid = (lo+hi)//2
            if max(j for _, j, _ in bound(F(mid, 1000), omega, squares, box)) < box[3]:
                lo = mid
            else:
                hi = mid
        all_boxes[str(widths)] = dict(J_required=float(box[3]), windows=windows,
                                    all_branch_radius_degrees=lo/1000)
    assert all_boxes["2"]["windows"]["1.0"]["positive_gap_branch_count"] == 8
    # Shared |Vcb| cancels from the lower comparison after normalizing the quartet.
    # This also tests finite common transport of the two b-current moduli.
    normalized = {}
    for widths, transfer in itertools.product((1, 2), (False, True)):
        box = weak_box(widths, transfer)
        target = box[-1]
        lo, hi = 0, 5000
        while hi-lo > 1:
            mid = (lo+hi)//2
            rows = bound(F(mid,1000), omega, squares, box, scaled=True)
            if max(j for _, j, _ in rows) < target:
                lo = mid
            else:
                hi = mid
        at_one = max(bound(F(1), omega, squares, box, scaled=True), key=lambda r:r[1])
        normalized[f"widths_{widths}_transport_{transfer}"] = dict(
            all_branch_radius_degrees=lo/1000, required_J_over_cb_squared=float(target),
            one_degree_upper=float(at_one[1]), one_degree_required_distance=float(at_one[2]),
            one_degree_worst_branch=at_one[0])
    return dict(parent="1515 central spectrum at 1e16; all eight real-root assignments",
                root_key="0 smaller magnitude, 1 larger magnitude",
                phase_interval="both signs and zero included",
                widths=all_boxes, J_normalized=normalized)


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = certificate()
    if args.write:
        with OUTPUT.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
            stream.write("\n")
    elif json.loads(OUTPUT.read_text(encoding="utf-8")) != result:
        raise SystemExit("saved interval certificate disagrees")
    print({k: (v["J_required"], v["all_branch_radius_degrees"]) for k, v in result["widths"].items()})


if __name__ == "__main__":
    main()
