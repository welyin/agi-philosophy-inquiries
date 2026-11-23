"""Bound the common charged-scalar readout in two neutrino-mass mechanisms."""
import argparse
from itertools import product
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def interval(values, widths):
    return (values[0] - widths*values[1], values[0] + widths*values[-1])


def electron_row(solar, reactor):
    return np.array([(1-solar)*(1-reactor), solar*(1-reactor), reactor])


def limits(items):
    values = list(items)
    return [float(min(values)), float(max(values))]


def solve():
    data = read(BASE/"inputs.json")
    count = data["widths"]
    splitting = read(BASE/data["splittings"])["juno_207_days"]
    t = interval(data["solar_sin2"], count)
    amp = read(BASE/data["reactor_2theta"])["angles"]["reactor_amplitude"]
    def invert(z):
        return (1-math.sqrt(1-z))/2
    z = tuple(map(invert, interval(amp, count)))
    d = interval(splitting["solar_21"], count)
    a = interval(splitting["normal_atmospheric_31"], count)
    a_inv = interval(splitting["inverted_abs_31"], count)
    f = [(1-electron_row(s,b)[0])/2 for s,b in product(t,z)]
    corners = list(product(t,z,d,a))
    shared, separate, floor_beta_sq = [], [], []
    for s,b,ds,atm in corners:
        p = electron_row(s,b)
        q = p[1]*ds+p[2]*atm
        floor_beta_sq.append(q)
        single = (1-p[0])/2
        for fraction in limits(f):
            separate.append((fraction*(ds+atm)-q)/(1-3*fraction))
        x = (p[1]-p[2])*(atm-ds)/(2-3*(p[1]+p[2]))
        shared.append(x)
        assert math.isclose((x+q)/(3*x+ds+atm),single)
    def crossing(xs):
        xl,xh = limits(xs)
        # Conservative enclosures; endpoints need not use identical input corners.
        return {"lightest_eV":list(map(math.sqrt,(xl,xh))),
                "sum_eV":[np.sqrt([v,v+solar,v+atm]).sum()
                          for v,solar,atm in ((xl,d[0],a[0]),(xh,d[1],a[1]))],
                "beta_eV":[math.sqrt(xl+min(floor_beta_sq)),math.sqrt(xh+max(floor_beta_sq))]}
    io_singlet = limits((1-b)/2 for b in z)
    io_zero = []
    gap = []
    for s,b,ds,atm in product(t,z,d,a_inv):
        p = electron_row(s,b)
        frac = ((1-b)*atm+p[1]*ds)/(2*atm+ds)
        io_zero.append(frac)
        gap.append((1-b)/2-frac)
    # Exact witness with identical mu and tau rows and the same oscillation inputs.
    mean = electron_row(data["solar_sin2"][0],invert(amp[0]))
    labels = np.vstack((mean,(1-mean)/2,(1-mean)/2))
    offset = np.array([0.,splitting["solar_21"][0],splitting["normal_atmospheric_31"][0]])
    single = (1-labels[:,0])/2
    x = (single[0]*offset.sum()-mean@offset)/(1-3*single[0])
    triple = labels@(x+offset)/(3*x+offset.sum())
    np.testing.assert_allclose(single,triple,atol=1e-15)
    mass_sum = float(np.sqrt(x+offset).sum())
    assert mass_sum > data["cosmology"]["sum_limits_eV"][0]
    assert mass_sum < data["cosmology"]["sum_limits_eV"][1]
    return {"NO_singlet_e":limits(f),"NO_triplet_massless_e":limits(
            electron_row(s,b)@np.array([0,ds,atm])/(ds+atm) for s,b,ds,atm in corners),
            "NO_equal_e_shared_inputs":crossing(shared),"NO_equal_e_separate_inputs":crossing(separate),
            "NO_singlet_beta_eV":list(map(math.sqrt,limits(floor_beta_sq))),
            "NO_exact_three_flavour_witness":{"triplet_lightest_eV":math.sqrt(x),
                "triplet_sum_eV":mass_sum,"fractions":single.tolist()},
            "IO_singlet_e":io_singlet,"IO_max_triplet_e_at_zero":max(io_zero),
            "IO_min_gap_shared_inputs":min(gap),"IO_min_gap_separate_inputs":io_singlet[0]-max(io_zero)}


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    mode = cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = cli.parse_args()
    result = solve()
    out = BASE/"results.json"
    if args.write_new:
        with out.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    else:
        assert result == read(out)
        print("1485: profiled flavour fractions, ordering gaps and absolute-mass crossing verified.")
