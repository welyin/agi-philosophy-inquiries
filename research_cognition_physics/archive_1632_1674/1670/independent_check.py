"""Verify the two projected IceCube tests by SVD and a rectangular outer cover."""
import argparse
import hashlib
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ARCH = HERE.parent
SPEC = runpy.run_path(str(ARCH / "1661/independent_check.py"))


def read(name):
    return json.loads((ARCH / name).read_bytes())


def load_grid():
    first = read("1555/inputs.json")["icecube"]
    rest = read("1670/inputs.json")
    xs = np.array(first["grid_x"] + rest["additional_x"])
    ms = np.array(first["mass_squared_eV2"])
    qs = np.array(first["q_rows"] + rest["additional_q_rows"])
    assert np.all(np.diff(xs) > 0) and np.all(np.diff(ms) > 0)
    assert qs.shape == (39, 29)
    return xs, ms, qs


def at(grid, angle, mass):
    xs, ms, qs = grid
    row = np.array([np.interp(mass, ms, q) for q in qs])
    return float(np.interp(angle, xs, row))


def rectangle(grid, xlim, mlim):
    xs, ms, qs = grid
    def clip(axis, limits):
        a, b = limits
        assert axis[0] <= a <= b <= axis[-1]
        return [a, *axis[(axis > a) & (axis < b)], b]
    vals = [at(grid, x, m) for x in clip(xs, xlim) for m in clip(ms, mlim)]
    xids = np.flatnonzero((xs[:-1] <= xlim[1]) & (xs[1:] >= xlim[0]))
    mids = np.flatnonzero((ms[:-1] <= mlim[1]) & (ms[1:] >= mlim[0]))
    verts = qs[xids[0]:xids[-1]+2, mids[0]:mids[-1]+2]
    return min(vals), float(verts.min()), float(verts.max())


def prob(e, phases, sign, ratio):
    legs = e[:, :2] * np.exp(1j*np.array([0, *phases]))[:, None]
    heavy = SPEC["heavy"](legs, sign, ratio, SPEC["roots"](ratio)[0])
    p1 = abs(SPEC["NULL"] @ legs[:, 0])**2
    return p1, heavy[0], heavy[1]


def check_branch(e, name, sign, phases, grid, angles, solar, reactor, saved):
    lo, hi = saved["with_JUNO_two_widths"]
    edges = np.linspace(lo - 1.e-9, hi + 1.e-9, 2049)
    mid = (edges[:-1] + edges[1:])/2
    p1, pe, pm = prob(e, phases, sign, mid)
    # Round 1661 bounds each heavy probability derivative by 2 on the parent box.
    pad = (edges[1] - edges[0]) + 3.e-9
    xmin = min((pm-pad)/(1-pe+pad))
    xmax = max((pm+pad)/(1-pe-pad))
    # d32 increases with R at each fixed d21 or d31.
    mlo = max(solar[0]*(edges[0]**2-1), reactor[0]*(1-edges[0]**-2)) - 1.e-12
    mhi = min(solar[1]*(edges[-1]**2-1), reactor[1]*(1-edges[-1]**-2)) + 1.e-12
    lower, corner_min, corner_max = rectangle(grid, (xmin, xmax), (mlo, mhi))
    fixed = saved["witness"]
    r = fixed["R"]
    ds = fixed["d21"]
    p1, p3, m3 = prob(e, phases, sign, np.array([r]))
    x = float((m3/(1-p3))[0])
    three = [float((4*p3*(1-p3))[0]), float((1-p1/(1-p3))[0]), x]
    assert all(lower <= v <= upper for v, (lower, upper) in zip(three, angles))
    assert solar[0] <= ds <= solar[1] and reactor[0] <= r*r*ds <= reactor[1]
    witness_mass = (r*r-1)*ds
    qfixed = at(grid, x, witness_mass)
    # Cross-check 1668's angles against the alternate, SVD-based projection.
    np.testing.assert_allclose(three, fixed["three_angles"], rtol=0, atol=1.e-9)
    return dict(R_bounds=[lo, hi], sin2theta23_outer=[xmin, xmax],
                dm32_outer=[mlo, mhi], full_rectangle_q_lower=lower-1.e-7,
                cell_vertex_range=[corner_min, corner_max], explicit_q=qfixed,
                explicit_R=r, explicit_d21=ds, explicit_d32=witness_mass,
                explicit_three_angles=three)


def main():
    grid = load_grid()
    e, _ = SPEC["OLD"]["charged"](read("1657/independent_results.json")["extended_box_GJ_test"])
    angles = read("1659/results.json")["joint_windows"]
    cfg = read("1668/inputs.json")
    solar = SPEC["window"](cfg["juno_NO_dm21"])
    reactor = SPEC["window"](cfg["juno_NO_dm31"])
    prev = read("1668/results.json")["branches"]
    phi = -math.pi/6
    cases = [
        ("third_leg_3_to_2", 1, (2*phi, 3*phi + np.angle(3+2*np.exp(4j*phi)))),
        ("second_leg_counterpart", -1, (10*phi + math.radians(11.5), 7*phi)),
    ]
    out = {name: check_branch(e, name, sign, phases, grid, angles, solar, reactor, prev[name])
           for name, sign, phases in cases}
    third, second = out.values()
    assert third["cell_vertex_range"][1] < 1 and second["cell_vertex_range"][0] > 4
    # A post hoc counterexample to extending the second-leg exclusion to all
    # real two-path coefficients. Its rational ratio has no fitted precision.
    ratio, ds = 5.8, 7.35e-5
    psi = np.angle(2 - np.exp(4j*phi))
    p1, pe, pm = prob(e, (10*phi+psi, 7*phi), -1, np.array([ratio]))
    alt = np.array([4*pe*(1-pe), 1-p1/(1-pe), pm/(1-pe)])[:, 0]
    d31 = ds*ratio**2
    assert all(lower < val < upper for val, (lower, upper) in zip(alt, angles))
    assert solar[0] < ds < solar[1] and reactor[0] < d31 < reactor[1]
    val = at(grid, alt[2], d31-ds)
    assert val < 1
    out["modified_second_leg"] = dict(amplitudes=[2, -1], R=ratio, d21=ds, d31=d31,
                                      d32=d31-ds, psi_degrees=math.degrees(psi),
                                      three_angles=alt.tolist(), q=val)
    # Same phase transfer and R, tighter one-width frequency checks at another d21.
    ds = 7.4e-5
    d31 = ds*ratio**2
    d32 = d31-ds
    for value, key in zip((ds, d31, d32), ("juno_NO_dm21", "juno_NO_dm31", "accelerator_NO_dm32")):
        center, left, right = cfg[key]
        assert center-left < value < center+right
    vals = SPEC["OLD"]["mass_source"](SPEC["roots"](ratio)[0], -1, ds)
    phase = np.exp(-1j*np.array([0, 10*phi+psi, 7*phi]))
    mix = e.conj().T @ (phase[:, None]*vals[0])
    z = abs(mix[0, 2])**2
    s, t = alt[1:]
    j = (mix[0, 0]*mix[1, 1]*mix[0, 1].conjugate()*mix[1, 0].conjugate()).imag
    norm = math.sqrt(s*(1-s)*t*(1-t)*z)
    sd = j/((1-z)*norm)
    cd = (abs(mix[1, 0])**2-s*(1-t)-(1-s)*t*z)/(2*norm)
    assert abs(sd*sd+cd*cd-1) < 1.e-11
    delta = math.atan2(sd, cd)/math.pi
    assert cfg["accelerator_NO_delta_over_pi"][0] < delta < cfg["accelerator_NO_delta_over_pi"][1]
    out["modified_second_leg_one_width"] = dict(R=ratio, d21=ds, d31=d31, d32=d32, delta_over_pi=delta,
                                               q=at(grid, alt[2], d32))
    return out


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check", action="store_true")
    cli.add_argument("--source", type=Path, help="Optionally validate the archived official .npy")
    args = cli.parse_args()
    if args.source:
        sha = hashlib.sha256(args.source.read_bytes()).hexdigest()
        assert sha == "7aadd4ae1ca4624fe6e0909ebf17b45f3f0b8bdc4270d1fa66585fe2c3b5c608"
        raw = np.load(args.source, allow_pickle=True).item()
        xs, ms, qs = load_grid()
        np.testing.assert_array_equal(raw["sin2_theta23"], np.tile(xs, len(ms)))
        np.testing.assert_array_equal(raw["deltam232"], np.repeat(ms, len(xs)))
        np.testing.assert_array_equal(raw["deltallh"].reshape(len(ms), len(xs)).T, qs)
    out = main()
    target = HERE / "independent_results.json"
    if args.check:
        old = json.loads(target.read_bytes())
        for name, result in out.items():
            for key, val in result.items():
                np.testing.assert_allclose(val, old[name][key], rtol=1.e-11, atol=1.e-13)
        print("PASS: independent SVD, complete outer rectangle, grid cells and feasible points")
    else:
        with target.open("x", encoding="utf8") as f:
            f.write(json.dumps(out, indent=2) + "\n")
