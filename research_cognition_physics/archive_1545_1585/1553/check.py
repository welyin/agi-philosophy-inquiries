"""Right-handed doubly charged dilepton fractions from a common symmetric mass tensor."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def fractions(leptons, down):
    f = down - np.diag(leptons)
    norm = np.vdot(f, f).real
    return abs(f[2, 2])**2/norm, 2*np.vdot(f[:2, 2], f[:2, 2]).real/norm


def near_edge(row, z):
    d, s, b, e, mu, tau = row
    assert s-1e-15 <= z <= b+1e-15
    cosine = math.sqrt(max(0, (z-s)/(b-s)))
    sine = math.sqrt(max(0, 1-cosine**2))
    r = np.array([[1, 0, 0], [0, cosine, sine], [0, -sine, cosine]])
    u = np.diag([1j,1j,1]) @ r
    y = u @ np.diag([d,s,b]) @ u.T
    np.testing.assert_allclose(np.linalg.svd(y, compute_uv=False), [b,s,d])
    return fractions([e,mu,tau], y)


def exact(row):
    d, s, b, e, mu, tau = row
    a = tau+mu
    l = (b+mu)**2 + (s+mu)**2 + (d+e)**2
    k = l + tau**2-mu**2
    if tau > b:
        z = l/a - mu
        assert s <= z <= b
        min_same = 1-l/a**2
        min_mixed = None
        np.testing.assert_allclose(near_edge(row,z)[0], min_same, atol=1e-13)
    else:
        assert s < tau <= b
        z = tau
        min_same = 0.
        min_mixed = 2*(b-tau)*(tau-s)/(k-2*a*tau)
        same, mixed = near_edge(row,z)
        assert same < 1e-25
        np.testing.assert_allclose(mixed,min_mixed,atol=1e-12)
    return min_same, min_mixed


def box(rows, mult, seed):
    x = np.asarray(rows)
    lo, hi = x[:,0]-mult*x[:,1], x[:,0]+mult*x[:,1]
    center = x[:,0]
    cen_same, cen_mixed = exact(center)
    out = {"central_min_same":cen_same, "central_min_mixed_if_same_zero":cen_mixed}
    # E = W diag(e,mu,tau) after putting only the right-handed masses on physical axes.
    # D-E remains symmetric. The two nonzero singular values of 3x3 skew D-D.T agree.
    skew_ratio = (hi[2]+hi[1])/(lo[5]-hi[4])
    abs_tau = max(0., lo[5]*math.sqrt(max(0.,1-skew_ratio**2))-hi[2])
    out["symmetric_difference_only_two_width_same_floor"] = abs_tau**2/np.sum((hi[:3]+hi[3:])**2)
    if lo[5] > hi[2]:
        d, s, b, e, mu, tau = lo
        D, S, B, E, MU, TAU = hi
        # z+mu=L/(tau+mu); ensure stationary point is attainable throughout.
        assert (b+mu)**2+(s+mu)**2+(d+e)**2 >= (S+MU)*(TAU+MU)
        assert (B+MU)**2+(S+MU)**2+(D+E)**2 <= (b+mu)*(tau+mu)
        assert (b+mu)*(tau-B)+(s+mu)*(tau-S) > (D+E)**2
        worst = np.array([D,S,B,E,MU,tau])
        beta = exact(worst)[0]
        out["two_width_min_same_exact"] = beta
        ref = 1/6
        out["relative_extra_tensor_to_six_equal_at_least"] = max(0.,math.sqrt(beta*(1-ref))-math.sqrt((1-beta)*ref))
    else:
        assert lo[2] > hi[5] and lo[5] > hi[1]
        d,s,b,e,mu,tau = lo
        D,S,B,E,MU,TAU = hi
        # Denominator K-2(tau+mu)tau = L-(tau+mu)^2.
        den = (B+MU)**2+(S+MU)**2+(D+E)**2-(tau+mu)**2
        out["two_width_min_same_exact"] = 0.
        out["two_width_mixed_if_same_zero_certified_floor"] = 2*(b-TAU)*(tau-S)/den
        eps = .01*np.linalg.norm(hi[:3]+hi[3:])
        coefficient = max(abs(B+S-2*tau),abs(b+s-2*TAU))
        numerator = (b-TAU)*(tau-S)-coefficient*eps-eps**2
        out["two_width_mixed_if_same_at_most_one_in_10000"] = max(0.,2*numerator/(den+2*(TAU+MU)*eps))
    for corner in itertools.product((0,1),repeat=6):
        row = np.where(corner,hi,lo)
        same, mixed = exact(row)
        assert same >= out["two_width_min_same_exact"]-1e-12
        if mixed is not None:
            assert mixed >= out["two_width_mixed_if_same_zero_certified_floor"]-1e-12
    rng = np.random.default_rng(seed)
    for _ in range(80):
        u, _ = np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
        same, _ = fractions(center[3:], u @ np.diag(center[:3]) @ u.T)
        assert same >= cen_same-1e-13
    return out


def evaluate():
    config = read(HERE/"inputs.json")
    low = read(HERE/config["low_input"])["light_scales"]
    high = read(HERE/config["upper_input"])["scales"]
    chosen = [(k,low[k]) for k in config["low_scales"]]
    chosen += [(k,high[k]) for k in config["upper_scales"]]
    results = {k:box(rows,config["widths"],config["random_seed"]) for k,rows in chosen}
    # Direct Poisson upper limit on the accepted signal, using no background subtraction.
    a,b = 0.,12.
    for _ in range(80):
        z = (a+b)/2
        if math.exp(-z)*(1+z) > .05:
            a = z
        else:
            b = z
    beta = results["1e16"]["two_width_min_same_exact"]
    results["DELPHI_prompt_90GeV_conditional"] = {
        "accepted_signal_events_one_sided_upper": b,
        "reference_detector_eff_at_206p7": .448,
        "working_efficiency_floor": .40, "tau_shape_transfer_floor": .99,
        "produced_pairs_times_Bll_squared_upper_at_206p7": b/(.40*(.99*beta)**2)}
    return results


def compare(out,saved):
    assert set(out) == set(saved)
    for key,val in out.items():
        if isinstance(val, dict):
            compare(val,saved[key])
        elif val is None:
            assert saved[key] is None
        else:
            np.testing.assert_allclose(val,saved[key],rtol=1e-10,atol=1e-14)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args = parser.parse_args()
    result = evaluate()
    path = HERE/"results.json"
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        compare(result,read(path))
    print("PASS 1553 charged spectra, sharp tau fractions, two-width windows and conditional four-tau count")
