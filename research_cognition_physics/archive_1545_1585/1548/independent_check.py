"""Two-pole flavor budget and independent imperfect-alignment cross-check."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def unitary(n, rng):
    q, _ = np.linalg.qr(rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n)))
    return q


def pair(d, e, left, right):
    h = d[:, None]*left-right*e
    g = d[:, None]*right-left*e
    power = np.vdot(h,h).real
    np.testing.assert_allclose(power,np.vdot(g,g).real,rtol=1e-13)
    return abs(h)**2/power,abs(g)**2/power


def evaluate():
    cfg = load(HERE/"inputs.json")
    ds = np.array(load(HERE/cfg["down_anchor"])["down_yukawa_dr_1TeV"])
    es = np.array(cfg["charged_leptons_DR_1TeV"])
    d, e = ds[:,0], es[:,0]
    rng = np.random.default_rng(821548)
    # Scale each row or column by its smallest denominator for conditioning.
    denom = d[:,None]**2-e**2
    max_null = 0.
    for _ in range(90):
        h, g = pair(d,e,unitary(3,rng),unitary(3,rng))
        c2 = rng.uniform()
        pole1, pole2 = c2*h+(1-c2)*g, (1-c2)*h+c2*g
        signed = (pole1-pole2)/denom
        for axis in (0,1):
            residual = signed.sum(axis=axis)*abs(denom).min(axis=axis)
            max_null = max(max_null, float(np.max(abs(residual))))
    assert max_null < 1e-12
    opt = np.eye(3)[[0,2,1]]
    lower = 1.
    coeff = 1.
    for signs in itertools.product((-2,2), repeat=6):
        dx = d+ds[:,1]*signs[:3]
        ex = e+es[:,1]*signs[3:]
        h,_ = pair(dx,ex,opt,np.diag([-1,-1,1])@opt)
        lower = min(lower,float(h[2,1]))
        coeff = min(coeff,float((dx[2]**2-ex[2]**2)/(dx[2]**2-ex[1]**2)))
    reported = load(HERE/"results.json")
    np.testing.assert_allclose(lower,reported["two_width_sharp_lower"])
    # A unitary right-hand partner cancels the b mu cell of the second pole.
    c = e[1]/d[2]
    s = np.sqrt(1-c*c)
    h, g = pair(d,e,opt,np.array([[1,0,0],[0,c,-s],[0,s,c]])@opt)
    delta = float(h[2,1]-g[2,1])
    assert g[2,1]<1e-25 and np.sum(g[2,[0,2]]) >= coeff*delta
    # A fictitious reassignment within one pole keeps total normalization but
    # fails a sum rule. Use the directly measured bottom-row sign change.
    modified = h.copy()
    modified[2,1] -= 0.010
    modified[0,1] += 0.010
    assert modified.min()>=0 and abs(modified.sum()-1)<1e-13
    mismatch = abs(((modified-g)/denom)[2].sum())*abs(denom[2,1])
    np.testing.assert_allclose(mismatch,0.010)
    return {
        "two_width_exclusive_sharp": lower,
        "two_width_compensating_branch_ratio": coeff,
        "companion_be_or_btau_sum_if_pure_other_bmu_dark_min": coeff*lower,
        "central_dark_example_other_be_btau": float(g[2,[0,2]].sum()),
        "mixed_poles_weighted_null_max": max_null,
        "deliberate_one_percent_reassignment_residual": float(mismatch),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args = parser.parse_args()
    result = evaluate()
    dest = HERE/"independent_results.json"
    if args.write_new:
        with dest.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        saved = load(dest)
        assert saved.keys()==result.keys()
        for name in saved:
            np.testing.assert_allclose(result[name],saved[name],rtol=1e-10,atol=1e-13)
        print("1548: independent saturation, two-pole flavor closure and migration checks passed")
