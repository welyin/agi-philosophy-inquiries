"""Predict the light down-type Yukawa using an identified Hermitian port."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"))


def box(row, n):
    value, *errors = row
    half = n*sum(errors)
    return (value-half, value+half)


def overlap(first, second):
    return bool(max(first[0], second[0]) <= min(first[1], second[1]))


def probabilities(s12, s13):
    return np.array([(1-s12*s12)*(1-s13*s13), s12*s12*(1-s13*s13), s13*s13])


def scalar_outputs(down, p):
    return {tag: abs((p[1]*down[1]+sign*p[2]*down[2])/p[0])
            for tag, sign in [("equal_sb", 1), ("opposite_sb", -1)]}


def fit_box(rows, angles, n, sine_input):
    down_target, strange, bottom = [box(x,n) for x in rows]
    theta12, theta13 = [box(x,n) for x in angles[:2]]
    vals = {key: [] for key in ("equal_sb", "opposite_sb")}
    for s, b, a, t in itertools.product(strange, bottom, theta12, theta13):
        sine12, sine13 = (a,t) if sine_input else np.sin([a,t])
        for key, val in scalar_outputs([0,s,b], probabilities(sine12,sine13)).items():
            vals[key].append(val)
    return {key: {"prediction":[min(v),max(v)],
                  "overlap":overlap([min(v),max(v)], down_target)} for key,v in vals.items()}


def named_sources(point, labels, target):
    v = dict(zip(labels, point["values"]))
    y = np.array([v[k] for k in ["yd","ys","yb"]])
    s12,s13 = np.sin([v["theta12_q"],v["theta13_q"]])
    p = probabilities(s12,s13)
    # A full unitary basis with this first row suffices to test T=T† and <u|T|u>.
    first = np.sqrt(p)
    ortho = np.array([-first[1],first[0],0.])
    ortho /= np.linalg.norm(ortho)
    mixing = np.vstack([first,ortho,np.cross(first,ortho)])
    out = {}
    for signs in itertools.product((-1,1), repeat=2):
        signed_masses = y*np.array([1,*signs])
        cof = np.prod(signed_masses)/signed_masses
        c = -(p@signed_masses)/(p@cof)
        post = signed_masses+c*cof
        mass = mixing@np.diag(signed_masses)@mixing.T
        poly = np.linalg.det(mass)*np.linalg.inv(mass)
        full = mass+c*poly
        assert np.allclose(np.linalg.eigvalsh(full),sorted(post),atol=3e-15,rtol=1e-13)
        assert abs(full[0,0]) < 2e-17
        assert np.all(np.diff(abs(post)) > 0)
        # For fixed observed |first mass| minimize the return residual. It is affine
        # in the signed first output and its sole zero is given above.
        residuals = []
        for isign in (-1,1):
            for dmeas in target:
                c_target = (isign*dmeas-signed_masses[0])/cof[0]
                residuals.append(abs(p@(signed_masses+c_target*cof)))
        assert not overlap(target, [abs(post[0]),abs(post[0])])
        out[str(signs)] = {"c":c,"signed_after":post.tolist(),
                           "minimum_abs_return_for_target":min(residuals)}
    return out


def low_energy_tree(cfg,n):
    tree = load(cfg["tree"])
    belle = load(cfg["belle"])
    qcd = cfg["lattice"]
    lf, uf = box(tree["K_form"],n)
    lp, up = box(tree["K_product"],n)
    v_us = [lp/uf,up/lf]
    cb = box(belle["BelleII_Vcb_BGL"],n)
    rho = box(tree["ub_over_cb"],n)
    v_ub = [cb[0]*rho[0],cb[1]*rho[1]]
    sl = box(qcd["strange_to_light"],n)
    ud = box(qcd["up_to_down"],n)
    ds = [2/((1+ud[1])*sl[1]),2/((1+ud[0])*sl[0])]
    b_over_s = box(qcd["bottom_to_strange"],n)
    pred = {key:[] for key in ("equal_sb","opposite_sb")}
    for s,b,ratio in itertools.product(v_us,v_ub,b_over_s):
        p = [1-s*s-b*b,s*s,b*b]
        for key,val in scalar_outputs([0,1,ratio],p).items():
            pred[key].append(val)
    return {"lattice_d_over_s":ds,
            "predictions":{key:[min(v),max(v)] for key,v in pred.items()},
            "overlap":{key:overlap([min(v),max(v)],ds) for key,v in pred.items()}}


def compute():
    cfg = load("inputs.json")
    hi = load(cfg["high"])["high_scale"]
    lo = load(cfg["low"])
    tables = {}
    for name, rows, ang, sine in [
            ("MZ",lo["yukawa"]["down"],lo["quark_sines_and_phase"],True),
            ("higher",hi["down"],hi["angles_theta12_theta13_theta23_delta"],False)]:
        tables[name] = {str(n):fit_box(rows,ang,n,sine) for n in cfg["widths"]}
    tree = {str(n):low_energy_tree(cfg,n) for n in cfg["widths"]}
    seeds = load(cfg["seeds"])
    fixed = {name:named_sources(seeds["candidates"][name],seeds["labels"],box(hi["down"][0],2))
             for name in cfg["seeds_used"]}
    return {"scale_boxes":tables, "fixed_signed_seeds":fixed, "lattice_tree":tree}


def verify(got, expected):
    if isinstance(got, dict):
        assert got.keys() == expected.keys()
        for k in got:
            verify(got[k],expected[k])
    else:
        np.testing.assert_allclose(got,expected,rtol=7e-11,atol=3e-17)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    sel = cli.add_mutually_exclusive_group(required=True)
    sel.add_argument("--write-new",action="store_true")
    sel.add_argument("--check",action="store_true")
    arg = cli.parse_args()
    result = compute()
    path = HERE/"results.json"
    if arg.write_new:
        payload = json.dumps(result,indent=2,allow_nan=False) + "\n"
        with path.open("x",encoding="utf-8") as f:
            f.write(payload)
    else:
        verify(result,load(path))
        print("1565 signed return, cofactor matching and independent spectra verified.")
