"""Charged-row preserving singlet completions and neutral-charge tests."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def root_on_interval(f, low, high):
    assert f(low) >= 0 >= f(high)
    for _ in range(70):
        mid = (low+high)/2
        if f(mid) > 0:
            low = mid
        else:
            high = mid
    return (low+high)/2


def rotation(i, j, p, size=4):
    mat = np.eye(size)
    mat[i,i] = mat[j,j] = math.sqrt(1-p)
    mat[j,i] = math.sqrt(p)
    mat[i,j] = -mat[j,i]
    return mat


def with_first_row(row):
    row = np.array(row, dtype=float)
    assert abs(np.linalg.norm(row)-1) < 2e-13
    _, _, right = np.linalg.svd(row[None, :], full_matrices=True)
    return np.vstack((row, right[1:]))


def complete(charged, deficit, x, down_index):
    obs = np.asarray(charged)
    y = (deficit-x)/(obs[down_index]**2 + deficit-x)
    bare = obs/math.sqrt(1-x)
    bare[down_index] /= math.sqrt(1-y)
    w = np.zeros((4,4))
    w[:3,:3] = with_first_row(bare)
    pu, pd = rotation(0,3,x), rotation(down_index,3,y)
    w = pu.T @ w @ pd
    expected_u = pu.T @ np.diag([1.,1.,1.,0.]) @ pu
    expected_d = pd.T @ np.diag([1.,1.,1.,0.]) @ pd
    assert np.allclose(w[0,:3], obs, atol=1e-13, rtol=0)
    assert np.allclose(w @ w.T, expected_u) and np.allclose(w.T @ w, expected_d)
    return float(y)


def certificate():
    a = read(HERE / "inputs.json")
    ch = read(HERE / a["charged_inputs"])
    history = read(HERE / a["charged_result"])
    charges = [ch["ud_beta"][0], ch["K_product"][0]/ch["K_form"][0],
               ch["ub_over_cb"][0]*ch["F_times_cb"][0]/ch["F_example"]]
    d = 1-sum(i*i for i in charges)
    t = charges[0]**2
    assert math.isclose(1-d, history["row_nominal"], abs_tol=1e-12)
    p, cs = a["proton"], a["caesium"]
    ep = math.hypot(p["pv_electron_scattering_only"][1],p["SM"][1])
    ec = math.hypot(*cs["value_and_errors_experiment_atom"][1:],cs["SM"][1])
    op = p["pv_electron_scattering_only"][0] - p["SM"][0]
    oc = cs["value_and_errors_experiment_atom"][0] - cs["SM"][0]
    z, n = cs["Z"], cs["N"]
    l, r = 2*z+n, z+2*n

    def shifts(x):
        y = (d-x)/(t+d-x)
        proton, nucleus = -2*x+y, -l*x+r*y
        assert abs((nucleus-r*proton)/(3*n)-x) < 1e-14
        assert abs((2*nucleus-l*proton)/(3*n)-y) < 1e-14
        assert abs(x+t*y/(1-y)-d) < 1e-14
        return y, proton, nucleus

    def report(x):
        y, sp, sc = shifts(x)
        assert math.isclose(complete(charges,d,x,0),y,abs_tol=1e-13)
        return {"up_leak":x, "down_leak":y, "proton":sp, "caesium":sc,
                "residual_to_historical_in_widths":[(sp-op)/ep,(sc-oc)/ec]}

    x_proton_zero = root_on_interval(lambda x:shifts(x)[1],0,d)
    x_minimax = root_on_interval(lambda x:shifts(x)[1]/ep+shifts(x)[2]/ec,0,d)
    _, sp, sc = shifts(x_minimax)
    scale = max(abs(sp)/ep,abs(sc)/ec)
    # At minimax y = k*x. Bounding d below and |Vud|^2 above bounds the
    # optimum for every point of the earlier one-width charged-current box.
    k = (2*ec+l*ep)/(ec+r*ep)
    smallest_d = 1 - history["bounds"]["1"]["official_mixed"]["row"][1]
    largest_t = (ch["ud_beta"][0] + ch["ud_beta"][1])**2
    x_box = root_on_interval(lambda x:smallest_d-x-largest_t*k*x/(1-k*x),0,smallest_d)
    box_scale = (2-k)*x_box/ep
    assert math.isclose((2-k)*x_minimax/ep,scale,rel_tol=2e-12)

    def actual_box(k):
        lo, hi = 0., d
        for i, center, width in ((1,op,ep),(2,oc,ec)):
            func = lambda x:shifts(x)[i]
            if func(lo) < center-k*width or func(hi)>center+k*width:
                return []
            if func(lo)>center+k*width:
                lo=root_on_interval(lambda x:func(x)-center-k*width,lo,hi)
            if func(hi)<center-k*width:
                hi=root_on_interval(lambda x:func(x)-center+k*width,lo,hi)
        return [lo,hi]

    # A true fourth sequential doublet has T3 proportional to the identity per charge.
    unitary = with_first_row([*charges, math.sqrt(d)])
    assert np.allclose(unitary @ unitary.T, np.eye(4), atol=1e-13)
    assert np.allclose(unitary.T @ unitary, np.eye(4), atol=1e-13)
    return {"visible_row":charges,"D":d,"Qp_work_width":ep,"QCs_work_width":ec,
            "data_minus_baseline":[op,oc], "pure_u":report(d),"pure_d":report(0),
            "proton_blind":report(x_proton_zero), "null_minimax_point":report(x_minimax),
            "common_width_multiplier_for_null":scale,"improvement_factor":1/scale,
            "required_null_widths":[ep*scale,ec*scale],
            "one_charged_box_null_guarantee":{"D_lower":smallest_d,"ud_squared_upper":largest_t,
                "width_scale_lower":box_scale,"sufficient_improvement":1/box_scale},
            "pure_s_with_both_valence_targets_blind":complete(charges,d,0,1),
            "sequential_heavy_column_amplitude":float(unitary[0,3]),
            "admissible_x_for_historical_centered_boxes":{"one":actual_box(1),"two":actual_box(2)}}


def agrees(a,b):
    if isinstance(a,dict):
        return isinstance(b,dict) and a.keys()==b.keys() and all(agrees(v,b[k]) for k,v in a.items())
    if isinstance(a,list):
        return isinstance(b,list) and len(a)==len(b) and all(agrees(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=2e-11,abs_tol=2e-12) if isinstance(a,float) else a==b


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write",action="store_true")
    out = certificate()
    path = HERE / "results.json"
    if ap.parse_args().write:
        with path.open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2)
            f.write("\n")
    else:
        assert agrees(out,read(path))
        print("1385: weak generators, complete row, neutral targets and precision verified")
