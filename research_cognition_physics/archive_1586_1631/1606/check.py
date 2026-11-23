"""CP triangle transported by two simultaneous real-axis tags. --write-new or --check."""
import argparse
import cmath
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def box(data, n):
    return [(m - n*w, m + n*w) for m, w in data]


def delta(data):
    return math.prod(data[j]**2-data[i]**2 for i, j in itertools.combinations(range(3), 2))


def delta_bounds(intervals):
    pairs = list(itertools.combinations(range(3), 2))
    low = math.prod(intervals[j][0]**2-intervals[i][1]**2 for i, j in pairs)
    high = math.prod(intervals[j][1]**2-intervals[i][0]**2 for i, j in pairs)
    assert low > 0
    return low, high


def rectangle(down, lep, Jlo, Jhi):
    dl, dh = delta_bounds(down)
    el, eh = delta_bounds(lep)
    return {
        "J_down_lower":Jlo, "J_down_upper":Jhi,
        "I_down_abs_lower":dl*Jlo, "I_down_abs_upper":dh*Jhi,
        "J_right_over_k6_abs_lower":9*dl/eh*Jlo,
        "J_right_over_k6_abs_upper":9*dh/el*Jhi
    }


def untouched_top(down, lep, vtb_lo):
    # The positive spectral polynomial (he-mu²)(he-tau²) is tested on e_t.
    # Its expectation divided by k^0 is bounded above by this quadratic in k².
    a = down[2][1]**4
    b = (lep[1][0]**2 + lep[2][0]**2) * (vtb_lo*down[2][0])**2
    c = (lep[1][1] * lep[2][1])**2
    disc = math.sqrt(b*b - 4*a*c)
    return {"k_excluded_open_interval_from_two_moments":
            [math.sqrt(2*c/(b+disc)), math.sqrt((b+disc)/(2*a))],
            "k_upper_from_first_moment":lep[2][1]/(down[2][0]*vtb_lo),
            "polynomial_upper_k1":a-b+c}


def sine_products(a):
    si = [math.sin(x[0]) for x in a[:3]]
    co = [math.cos(x[1]) for x in a[:3]]
    low = si[0]*si[1]*si[2]*co[0]*co[1]**2*co[2]
    upper = math.prod(math.sin(x[1]) for x in a[:3])
    return low, upper


def triangle(h):
    return float((h[0,1]*h[1,2]*h[2,0]).imag)


def build(mother, tags, gains, clebsch):
    columns = []
    for j in tags:
        v = np.linalg.solve(mother, np.eye(3)[:,j])
        v /= np.linalg.norm(v)
        columns.append(v)
    insert = sum(g*np.outer(v,v) for g,v in zip(gains,columns))
    return mother @ (np.eye(3) + clebsch*insert), columns


def invariants(matrix):
    h = matrix @ matrix.conj().T
    spectrum = np.linalg.svd(matrix,compute_uv=False)[::-1]
    return triangle(h), spectrum


def witness(ins):
    mother = np.array(ins["example_mother"])
    z = ins["example_source_strength"] * cmath.exp(1j*math.radians(ins["example_vacuum_phase_deg"]))
    k = ins["example_k"]
    d, (a,b) = build(mother, [0,2], [z,z.conjugate()], 1)
    e = k * build(mother, [0,2], [z,z.conjugate()], -3)[0]
    c = float(a @ b)
    det_expected = 1+2*z.real+(1-c*c)*abs(z)**2
    det_ratio = np.linalg.det(d)/np.linalg.det(mother)
    assert abs(det_ratio-det_expected) < 1.e-14 and det_ratio.real>0
    idown, yd = invariants(d)
    ilep, yl = invariants(e)
    assert abs(idown)>1e-5 and abs(ilep-9*k**6*idown)<1e-12
    single_phase = build(mother, [0,2], [z,-0.55*z], 1)[0]
    double_phase = build(mother, [0,2], [0.28j,-0.41j], 1)[0]
    assert max(abs(invariants(m)[0]) for m in [single_phase,double_phase]) < 2e-15
    # An independent neutral mass orientation changes PMNS, holding charged masses fixed.
    u_right, _, vdag = np.linalg.svd(e)
    polar = u_right @ vdag
    lepton_left = vdag.conj().T
    neutral_aligned = polar.conj().T
    W = lepton_left.conj().T @ neutral_aligned
    assert np.allclose(W,u_right.conj().T)
    neutral_competitor = lepton_left
    assert np.allclose(lepton_left.conj().T @ neutral_competitor,np.eye(3))
    # The third tagged edge is a genuine operator read by both sectors.
    zd = [z, 0.17+0.09j, z.conjugate()]
    three_d = build(mother, [0,1,2],zd,1)[0]
    three_e = k*build(mother,[0,1,2],zd,-3)[0]
    failure = invariants(three_e)[0] - 9*k**6*invariants(three_d)[0]
    assert abs(failure)>1e-5
    return {"tag_overlap":c,"detD_over_detMother":float(det_ratio.real),
            "theta_down_rad":float(np.angle(det_ratio)), "I_down":idown, "I_charged10":ilep,
            "J_down_abs":abs(idown/delta(yd)), "J_charged10_abs":abs(ilep/delta(yl)),
            "I_common_phase_zero_checked":1, "three_tag_identity_residual":failure}


def compute():
    ins=read("inputs.json")
    n=ins["widths"]
    high=read(ins["high_masses"])["scales"]["1e16"]
    angles=box(read(ins["high_ckm"])["high_scale"]["angles_theta12_theta13_theta23_delta"],n)
    low_factor, high_factor = sine_products(angles)
    hi=rectangle(box(high[:3],n),box(high[3:],n),low_factor*math.sin(angles[3][0]),high_factor)
    tree=read(ins["weak_tree"])["two_width_box"]
    us, cb, ratio=[tree[key][1] for key in ("Vus","Vcb","ratio_ub_cb")]
    ub=cb*ratio
    # Exact J <= |Vus Vcb Vub| (unitary PDG, c12 c23 <= 1).
    weak=read(ins["weak_masses"])["yukawa"]
    low=rectangle(box(weak["down"],n),box(weak["charged_leptons"],n),
                  tree["absolute_J_lower"],us*cb*ub)
    # |Vtb|_tree = sqrt(1-|Vcb|²-|Vub|²) from three-generation unitarity.
    top = {"high":untouched_top(box(high[:3],n),box(high[3:],n),
                               math.cos(angles[1][1])*math.cos(angles[2][1])),
           "weak":untouched_top(box(weak["down"],n),box(weak["charged_leptons"],n),
                               math.sqrt(1-cb**2-ub**2))}
    osc=read(ins["neutrino"])
    theta23=read(ins["atmospheric"])["theta23_sin2"]
    phase=read(ins["phase"])["phase_credible_arcs_pi"]
    predictions={}
    for order in ("NO","IO"):
        _, xlo, xhi = osc[order]["s12"]
        _, zlo, zhi = osc[order]["s13"]
        _, vlo, vhi = theta23[order]
        # sqrt[x(1-x) v(1-v) z]*(1-z) separately bounded from below.
        vprod = min(vlo*(1-vlo),vhi*(1-vhi))
        factor = math.sqrt(xlo*(1-xlo)*vprod*zlo)*(1-zhi)
        predictions[order]= {}
        for label,reg in (("high",hi),("weak",low)):
            jmax=reg["J_right_over_k6_abs_upper"]
            sine_upper=jmax/factor
            row={"leptonic_CP_even_factor_lower":factor,"abs_sin_delta_upper_k1":sine_upper,
                 "distance_to_piZ_degrees_upper_k1":math.degrees(math.asin(sine_upper))}
            if order=="IO":
                guard=ins["IO_rounding_pi"]
                a,b=phase["IO_99.73"]
                floor=min(abs(math.sin(math.pi*(a-guard))),abs(math.sin(math.pi*(b+guard))))
                wanted=factor*floor
                row.update({"IO_phase_abs_sine_lower_guarded":floor,"observed_J_abs_lower":wanted,
                            "necessary_k_for_IO":(wanted/jmax)**(1/6),
                            "required_additive_abs_J_k1":wanted-jmax,
                            "relative_unitary_error_op_lower_k1":(wanted-jmax)/4})
            predictions[order][label]=row
    return {"scales":{"high":hi,"weak":low},"untouched_top":top,"conditional_neutral_readout":predictions,
            "same_tree_strong_and_weak_example":witness(ins)}


def equal(a,b):
    if isinstance(a,dict):
        return isinstance(b,dict) and a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return isinstance(b,list) and len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    if a == 0 or b == 0:
        return abs(a-b) < 1.e-15
    return math.isclose(a,b,abs_tol=0,rel_tol=2.e-10)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    opts=parser.add_mutually_exclusive_group(required=True)
    opts.add_argument("--write-new",action="store_true")
    opts.add_argument("--check",action="store_true")
    args=parser.parse_args()
    cert=compute()
    if args.write_new:
        with (HERE/"results.json").open("x", encoding="utf-8") as stream:
            json.dump(cert,stream,ensure_ascii=False,indent=2,allow_nan=False)
            stream.write("\n")
    elif not equal(cert,read("results.json")):
        raise SystemExit("1606 comparison failed.")
    print("1606: two-tag CP transport, charged and neutral conditions checked")
