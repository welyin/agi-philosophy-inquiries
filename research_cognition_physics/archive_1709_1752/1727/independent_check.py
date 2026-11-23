"""Invert four weak transfer poles, allowing unequal bare masses and a complex direct term."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def invariants(poles):
    m = np.sort(np.asarray(poles, float))
    a = m[1]**2
    # Normalize before factoring the physical characteristic polynomial.
    roots = (m[[0,2,3]]/m[1])**2
    pa = np.prod(1-roots)
    p2a = np.prod(2-roots)
    trace4 = sum((m/m[1])**2)-4
    return a, pa, p2a, trace4


def masses(a,b,p,eps):
    s = np.array([1.,0,0])
    u = np.array([math.sqrt(p),math.sqrt(1-p),0])
    mat = np.zeros((4,5), complex)
    mat[:3,0] = -math.sqrt(a)*s
    mat[:3,1:4] = math.sqrt(a)*np.eye(3)
    mat[3] = math.sqrt(b)*np.r_[eps,-u,1]
    return np.sort(np.linalg.svd(mat,compute_uv=False))


def reconstruct(poles,p):
    a,pa,p2a,tr = invariants(poles)
    b = pa/(1-p)
    e = tr/b-2
    re_z_eps = (-p2a/b-p-e)/2
    return a,b*a,e,re_z_eps,p+e-2*re_z_eps


def run():
    point = read("../1720/results.json")["central_conditional_combined_point"]
    poles = point["four_weak_poles_GeV"]
    a,pa,p2a,tr = invariants(poles)
    b0 = tr/2
    p0 = 1-pa/b0
    np.testing.assert_allclose(p0,-p2a/b0,atol=2e-9)
    np.testing.assert_allclose(p0,point["weak_overlap_from_four_poles"],rtol=1e-8)
    tmax = 9/(9-8*p0)
    def point_at(t, b_ref=b0):
        p = 1-t*(1-p0)
        real_eps = -(t-1)/(2*math.sqrt(p))
        imag_eps = math.sqrt(max(0.,2*(t-1)-real_eps**2))
        eps = complex(real_eps,imag_eps)
        b = b_ref/t
        m = masses(1.,b,p,eps)
        np.testing.assert_allclose(m,masses(1.,b_ref,p0,0.),atol=2e-13)
        return p,eps,m

    for b_ref in (.3,1.,2.5):
        for t in np.linspace(1,tmax,12):
            p,eps,m = point_at(t,b_ref)
            aa,bb,e,c,w = reconstruct(m,p)
            np.testing.assert_allclose([aa,bb,e,c,w],
                    [1,b_ref/t,abs(eps)**2,math.sqrt(p)*eps.real,abs(math.sqrt(p)-eps)**2],atol=3e-13)
    z = complex(*read("../1723/results.json")["exact_and_conditional_implementations"]["second_complex_root"]["z"])
    # External Vub range is a deliberately generous matching budget, rather than a fit.
    vmax = .01
    floor = math.sqrt(1-p0-vmax**2)
    v_required = math.sqrt(1-p0-abs(z)**2)
    assert floor > abs(z) and v_required > vmax
    vb_example = .004
    requested = 1-abs(z)**2-vb_example**2
    t = (2+requested)/(2+p0)
    assert 1<t<tmax
    p,eps,m = point_at(t)
    np.testing.assert_allclose(abs(math.sqrt(p)-eps)**2,requested)
    return dict(no_direct_term_reconstruction_source_scale_GeV2=a,
                no_direct_term_b_over_a=b0,no_direct_term_p_from_all_four_poles=p0,
                all_isospectral_true_p_range=[p0/(9-8*p0),p0],
                all_isospectral_transfer_probability_range=[p0,tmax*(2+p0)-2],
                two_arm_root_abs=abs(z),physical_Vud_min_if_abs_Vub_at_most_0p01=floor,
                necessary_abs_Vub_for_equal_root_and_true_Vud=v_required,
                example_abs_Vub=vb_example,matching_complement_of_transmitted_amplitude={
                    "true_abs_Vus":math.sqrt(p), "absolute_direct_mixing":abs(eps),
                    "real_and_imag_direct_mixing":[eps.real,eps.imag],
                    "transmitted_abs_amplitude":abs(math.sqrt(p)-eps),
                    "poles_over_bare_first_mass":m.tolist()})


if __name__ == "__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    value=json.loads(json.dumps(run()),parse_float=lambda s:float(f"{float(s):.10g}"))
    file=HERE/"independent_results.json"
    if args.check:
        assert json.loads(file.read_bytes())==value
        print("1727 verified exact inverse-pole and sharp isospectral bounds including a complex direct term.")
    else:
        with file.open("x",encoding="utf8") as stream:
            stream.write(json.dumps(value,indent=2,allow_nan=False)+"\n")
