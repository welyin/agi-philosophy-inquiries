"""Three real Hermitian sources shared by two 3-state mass sectors."""
import argparse
import json
from itertools import product
from math import asin, cos, pi, sin, sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
I = np.eye(3)


def load(path):
    return json.loads((HERE / path).read_text(encoding="utf-8"))


def unit(h):
    h = h - I*np.trace(h)/3
    return h / np.linalg.norm(h)


def pair(h, k):
    a = unit(h)
    k = unit(k)
    b = k - np.vdot(a, k).real * a
    return a, b / np.linalg.norm(b)


def root_frame(plane):
    a, b = plane

    def direction(theta):
        return cos(theta)*a + sin(theta)*b

    left, right = 0., pi
    sign = np.linalg.det(a).real
    if abs(sign) < 1.e-14:
        theta = left
    else:
        for _ in range(62):
            mid = (left+right)/2
            if np.linalg.det(direction(mid)).real * sign > 0:
                left = mid
            else:
                right = mid
        theta = (left+right)/2
    h = direction(theta)
    eig, basis = np.linalg.eigh(h)
    assert np.allclose(eig, [-1/sqrt(2), 0, 1/sqrt(2)], atol=1.e-12)
    return h, basis, -sin(theta)*a + cos(theta)*b


def mix(s12, s23, s13, delta):
    def rotation(i, j, s, phase=0.):
        r = np.eye(3, dtype=complex)
        r[i, i] = r[j, j] = sqrt(1-s*s)
        r[i, j] = s*np.exp(-1j*phase)
        r[j, i] = -np.conj(r[i, j])
        return r
    return rotation(1, 2, s23) @ rotation(0, 2, s13, delta) @ rotation(0, 1, s12)


def real_vector(h):
    return np.r_[h.real.ravel(), h.imag.ravel()]


def group_test(qpair, lpair, pmns):
    q_shared, q_frame, q_other = root_frame(qpair)
    l_shared, l_frame, l_other = root_frame(lpair)
    w = q_frame @ l_frame.conj().T
    assert np.linalg.norm(w @ l_shared @ w.conj().T - q_shared) < 1.e-12
    sources = np.stack([real_vector(s) for s in (q_shared, q_other, w @ l_other @ w.conj().T)], axis=1)
    objects = [*qpair, *(w @ s @ w.conj().T for s in lpair)]
    errors = []
    for m in objects:
        target = real_vector(m)
        coeff = np.linalg.lstsq(sources, target, rcond=None)[0]
        errors.append(np.linalg.norm(sources @ coeff-target))
    assert max(errors) < 2.e-12
    # Any diagonal phase in the common-source basis preserves the calibration.
    terms = q_frame[0, :] * l_frame[0, :].conj()
    lengths = np.abs(terms)
    cross_interval = [max(0., 2*max(lengths)-sum(lengths)), sum(lengths)]
    aligned = q_frame @ np.diag(np.exp(-1j*np.angle(terms))) @ l_frame.conj().T
    assert abs(abs(aligned[0,0])-cross_interval[1]) < 1.e-12
    j = (pmns[0,0]*pmns[1,1]*pmns[0,1].conjugate()*pmns[1,0].conjugate()).imag
    return {"rank":int(np.linalg.matrix_rank(sources, tol=1.e-9)), "lepton_J":float(j),
            "one_cross_reference_overlap":list(map(float,cross_interval)),
            "all_objects_reconstructed":bool(max(errors) < 2.e-12)}


def gap_box(marginals):
    ratios = []
    for signs in product((-2,2), repeat=3):
        x,y,z = [(value+sign*width)**2 for (value,width),sign in zip(marginals,signs)]
        ratios.append((y-x)/(z-x))
    return [min(ratios), max(ratios)]


def main():
    info = load("inputs.json")
    light, uc, heavy = map(load, info["charged_yukawas"])
    light_mz = light["points"]["MZ"]
    light_inputs = [(v*u,e*u) for v,e,u in zip(light_mz["central"], light_mz["half_width"], light["units"])]
    up = [uc["MZ"]["up"], uc["MZ"]["charm"], heavy["yt"]]
    down = [*light_inputs[:2], heavy["yb"]]
    lepton = [*light_inputs[2:], heavy["ytau"]]

    def spectrum(values):
        return np.diag([v**2 for v,_ in values])
    v = mix(*info["quark_S12_S23_S13_delta"])
    qpair = pair(spectrum(up), v @ spectrum(down) @ v.conj().T)
    a = load(info["reactor"])
    s13 = sin(0.5*asin(sqrt(a["daya_bay"]["sin_sq_2theta13"])))
    s12 = sqrt(a["juno_primary"]["sin_sq_theta12"])
    s23 = sqrt(info["T2K_normal_s23_squared"])
    delta21 = info["delta21_eV2"]
    nu = np.diag([0., delta21, delta21+info["T2K_normal_delta32_eV2"]])
    examples = []
    for delta in info["phase_examples_radian"]:
        u = mix(s12,s23,s13,delta)
        lpair = pair(spectrum(lepton), u @ nu @ u.conj().T)
        examples.append(group_test(qpair,lpair,u))
    # Four distinct fixed dictionary directions, with and without a freely chosen identification.
    def off(i,j,value):
        h = np.zeros((3,3),complex)
        h[i,j] = value
        h[j,i] = value.conjugate()
        return h
    modes = [off(0,1,1),off(0,1,1j),off(0,2,1),off(0,2,1j)]
    fixed_rank = np.linalg.matrix_rank(np.stack([real_vector(h) for h in modes]))
    assert fixed_rank == 4
    def bracket(h,k):
        a,b=pair(h,k)
        return np.linalg.norm(a@b-b@a)**2
    bq, bl = bracket(modes[0],modes[1]), bracket(modes[0],modes[2])
    assert np.allclose([bq,bl],[2,.5])
    up_box, lepton_box = gap_box(up), gap_box(lepton)
    assert up_box[1] < lepton_box[0] < lepton_box[1] < .5
    return {"examples":examples, "fixed_dictionary_rank":int(fixed_rank),
            "analytic_sharpness_bracket_pair":[float(bq),float(bl)],
            "affinely_named_up_lepton_ordered_gap_two_width_boxes":{"up":up_box,"lepton":lepton_box}}


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check", action="store_true")
    options = cli.parse_args()
    result = main()
    if options.check:
        def compare(x,y):
            if isinstance(x,dict):
                assert x.keys()==y.keys()
                for k in x:
                    compare(x[k],y[k])
            elif isinstance(x,list):
                assert len(x)==len(y)
                for m,n in zip(x,y):
                    compare(m,n)
            elif isinstance(x,float):
                assert abs(x-y)<2.e-9
            else:
                assert x==y
        compare(result, load("results.json"))
        print("1388 common-source synthesis, dictionary and named-source checks passed.")
    else:
        print(json.dumps(result,indent=2))
