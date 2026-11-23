"""Independent pole-response, flavor-trace and LFV tests for round 1595."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def loop_example(left, right, mixing, fermion, higgs, order):
    x, w = np.polynomial.legendre.leggauss(order)
    t, wt = (1 + x) / 2, w / 2
    y, z = t[:, None] * t, t[:, None] * (1 - t)
    weight = (wt * t)[:, None] * wt
    off = mixing * left * right

    def mass_integral(b):
        squares = np.linalg.eigvalsh([[left**2, b], [b, right**2]])
        gap = squares[1] - squares[0]
        return np.dot(wt, np.log1p(t * gap / ((1 - t) * fermion**2 + t * squares[0]))) / gap

    vals, modes = np.linalg.eigh([[left**2, off], [off, right**2]])
    assert 2 * math.sqrt(vals[0]) > higgs
    source = np.outer(modes[0], modes[1])
    weights = (source + source.T)**2 / 2
    assert np.isclose(weights.sum(), 1)

    def pole(q2):
        response = 0.
        for i in range(2):
            for j in range(2):
                den = (1 - t[:, None]) * fermion**2 + y * vals[i] + z * vals[j] - q2 * y * z
                response += weights[i, j] * (weight / den).sum()
        return response / mass_integral(off)

    eps = 1.e-5
    soft_difference = (math.log((1+eps)*mass_integral(off*(1+eps))) -
                       math.log((1-eps)*mass_integral(off*(1-eps)))) / math.log((1+eps)/(1-eps))
    soft, on = pole(0.), pole(higgs**2)
    assert math.isclose(soft, soft_difference, abs_tol=5.e-9) and on > soft > 1
    return dict(soft=soft, on_shell=on, finite_difference=soft_difference,
                scalar_poles_GeV=np.sqrt(vals).tolist(), triangle_weights=weights.tolist())


def example_flavor(c1, c2):
    # An IO null direction (third column); the charged-lepton columns can be complex.
    x, y = .02231, .55
    n = np.array([math.sqrt(x)*np.exp(.7j), math.sqrt((1-x)*y), math.sqrt((1-x)*(1-y))])
    _, _, vh = np.linalg.svd(n.conj()[None, :], full_matrices=True)
    u, w = vh.conj().T[:, 1:].T
    third = n + (.8+1j)*u + (.3-.8j)*w
    left = np.column_stack([u + .2j*w, w, third])
    masses = np.diag([.000510999, .105658, 1.77686])
    right_adj = np.linalg.solve(left, masses)
    assert np.allclose(left @ right_adj, masses) and np.linalg.norm(n.conj() @ left[:, :2]) < 1.e-12

    def response(slopes):
        yh_v = left @ np.diag(slopes) @ right_adj
        normalized = yh_v @ np.diag(1 / np.diag(masses))
        assert np.allclose(n.conj() @ (yh_v - masses), 0, atol=1.e-13)
        assert np.allclose(normalized.trace(), sum(slopes))
        return normalized

    different = response([c1, c2, 1])
    common = response([c1, c1, 1])
    ratio = n[2].conjugate() / n[0].conjugate()
    assert np.allclose(common[0, 2], (common[0, 0]-c1)*ratio)
    assert abs(different[0, 2]-(different[0, 0]-c1)*ratio)>1.e-5
    return {"unequal_responses_Im_diagonal": np.diag(different).imag.tolist(),
            "unequal_responses_Re_diagonal": np.diag(different).real.tolist(),
            "equal_responses_Im_diagonal": np.diag(common).imag.tolist()}


def compute():
    angles = read(HERE.parent/"1370"/"inputs.json")
    s23 = read(HERE.parent/"1580"/"inputs.json")["theta23_sin2"]
    results = {}
    loops = [loop_example(*args, 125.2, 64) for args in [(300, 500, .15, 800), (350, 260, .6, 610)]]
    results["loops"] = loops
    results["physical_flavor_example"] = example_flavor(*[r["on_shell"] for r in loops])
    mh, v, minimum_tau = 124.8, 246.22, 1.775
    # Two distinct transfer scenarios: combined published rate kappas vs LFV/absolute width.
    trace = {}
    for scenario, tau, mu in [(2, (.96,.06), (1.09,.24)), (3, (.95,.10), (1.08,.26))]:
        lower = 3 - tau[0] - 2*tau[1] - mu[0] - 2*mu[1]
        trace[str(scenario)] = [lower, lower-2*.05]
    results["Re_kappa_e_lower_from_trace_unpadded_padded"] = trace
    lfvcap = .002
    width = .0101
    eps = math.sqrt(8*math.pi*width*lfvcap/(.8*.999*mh))*v/minimum_tau
    disk = {}
    for ordering in ["NO", "IO"]:
        a, b = angles[ordering]["s12"][1], angles[ordering]["s13"][2]
        if ordering == "NO":
            c = s23[ordering][1]
            ratio = math.sqrt((1-a)*(1-b))/(math.sqrt(a*c)-math.sqrt((1-a)*(1-c)*b))
        else:
            c = s23[ordering][2]
            ratio = math.sqrt(b/((1-b)*(1-c)))
        disk[ordering] = dict(e_over_tau_max=ratio, radius=eps*ratio,
                              Re_kappa_e_min_using_c_above_one=1-eps*ratio)
    results["common_spectrum_LFV_disk"] = disk
    # Independent maximum of the author's weighted null-direction test.
    tau_floor = 2.94 - 2*math.hypot(.21,.32)
    lfv_caps = np.sqrt(np.array([.002,.0018])*53.6/(.8*.995*tau_floor))
    for ordering, expected in [("NO",1.3074935393), ("IO",.3403541933)]:
        a, b = angles[ordering]["s12"][1], angles[ordering]["s13"][2]
        c = s23[ordering][1 if ordering == "NO" else 2]
        if ordering == "NO":
            col = [math.sqrt((1-a)*(1-b)),
                   math.sqrt(a*(1-c))+math.sqrt((1-a)*b*c),
                   math.sqrt(a*c)-math.sqrt((1-a)*b*(1-c))]
        else:
            col = [math.sqrt(b),math.sqrt((1-b)*c),math.sqrt((1-b)*(1-c))]
        assert np.isclose(np.dot(col,col),1)
        assert math.isclose(np.dot(lfv_caps,col[:2])/col[2],expected,abs_tol=1e-9)
    results["B_ee_lower_if_10p1_MeV_and_trace_scenario2_padded"] = (
        mh / (8*math.pi*width) * (.00051099/v)**2 * trace["2"][1]**2
        * (1-(2*.000511/mh)**2)**1.5)
    return results


def equal(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=2e-8,abs_tol=1e-10)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = compute()
    path = HERE/"independent_results.json"
    if args.write_new:
        with path.open("x", encoding="utf-8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2,allow_nan=False)
            stream.write("\n")
    elif not equal(result,read(path)):
        raise SystemExit("stored independent result differs")
    print(json.dumps(result, ensure_ascii=False))
