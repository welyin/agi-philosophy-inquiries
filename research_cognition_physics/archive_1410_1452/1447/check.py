"""Check the covariant closure and the aligned Higgs-rate holdout."""
import argparse
import json
from math import cos, sin, pi
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent

def read(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8-sig"))

def cubic(y):
    return y @ y.conj().T @ y

def t_from_k(k):
    return (k-1)/(3-k)

def k_from_t(t):
    return (1+3*t)/(1+t)

def invert_g(value, lo, hi, increasing):
    def g(t):
        return t*(1+t)**2
    assert min(g(lo),g(hi)) <= value <= max(g(lo),g(hi))
    for _ in range(80):
        mid = (lo+hi)/2
        if (g(mid)<value) == increasing:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2

def run():
    data = read("inputs.json")
    branches = read(data["hs_branches"])["branches"]
    competition = []
    rng = np.random.default_rng(1447)
    def unitary():
        return np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
    left, right = unitary(), unitary()
    for row in branches:
        l0,b,A,cu,cd = row["parameters_lambda0_b_A0_cu_cd"]
        values = []
        for z,c in [(-1j*l0*sin(pi/8),cu), (l0*cos(pi/8),cd)]:
            y = np.array([[c*abs(z)**4, b*z**3, 0],
                          [b*z.conjugate()**3, b*abs(z)**2, A*z**2],
                          [0,A*z.conjugate()**2,1]])
            product = cubic(y)
            assert np.allclose(cubic(left@y@right.conj().T),left@product@right.conj().T)
            assert np.isclose(product[0,2], b*A*z**5*np.trace(y),atol=1e-18)
            assert abs(product[0,2]) > 1e-10
            _, before = np.linalg.eigh(y)
            _, after = np.linalg.eigh(y+.1*product)
            assert np.allclose(abs(before.conj().T@after), np.eye(3), atol=1e-9)
            values.append(round(float(abs(product[0,2])),12))
        competition.append({"branch":row["signs_middle_down_and_light_up_over_charm_down_over_strange"],
                            "abs_cubic_13_up_down":values})
    bottom = read(data["direct_ratio"])["direct_ATLAS_ratio_boxes"]["bb_ZZ"]["2"]["allowed_abs_kf_over_kZ"]
    _, strange, beauty = read(data["mass_calibration"])["yukawa"]["down"]
    ratio = (strange[0]+2*strange[1])/(beauty[0]-2*beauty[1])*(1+data["max_relative_transport_mZ_to_h"])
    strange_k = []
    for k in bottom:
        t = t_from_k(k)
        target = ratio**2*t*(1+t)**2
        t_s = invert_g(target,-1/3,2.,True)
        strange_k.append(k_from_t(t_s))
    t_b = data["wide_competitor_t_bottom"]
    test_ratio = data["wide_competitor_mass_ratio"]
    same_target = test_ratio**2*t_b*(1+t_b)**2
    low = invert_g(same_target,-1/3,0,True)
    folded = invert_g(same_target,-1,-1/3,False)
    return {"seed_cubic_13":competition,"absolute_bottom_kappa":bottom,
            "same_scale_strange_bottom_mass_ratio_max":ratio,"normal_branch_strange_kappa":strange_k,
            "same_bottom_kappa_competition":k_from_t(t_b),
            "same_mass_ratio":test_ratio,"strange_kappa_normal_folded":[k_from_t(low), k_from_t(folded)]}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    options = parser.add_mutually_exclusive_group()
    options.add_argument("--write-new",action="store_true")
    options.add_argument("--check",action="store_true")
    args = parser.parse_args()
    results = run()
    path = ROOT/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf8") as out:
            json.dump(results,out,indent=2)
            out.write("\n")
    elif args.check:
        assert results == read("results.json")
        print("1447 flavour closure and aligned two-Higgs-rate contract verified")
    else:
        print(json.dumps(results,indent=2))
