"""Independent on-shell decay, lifetime and CMS response check for round 1609."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"
VEV = 246.22
HBARC = 1.973269804e-13  # GeV mm
MW, MZ, MH = 80.38, 91.1876, 125.25
MB, MT = 2.8354, 172.5
M = 1100.0


def coefficients():
    data = json.loads((HERE.parent / "1386" / "inputs.json").read_text("utf-8"))
    _, s13, s23, _ = [x[0] for x in data["quark_sines_and_phase"]]
    weights = np.array([s13**2, (1-s13**2)*s23**2, (1-s13**2)*(1-s23**2)])
    assert abs(sum(weights)-1) < 1e-14
    return weights


def phase(x, y):
    return math.sqrt((1-x-y)**2 - 4*x*y)


def vector(m, boson):
    x, y = (m/M)**2, (boson/M)**2
    return phase(x, y) * ((1-x)**2 + y*(1+x) - 2*y*y)


def scalar(m, boson):
    x, y = (m/M)**2, (boson/M)**2
    return phase(x, y) * ((1+x)*(1+x-y) + 4*x)


def widths(s):
    c2 = 1-s*s
    norm = M**3*s*s/(16*math.pi*VEV**2)
    modes = [vector(q, MW) for q in (0.0022, 1.27, MT)] * coefficients()
    return norm*np.array([*modes, c2*vector(MB, MZ)/2, c2*scalar(MB, MH)/2])


def spin_trace(m, boson, typ):
    """Trace in the Dirac representation, with couplings stripped of mixing."""
    pauli = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
             np.diag([1, -1])]
    identity, zero = np.eye(2), np.zeros((2, 2))
    time = np.diag([1, 1, -1, -1])
    gamma = [time, *[np.block([[zero, a], [-a, zero]]) for a in pauli]]
    chirality = 1j*gamma[0]@gamma[1]@gamma[2]@gamma[3]
    left, right = (np.eye(4)-chirality)/2, (np.eye(4)+chirality)/2
    eq = (M*M+m*m-boson*boson)/(2*M)
    q = math.sqrt(eq**2-m*m)
    outgoing = eq*time - q*gamma[3] + m*np.eye(4)
    incoming = M*(time + np.eye(4))
    if typ == "h":
        vertices = [(M*right+m*left)/VEV]
        pol = np.ones((1, 1))
    else:
        vertices = [a@left * boson/VEV*(math.sqrt(2) if typ == "W" else 1)
                    for a in gamma]
        k_lower = np.array([M-eq, 0, 0, q])
        pol = -time.copy().astype(float)
        pol[:] = np.diag([-1., 1., 1., 1.]) + np.outer(k_lower, k_lower)/boson**2
    amp = sum(pol[i,j]*np.trace(outgoing@a@incoming@time@b.conj().T@time).real/2
              for i,a in enumerate(vertices) for j,b in enumerate(vertices))
    return amp*q/(8*math.pi*M*M)


def solve_length(length):
    low, high = 0.0, .1
    target = HBARC/length
    for _ in range(100):
        mid = (low+high)/2
        if widths(mid).sum() < target:
            low = mid
        else:
            high = mid
    return (low+high)/2


def calc():
    norm = M**3/(16*math.pi*VEV**2)
    for m,b,typ,expected in [(MT,MW,"W",norm*vector(MT,MW)),
                             (MB,MZ,"Z",norm*vector(MB,MZ)/2),
                             (MB,MH,"h",norm*scalar(MB,MH)/2)]:
        assert math.isclose(spin_trace(m,b,typ),expected,rel_tol=2e-14)
    s = .035
    c = math.sqrt(1-s*s)
    r = np.eye(5)
    r[2:4,2:4] = [[c,s],[-s,c]]
    projector = r.T@np.diag([1.,1.,1.,0.,0.])@r
    assert np.allclose(projector[2:4,2:4],[[c*c,s*c],[s*c,s*s]]) and not projector[4].any()
    w = widths(s)
    br = w/w.sum()
    ratio = br[2:]/np.array([.50,.25,.25])
    # CMS arXiv:2209.07327 Table 8. These exclusive bins all use 138/fb.
    cms_yields = {1200:[5.,4.4,6.4,10.7,6.1,13.9,37.8,17.5,11.2],
                  1500:[1.37,.88,1.27,2.46,1.44,2.81,7.8,3.8,2.07]}
    return {"M_GeV":M, "weak_left_sine":s, "unchanged_B5_GeV":2100,
            "widths_GeV_uW_cW_tW_bZ_bh":w.tolist(), "branching_fractions":br.tolist(),
            "width_over_mass":sum(w)/M, "ctau_mm":HBARC/sum(w),
            "bottom_ZL_shift":s*s/2,"hbb_over_SM":c*c,
            "sine_for_ctau_mm":{str(d):solve_length(d) for d in (.1,1.,10000.)},
            "inclusive_pair_xsec_reference_fb":22.,
            "prompt_pair_2_known_third_generation_decays_fb_before_selection":22.*sum(br[2:])**2,
            "yield_multipliers_conditional_on_same_response":[min(ratio)**2,max(ratio)**2],
            "CMS_Table8_selected_template_fb":{str(m):sum(y)/138 for m,y in cms_yields.items()},
            "CMS_Table11_singlet_observed_expected_TeV":[1.47,1.37],
            "CMS_Run2_combination_observed_all_third_generation_TeV":1.49}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calc()
    if args.write:
        with OUT.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
    else:
        def same(a,b):
            if isinstance(a,dict):
                return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
            if isinstance(a,list):
                return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
            return math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-25)
        assert same(result,json.loads(OUT.read_text("utf-8")))
    print("1609 independent cross-check passed")
