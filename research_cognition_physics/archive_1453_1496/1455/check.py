"""Two symmetric sources: spectral distances, triplet dominance, and a six-state witness."""
import argparse
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
EV_PER_GEV = 1e9


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8"))


def minimax(up, left):
    ratios = up/left
    center_mix = (max(ratios)+min(ratios))/2
    spread = np.ptp(ratios)/2
    distance = max(abs(up[i]*left[j]-up[j]*left[i])/(left[i]+left[j])
                   for i, j in combinations(range(3), 2))
    center_dirac = max((up-distance)/left)
    assert np.max(abs(up-center_dirac*left)) <= distance+1e-10
    assert math.isclose(max(abs(ratios-center_mix)), spread, abs_tol=1e-10)
    return distance, spread


def uniform_type_i(up, left):
    u, _, t = up
    a, b, d = left
    cross = math.sqrt((u*u/a+t*t/b)/(b+d))
    optimal_modulus = max(cross, t/b)
    q = optimal_modulus
    return max(0.0, max(u*u/a-d*q*q, b*q*q-t*t/b)-2*t*q)


def structural_checks():
    rng = np.random.default_rng(1455)
    unitary = lambda: np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
    up = np.array([.7, 1.3, 2.7])
    left = np.array([.001, .4, 1.])
    dmin, thetamin = minimax(up, left)
    lower = uniform_type_i(up, left)
    for _ in range(48):
        a, b = unitary(), unitary()
        U, L = (a*up)@a.T, (b*left)@b.T
        c = rng.uniform(0, 80)*np.exp(1j*rng.uniform(0, 2*math.pi))
        D = U-c*L
        mixing = D@np.linalg.inv(L)
        correction = mixing@D.T
        assert np.linalg.norm(D, 2) >= dmin-1e-10
        assert np.linalg.norm(mixing, 2) >= thetamin-1e-10
        assert np.linalg.norm(correction, 2) >= lower-1e-9
    vector = np.array([1., 1j, 0.])
    nilpotent = np.outer(vector, vector)
    assert np.linalg.norm(nilpotent, 2)>1 and np.allclose(nilpotent@nilpotent.T, 0)
    assert np.linalg.matrix_rank(np.eye(3)+nilpotent)==3
    return {"complex_trials": 48, "positive_uniform_I_test_C_delta": lower,
            "nilpotent_D_norm": np.linalg.norm(nilpotent, 2)}


def diagonal_witness(up_gev, ds, da, m3_bare_gev):
    # eV throughout: no ill-conditioned diagonalization of a 10^25 eV spectrum.
    # Match three exact positive light roots, with the first/third portals zero.
    up = up_gev*EV_PER_GEV
    m3 = math.sqrt(da/(1-(up[0]/up[2])**2))
    m1 = up[0]*m3/up[2]
    m2 = math.sqrt(m1*m1+ds)
    C, c = m3_bare_gev*EV_PER_GEV/m3, up[2]/m3
    def light_equation(y):
        return y-m2-(up[1]-c*y)**2/(C*y-m2)
    lo, hi = m2, 1.2*m2
    assert light_equation(lo) <= 0 < light_equation(hi)
    for _ in range(90):
        mid = (lo+hi)/2
        if light_equation(mid) < 0:
            lo = mid
        else:
            hi = mid
    l2=(lo+hi)/2
    assert abs(light_equation(l2)) < 1e-15
    return {"light_eV": [m1, m2, m3], "triplet_diagonal_eV": [m1,l2,m3],
            "maximum_R_GeV": m3_bare_gev, "minimum_R_GeV": C*m1/EV_PER_GEV,
            "induced_block_D2_GeV": (up[1]-c*l2)/EV_PER_GEV,
            "relative_largest_light_shift": (l2-m2)/m3,
            "mix_norm": abs(up[1]-c*l2)/(C*l2)}


def compute():
    inp=read("inputs.json")
    quarks = read(inp["weak_scale_inputs"])["yukawa"]["up"]
    osc = read(inp["normal_order_gap_inputs"])["juno_207_days"]
    n=inp["marginal_half_width_multiples"]
    v, dv = inp["matched_vev_GeV"]
    top_min=(v-n*dv)*(quarks[2][0]-n*quarks[2][1])/math.sqrt(2)
    charm_max=(v+n*dv)*(quarks[1][0]+n*quarks[1][1])/math.sqrt(2)
    up_min=(v-n*dv)*(quarks[0][0]-n*quarks[0][1])/math.sqrt(2)
    ds, eds=osc["solar_21"]
    da=osc["normal_atmospheric_31"]
    # Source lists normal ordering -0.025,+0.027 (in 10^-3 eV^2).
    rho=math.sqrt((ds-n*eds)/(da[0]+n*da[2]))
    mass=inp["physical_reference_mass_GeV"]
    q_crit=(charm_max+mass)/(top_min-mass)
    boxes=[]
    for tau in inp["relative_spectral_transfer_to_triplet_L"]:
        q=max(0., (rho-tau)/(1+tau))
        d=max(0., (q*top_min-charm_max)/(1+q))
        r=max(0., (top_min-charm_max/q)/(2*inp["mixing_expansion_budget"]))
        boxes.append({"spectral_transfer":tau, "floor_triplet_second_over_third":q,
                      "full_six_largest_mass_floor_GeV":d, "bare_heaviest_for_mixing_budget_GeV":r})
    return {"two_width_min_max": {"up_min_GeV":up_min, "charm_max_GeV":charm_max, "top_min_GeV":top_min,
                                  "observable_second_over_third_floor":rho},
            "spectral_bounds":boxes,
            "minimum_tau_for_all_six_below_reference_GeV":(rho-q_crit)/(1+q_crit),
            "unified_up_triplet_and_neutral_blocks_witness":diagonal_witness(
                np.array([x[0] for x in quarks])*v/math.sqrt(2), ds, da[0], inp["positive_example_maximum_bare_R_GeV"]),
            "test_certificates":structural_checks()}


def same(a,b):
    if isinstance(a,dict) and isinstance(b,dict):
        return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
    if isinstance(a,list) and isinstance(b,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-11,abs_tol=2e-18)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    data=compute()
    out=HERE/"results.json"
    if args.write:
        with out.open("x",encoding="utf-8") as stream:
            stream.write(json.dumps(data,indent=2,allow_nan=False)+"\n")
    else:
        assert same(data,read("results.json"))
    print("PASS spectral minimax, oscillation windows, and exact neutral roots")
