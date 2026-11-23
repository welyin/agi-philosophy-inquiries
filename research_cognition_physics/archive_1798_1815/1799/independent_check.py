"""All Takagi phases: three down-flavour transitions and three dispersive mesons."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
PAIRS = {"K": (0, 1), "Bd": (0, 2), "Bs": (1, 2)}
MH = 125.


def read(path):
    return json.loads(path.read_bytes())


def modulus_floor(z):
    return max(0., 2*max(abs(z))-sum(abs(z)))


def modulus_ratio(a, b):
    # Fix the tiny up term, then optimize exactly over the charm phase.
    count = 8192
    eps = np.finfo(float).eps
    phi = np.arange(count)*(2*math.pi/count)
    n = a[2]+a[0]*np.exp(1j*phi)
    d = b[2]+b[0]*np.exp(1j*phi)
    p, q = abs(n)**2+abs(a[1])**2, abs(d)**2+abs(b[1])**2
    x, y = n.conjugate()*a[1], d.conjugate()*b[1]
    quad = q*q-4*abs(y)**2
    lin = p*q-4*(x*y.conjugate()).real
    const = p*p-4*abs(x)**2
    rad = np.sqrt(np.maximum(0.,lin*lin-quad*const))
    roots_lo, roots_hi = const/(lin+rad), (lin+rad)/quad
    nmax, dmax, dmin = sum(abs(a)), sum(abs(b)), modulus_floor(b)
    assert dmin > 0 and min(quad)>0 and min(roots_lo)>0
    derivative = 2*nmax*abs(a[0])/dmin**2 + 2*nmax**2*dmax*abs(b[0])/dmin**4
    slack = derivative*math.pi/count + 50*eps
    return {"squared_coupling_ratio_outer":[float(min(roots_lo)-slack),float(max(roots_hi)+slack)],
            "squared_coupling_ratio_sampled_extrema":[float(min(roots_lo)),float(max(roots_hi))],
            "outer_cover_slack":slack}


def meson_matching(legacy, name):
    had = np.array(legacy["CFG"]["mesons"][name]["matrix_elements"])
    width = legacy["CFG"]["matrix_element_relative_half_window"]
    floors, r_h = [], []
    for alpha, devs in itertools.product(legacy['CFG']['QCD_alpha_s_MZ_interval'],
                                         itertools.product([1-width,1+width],repeat=4)):
        lr, ss, _ = legacy["short"](MH,alpha,had*np.array(devs))
        cfg = dict(legacy["KC"])
        a163 = legacy['KA']['scale_alpha'](cfg,alpha)
        a125 = legacy['KA']['scale_alpha'](dict(cfg,mu_GeV=MH),alpha)
        floors.append((lr-abs(ss)) * (.95*(a163/a125)**(12/23))**2)
        r_h.append(abs(ss)/lr)
    # For M>=163 use the all-scale BMU leading-log weights in round 1796.
    adm = np.array([[6.,-40.],[1/6,16+34/3]])
    eig = np.sort(np.linalg.eigvals(adm))
    s_hi = (1+width)*max(abs(had[2]),abs(had[3])/(6*(eig[1]-eig[0])))
    l_lo = (1-width)*had[1]
    r = max(max(r_h),s_hi/l_lo)
    assert r < 1 and min(floors)>0
    return {"h_projection_GeV3_lower":min(floors),"max_abs_S_over_L_h_and_heavier":r}


def phase_error(r, epsilon):
    error = ((1+r)*epsilon+r*epsilon**2/2)/(1-r)
    return math.asin(min(1.,error))


def run():
    gen = runpy.run_path(str(HERE.parent/"1715/check.py"))
    flow, rep = gen["charged_background"](.1)
    up = np.array(rep["Y_up_at_10PeV"])
    v = flow["low_ckm"](rep["gamma_10PeV_deg"],rep["CKM_angles_at_10PeV_rad"])
    arms = {name: up*v[:,i]*v[:,j] for name,(i,j) in PAIRS.items()}
    ratios = {n:modulus_ratio(arms[n],arms["Bs"]) for n in ("K","Bd")}
    ratios["Bs"] = {"squared_coupling_ratio_outer":[1.,1.]}
    legacy = runpy.run_path(str(HERE.parent/"1795/check.py"))
    kernels = {name:meson_matching(legacy,name) for name in PAIRS}
    opposite = read(HERE.parent/"1581/independent_results.json")["distance_from_real_axis_mod_pi_floor_degrees"]
    lo, hi = 0., .05
    for _ in range(60):
        mid = (lo+hi)/2
        angle = sum(phase_error(row['max_abs_S_over_L_h_and_heavier'],mid) for row in kernels.values())
        if math.degrees(angle) < opposite:
            lo = mid
        else:
            hi = mid
    clocks = read(HERE.parent/"1795/results.json")["observed_norm_working_caps_GeV"]
    tops = read(HERE.parent/"1581/results.json")["all_delta_and_sign_boxes"]
    for name,k in kernels.items():
        coef = k["h_projection_GeV3_lower"]*ratios[name]["squared_coupling_ratio_outer"][0]/MH**2
        scale = 3*MH/(4*math.pi)
        k["single_h_mixing_floor_per_unit_bs_width_GeV_over_GeV"] = coef/scale
        k["tree_bs_width_cap_GeV_if_other_abs_le_measured"] = 2*clocks[name]*scale/coef
        if name != "K":
            box = tops[name.lower()]
            other = box["top_clock_upper_per_ps"]*6.582119569e-13/2+box["other_permitted_abs_over_measured"]*clocks[name]
            k["other_abs_budget_from_1581_inputs_GeV"] = other
            k["conditional_tree_bs_width_cap_GeV_with_top_box"] = (clocks[name]+other)*scale/coef
            k["example_branching_cap_at_4p07MeV_unit_width_matching"] = k["conditional_tree_bs_width_cap_GeV_with_top_box"]/.00407
    # Rephasing every down external line cancels around the oriented flavour triangle.
    rng = np.random.default_rng(1799)
    for _ in range(50):
        f = np.exp(2j*math.pi*rng.random(3))
        angles = rng.normal(size=3)
        factors = []
        for name,(i,j) in {"K":(0,1),"Bs":(1,2),"Bd":(2,0)}.items():
            g = arms[name]@f*np.exp(1j*.28)
            q = np.exp(1j*(angles[j]-angles[i]))
            a,b = q*g,q*g.conjugate()
            amp = a*b+.2*(a*a+b*b)/2
            factors.append(amp/abs(amp))
            if name in ("K","Bd"):
                ratio = abs((arms[name]@f)/(arms["Bs"]@f))**2
                rlo,rhi = ratios[name]["squared_coupling_ratio_outer"]
                assert rlo <= ratio <= rhi
        assert abs(np.prod(factors)-1)<2e-12
    return {"source_abs_lower_10PeV":{n:modulus_floor(a) for n,a in arms.items()},
            "ratios":ratios,"physical_mesons":kernels,"relative_complex_transpose_defect_for_1581_separation":lo,
            "three_phase_upper_deg_at_one_percent_defect": math.degrees(sum(phase_error(k['max_abs_S_over_L_h_and_heavier'],.01) for k in kernels.values())),
            "competing_hermitian_model_minimum_phase_deg":opposite}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--check',action='store_true')
    args = ap.parse_args()
    raw = json.dumps(run(),indent=2,allow_nan=False)+'\n'
    path = HERE/'independent_results.json'
    if args.check:
        assert json.loads(raw)==json.loads(path.read_bytes())
    else:
        path.write_bytes(raw.encode())
    print("1799 independent",hashlib.sha256(raw.encode()).hexdigest()[:12])
