"""Joint single-Higgs color and singlet curvature at one common MSbar boundary."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


def root(f, lo, hi):
    flo = f(lo)
    assert flo * f(hi) < 0
    for _ in range(85):
        mid = (lo+hi)/2
        if f(mid)*flo > 0:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def color_spectrum(z, k, u):
    """Eight physical squared singular values in units of the common bare mass."""
    masses = []
    for eigen in k:
        w = np.sqrt((1+eigen*np.array([1.,-1.]))/2)
        use = w > 0
        w = w[use]
        block = np.zeros((1+len(w), 1+len(w)))
        block[0, 1:] = math.sqrt(z)*w
        block[1:, 0] = -math.sqrt(u)*w
        block[1:, 1:] = np.diag(np.array([1.,-1.])[use])
        masses.extend(np.linalg.svd(block,compute_uv=False)**2)
    assert len(masses)==8
    return np.array(masses)


def differentiated_trace(k, u, L):
    def potential(z):
        a = color_spectrum(z,k,u)
        a = a[a>0]
        return np.sum(a*a*(np.log(a)+L-1.5))
    base = potential(0.)
    def quotient(z):
        return (potential(z)-base)/z
    # Extrapolate over z and z^2 log(z); f contains three massless eigenvalues at z=0.
    steps = 1e-4 * np.array([1.,.5,.25])
    coords = steps/steps[0]
    values = np.array([quotient(z) for z in steps])
    extrap = np.linalg.solve(np.array([np.ones(3),coords,coords*np.log(coords)]).T,values)[0]
    return extrap


def calculate():
    cfg = read("inputs.json")
    yb, ys, yd = read(cfg["color"])["tagged_yukawa_values_desc"]
    mass = read(cfg["color_inputs"])["mediator_matching_GeV"]
    vev = read(cfg["vev"])["v246_GeV"]/math.sqrt(2)
    k = np.array([-yd/yb,ys/yb,-1.])
    K = k@k
    factor = cfg["color_count"]*mass**2/(8*math.pi**2)
    q0 = mass*math.exp(cfg["common_log_Q_over_M"])
    B = 3*(2*math.log(mass/q0)-1)
    spec = read(cfg["neutral"])["shared_witness"]
    pole = np.array(spec["neutral_poles_GeV"])
    weights = pole**2*np.array(spec["neutral_pole_Yukawa_squared"])/(8*math.pi**2)
    S = weights.sum()
    def neutral(q):
        return -weights@(2*np.log(pole/q)-1)
    def color(g,q):
        u = yb**2/(g*g-yb**2)
        return -factor*g*g*(3*(2*math.log(mass/q)-1)+K*math.log1p(u))
    def scaled(s):
        return -s**3*(-neutral(q0)+2*S*math.log(s))
    def weak(g):
        u = yb**2/(g*g-yb**2)
        w = np.sqrt((1+k[0]*np.array([1.,-1.]))/2)
        mat = np.array([[0, g*vev/mass*w[0], g*vev/mass*w[1]],
                        [-math.sqrt(u)*w[0],1,0],[-math.sqrt(u)*w[1],0,-1]])
        left = np.linalg.svd(mat)[0]
        return np.sum(abs(left[1:,-1])**2)
    inputs = read(cfg["higgs"])["atlas"]["mass_GeV"]
    target = inputs[0]**2/2
    top = (inputs[0]+2*inputs[1])**2/2
    assert abs(S-2*target)<4e-5
    u_opt = root(lambda u:u-math.log1p(u)-B/K,.01,100)
    g_opt = yb*math.sqrt(1+1/u_opt)
    c_min = -color(g_opt,q0)
    fixed_room = top+neutral(q0)
    assert c_min>fixed_room>0
    s_opt = math.exp(-1/3+neutral(q0)/(2*S))
    pos_cap = 2*S*s_opt**3/3
    assert abs(pos_cap-scaled(s_opt))<1e-9
    gmax = root(lambda g:-color(g,q0)-top-pos_cap, g_opt,1)
    # A normalized light right zero vector at zero field is a variational test at v.
    lam_cap = (vev*yd/mass)**2
    odds = (gmax*vev/mass/(1-lam_cap))**2
    pmax = odds/(1+odds)
    trial = cfg["illustrative_g"]
    assert color(trial,q0)+scaled(s_opt)>-target
    s_w = root(lambda s:color(trial,q0)+scaled(s)+target,s_opt,1.)
    assert abs(scaled(s_w)+color(trial,q0)+target)<1e-7
    w_orig = read(cfg["physical_weak"])["d_singlet_fractions"]["1"]
    assert abs(weak(1)-w_orig)<2e-12
    assert weak(gmax) <= pmax*(1+1e-11) and pmax < w_orig/100
    scopes = {}
    for q in [q0,mass,cfg["common_soft_comparison_GeV"],cfg["flavor_cutoff_GeV"]]:
        c_d, c_n = color(1,q), neutral(q)
        local = -target-c_n-c_d
        # Transport the local boundary derived at Q0 by the negative of both loop slopes.
        moved = -target-neutral(q0)-color(1,q0) - (3*factor+S)*2*math.log(q/q0)
        assert abs(moved-local)<5e-8
        scopes[f"{q:.5f}"] = {"c_colored":c_d,"c_neutral":c_n,
                               "required_other_local_for_central_target":local}
    soft_row = read(cfg["neutral_soft"])["frozen_30TeV_universal_soft"]["soft_scale"]
    qsoft = cfg["common_soft_comparison_GeV"]
    assert abs(soft_row["Q_GeV"]-qsoft)<1e-5
    comparisons = []
    for g in [g_opt,trial,1.,3.]:
        u = yb**2/(g*g-yb**2)
        for q in [q0,mass,qsoft]:
            L = 2*math.log(mass/q)
            exact = -2*color(g,q)/(factor*g*g)
            calc = differentiated_trace(k,u,L)
            assert abs(calc-exact)<2e-6
            comparisons.append(abs(calc-exact))
    return {"q0_GeV":q0,"signed_down_compression":k.tolist(),"sum_k_squared":K,
            "one_loop_GeV2_per_log_Q_squared":{"g1_color":3*factor,"frozen_neutral":S},
            "target_central_and_two_width_abs_upper_GeV2":[target,top],
            "fixed_fermion_spectrum_common_boundary":{
                "neutral_coefficient_at_Q0":neutral(q0),
                "max_allowed_abs_color_at_two_widths":fixed_room,
                "sharp_minimum_abs_color_over_all_g":c_min,
                "g_at_sharp_minimum":g_opt,"kappa_at_minimum":math.sqrt(u_opt),
                "minimum_excess_above_two_width_target":c_min-fixed_room},
            "same_g1_physics_at_different_subtractions_GeV2":scopes,
            "supersymmetric_boundary_comparator":{
                "reused_neutral_soft_c_at_Q30TeV":soft_row["neutral_CW_coefficient_GeV2"],
                "original_other_local_at_Q30TeV":soft_row["required_other_and_tree_coefficient_GeV2"],
                "if_new_color_is_only_fermionic_required_other_local":
                    -target-soft_row["neutral_CW_coefficient_GeV2"]-color(1,qsoft)},
            "one_parameter_neutral_rescaling":{
                "s_at_maximum_positive_neutral_c":s_opt,"maximum_positive_neutral_c":pos_cap,
                "all_s_two_width_upper_g":gmax,
                "all_s_two_width_upper_d_singlet_fraction":pmax,
                "g_one_d_singlet_fraction":weak(1),
                "weak_defect_benchmark_over_upper":weak(1)/pmax,
                "witness_g":trial,"witness_s":s_w,"witness_color_c":color(trial,q0),
                "witness_neutral_c":scaled(s_w),"witness_weak_d_fraction":weak(trial),
                "witness_singlet_poles_GeV":(s_w*pole).tolist(),
                "witness_neutral_weak_trace":spec["total_weak_current_weight"]/s_w,
                "witness_leading_neutral_Yukawa_multiplier":math.sqrt(s_w)},
            "maximum_SVD_trace_derivative_absolute_error_in_dimensionless_units":max(comparisons)}


def tidy(x):
    if isinstance(x, dict):
        return {k:tidy(v) for k,v in x.items()}
    if isinstance(x, list):
        return [tidy(v) for v in x]
    return float(f"{x:.11g}")


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    table=tidy(calculate())
    if args.check:
        assert table==read("results.json")
        print("1719 one-loop color curvature, full spectra and neutral budget verified")
    else:
        (HERE/"results.json").write_text(json.dumps(table,indent=2)+"\n",encoding="utf-8")

