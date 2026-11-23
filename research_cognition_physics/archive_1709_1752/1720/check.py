"""Match the five colour masses to the four weak transfer masses and the Higgs curvature."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())
TAU = 2*math.pi


def read(key):
    return json.loads((HERE / CFG[key]).read_bytes())


def evaluate():
    gauge = runpy.run_path(str(HERE / CFG["gauge_inversion"]))
    weak = read("gauge_table")
    baseline = read("baseline")["shared_witness"]
    soft, cutoff = baseline["soft_GeV"], CFG["high_flavor_scale_GeV"]
    mz, bare = weak["mu_ref_GeV"], read("color_settings")["mediator_matching_GeV"]
    yb = read("color_model")["tagged_yukawa_values_desc"][0]
    old_c = np.array(baseline["colored_fermions_GeV"])
    old_w = np.array(baseline["weak_fermions_GeV"])
    motif = read("doublet_model")["same_1712_point"]
    ratios = np.array(motif["heavy_masses_GeV"]) / 3e7
    np.testing.assert_allclose(ratios.prod(), math.sqrt(2+abs(complex(*motif["B_over_H"]))**2), rtol=1e-10)
    b, _, offset = gauge["branch"](True, soft, mz)
    color, ew = gauge["COLOR"], gauge["WEAK"]
    ranges = np.array([weak[k] for k in ("g1_GUT","g2","g3")])
    curvature = read("higgs")
    pole = np.array(baseline["neutral_poles_GeV"])
    neutral_w = pole**2 * np.array(baseline["neutral_pole_Yukawa_squared"]) / (8*math.pi**2)
    target = curvature["target_central_and_two_width_abs_upper_GeV2"][0]
    q0 = curvature["q0_GeV"]

    def levels(U, fermion):
        return np.log(U/fermion).sum()*2/3 + np.log(U/np.hypot(fermion,soft)).sum()/3

    def envelope(a):
        length, difference = gauge["invert"](a, (b, 1, offset))
        U = mz*math.exp(length)
        # Old weak poles may be lowered to MZ independently; this maximizes their
        # total integrated loading. All old colour masses are fixed.
        D = difference + levels(U, old_c) - levels(U, np.full(2, mz))
        base_new = 10/3*math.log(U/bare)
        alpha3_before = a[2] - offset[2] - (b[2]*length+levels(U,old_c))/TAU
        # Any four complete new weak pairs, arbitrary splittings/framing, as long
        # as each partner has mass >= cutoff. Exact four level geometry is stronger.
        rmin_any = math.exp((D + base_new - 4*math.log(U/cutoff))/2)
        minpole_g1_upper = U*math.exp(-(D+base_new-2*math.log(r_at_one))/4)
        return U, alpha3_before, difference, rmin_any, minpole_g1_upper

    r_at_one = 1/math.sqrt(1-yb**2)
    a0 = 4*math.pi/ranges[:,0]**2
    U, inverse0, difference, _, _ = envelope(a0)
    np.testing.assert_allclose([U,inverse0], [baseline["U_GeV"],baseline["inverse_alpha_U"]], rtol=5e-9)
    oldc, oldw = levels(U,old_c), levels(U,old_w)
    shift = difference+oldc-oldw
    assert abs(shift)<1e-7

    def model(R):
        # Parameterize by physical mass ratio, avoiding cancellation near g=yb.
        g = yb / math.sqrt((1-1/R)*(1+1/R))
        mc = bare*np.array([1.,1.,R,R,R])
        bc = 2/3*np.log(U/mc).sum()
        weak_unit = U * math.exp(-(bc + shift + np.log(ratios).sum())/4)
        mw = weak_unit*ratios
        bw = np.log(U/mw).sum()
        inv3 = a0-offset-(b*math.log(U/mz)+color*(oldc+bc)+ew*(oldw+bw))/TAU
        np.testing.assert_allclose(inv3, inv3[2], rtol=0, atol=1e-11)
        return {"R":R, "g":g, "relative_g_excess_over_yb":g/yb-1,
                "five_color_poles_GeV":mc.tolist(), "four_weak_poles_GeV":mw.tolist(),
                "common_inverse_coupling":inv3.mean(), "weak_chain_unit_GeV":weak_unit,
                "weak_overlap_from_four_poles":6**4*np.prod((mw/weak_unit)**2)/np.sum((mw/weak_unit)**2)**4-2}

    # At fixed bare M, the entire weak spectrum scales as sqrt(R), exactly.
    g1 = model(r_at_one)
    r_flavor = r_at_one*(cutoff/min(g1["four_weak_poles_GeV"]))**2
    r_upper = cutoff/bare
    interior = model(CFG["illustrative_color_mass_ratio"])
    R, g = interior["R"], interior["g"]
    assert r_flavor<R<r_upper and min(interior["four_weak_poles_GeV"])>cutoff
    # Same zero-local boundary as R1719, with new heavy weak thresholds contributing
    # no H^dagger H term at this order under the paired spectral assumptions.
    K = curvature["sum_k_squared"]
    c_color = -3*bare**2*g**2/(8*math.pi**2)*(1.5+2*K*math.log(R))

    def c_singlet(s):
        return -s**3*np.sum(neutral_w*(np.log((s*pole/q0)**2)-1))

    lo = curvature["one_parameter_neutral_rescaling"]["s_at_maximum_positive_neutral_c"]
    hi = 1.
    assert (c_color+c_singlet(lo)+target)*(c_color+c_singlet(hi)+target)<0
    for _ in range(65):
        mid = (lo+hi)/2
        if c_color+c_singlet(mid)+target>0:
            lo = mid
        else:
            hi = mid
    s = (lo+hi)/2
    assert abs(c_color+c_singlet(s)+target)<1e-9
    corners = [envelope(4*math.pi/(ranges[:,0]+np.array(sign)*2*ranges[:,1])**2)
               for sign in itertools.product((-1,1),repeat=3)]
    c = np.array(corners)
    assert c[:,3].min()>r_at_one and c[:,4].max()<cutoff
    return {"baseline_unification_GeV":U, "old_g1_example_new_matched_spectrum":g1,
            "at_fixed_baseline_R_min_from_lightest_four_weak_poles":r_flavor,
            "at_fixed_baseline_R_max_if_coloured_poles_below_flavour_cutoff":r_upper,
            "central_conditional_combined_point":interior, "central_combined_Higgs_GeV2":{
                "color":c_color,"neutral_rescaling":s,"neutral":c_singlet(s),
                "target_negative_curvature":-target,"neutral_poles_after_rescaling_GeV":(pole*s).tolist()},
            "two_width_gauge_box_allow_old_weak_retuning_upper_bound_new_lightest_g1_GeV":
                float(c[:,4].max()),
            "two_width_gauge_box_necessary_min_R_even_free_weak_mass_ratios":float(c[:,3].min())}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=json.loads(json.dumps(evaluate()),parse_float=lambda s:float(f"{float(s):.11g}"))
    out=HERE/"results.json"
    if args.check:
        assert json.loads(out.read_bytes())==result
        print("1720 verified absolute one-loop unification, eight threshold-box corners and a common curvature point")
    else:
        with out.open("x",encoding="utf8") as fh:
            fh.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
