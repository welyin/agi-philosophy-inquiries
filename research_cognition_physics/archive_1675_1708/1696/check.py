"""Gauge unification with the two physical strong-angle partner thresholds."""
from pathlib import Path
from itertools import product
import argparse
import json
import math

import numpy as np


ROOT = Path(__file__).resolve().parent
TAU = 2*math.pi
SM = np.array([41/10, -19/6, -7])
MSSM = np.array([33/5, 1, -3])
COLOR = np.array([2/5, 0, 1])
WEAK = 1-COLOR
ADJ = np.array([0, 2, 3])
W = np.array([5, -3, -2])


def read(path):
    return json.loads((ROOT/path).read_bytes())


def bounds(values):
    return [min(values), max(values)]


def branch(mssm, soft, mz):
    if not mssm:
        return SM, 2/3, np.zeros(3)
    offset = ADJ/(12*math.pi)+(SM-MSSM)*math.log(soft/mz)/TAU
    return MSSM, 1, offset


def invert(a, theory):
    b, h, offset = theory
    effective = a-offset
    length = TAU*(W @ effective)/(W @ b)
    corrected = TAU*effective-b*length
    split = (corrected[1]-corrected[2])/h
    return length, split


def third(a, theory, split):
    b, h, offset = theory
    effective = a-offset
    slope = (b[2]-b[1])/(b[0]-b[1])
    strength = h*(1-2*slope/5)/TAU
    inv = offset[2]+effective[1]+slope*(effective[0]-effective[1])-strength*split
    return math.sqrt(4*math.pi/inv)


def actual_spectrum(raw, observed):
    scale = raw["cutoff_TeV"]*1e3
    radius = raw["scalar_magnitude_over_cutoff"]
    phase = raw["scalar_angle_over_pi"]*math.pi
    mix = scale*(radius*np.exp(1j*phase))**3
    bare = raw["colored_M_T_over_cutoff"]*scale
    mat = np.array([[0,mix],[-mix,bare]])
    spectrum = np.linalg.svd(mat, compute_uv=False)
    np.testing.assert_allclose(spectrum, np.array(observed)*1e3, atol=2e-7, rtol=0)
    return mat, spectrum


def direct_test(theory, a, length, mz, mD, mL):
    b, h, offset = theory
    upper = mz*math.exp(length)
    assert all(mz<mu<upper for mu in [*mD,*mL])
    nonuniversal = h*(COLOR*sum(math.log(upper/mu) for mu in mD)
                    + WEAK*sum(math.log(upper/mu) for mu in mL))/TAU
    at_u = a-offset-b*length/TAU-nonuniversal
    np.testing.assert_allclose(at_u, at_u[0], atol=1e-11, rtol=0)
    assert at_u[0]>1
    # Independently integrate up through each individual mass step.
    a_up = np.array(a, copy=True)
    if h == 1:
        # Exact scheme conversion and threshold contribution both live at MS.
        # offset was defined by (SM-MSSM)*ln(MS/MZ)/(2pi)+C/(12pi).
        dt = (offset[0]-ADJ[0]/(12*math.pi))*TAU/(SM[0]-MSSM[0])
        step = mz*math.exp(dt)
        events = [(step, MSSM-SM, -ADJ/(12*math.pi))]
    else:
        events = []
    events += [(mass, h*COLOR, 0) for mass in mD]
    events += [(mass, h*WEAK, 0) for mass in mL]
    running = SM.copy()
    scale = mz
    for mass, change, shift in sorted(events,key=lambda item:item[0]):
        a_up = a_up-running*math.log(mass/scale)/TAU+shift
        running += change
        scale = mass
    a_up -= running*math.log(upper/scale)/TAU
    np.testing.assert_allclose(a_up, at_u, atol=1e-11, rtol=0)
    return at_u[0]


def compute():
    inp = read("inputs.json")
    data = read(inp["couplings"])
    mz = data["mu_ref_GeV"]
    soft = inp["ordinary_mssm_decoupling_GeV"]
    texture = read(inp["color_texture"])
    mat, masses = actual_spectrum(texture, read(inp["color_mass_reference"])["conditional_masses_TeV_descending"])
    det = abs(np.linalg.det(mat))
    xmax = math.log(det/mz**2)
    bare = abs(mat[1,1])
    n = inp["electroweak_input_halfwidths"]
    rows = [data[key] for key in ("g1_GUT","g2","g3")]
    central = 4*math.pi/np.square([row[0] for row in rows])
    corners = [4*math.pi/np.square(point) for point in product(*[(v-n*dv,v+n*dv) for v,dv in rows])]
    mssm = branch(True, soft, mz)
    out = dict(color_masses_GeV=masses.tolist(), color_mass_product_GeV2=det,
               max_split_above_MZ=xmax, two_width_input_g3=[rows[2][0]-n*rows[2][1],rows[2][0]+n*rows[2][1]])
    assert W @ COLOR == W @ WEAK == 0
    for flag, name in ((False,"SM_Dirac"),(True,"MSSM_complete")):
        theory = branch(flag, soft, mz)
        b, h, offset = theory
        length, x = invert(central, theory)
        entries = [invert(a,theory) for a in corners]
        splits = [xx for _,xx in entries]
        slope = (b[2]-b[1])/(b[0]-b[1])
        strength = h*(1-.4*slope)/TAU
        limit = 0 if flag else xmax
        proj_budget = strength*(min(splits)-limit)
        model = dict(x_central=x,x_two_width=bounds(splits),
                     U_central_GeV=mz*math.exp(length),
                     U_two_width_GeV=bounds([mz*math.exp(l) for l,_ in entries]),
                     formal_weak_geometric_means_GeV=bounds([math.sqrt(det)*math.exp(-xx/2) for xx in splits]),
                     required_additional_projected_inverse_coupling_for_budget=proj_budget,
                     per_coupling_Linf_matching_floor=proj_budget/(abs(slope)+abs(1-slope)+1),
                     tested_split_upper_bound=limit,
                     holdout_g3_upper_with_budget=max(third(a,theory,limit) for a in corners),
                     equal_products_g3_two_width=bounds([third(a,theory,0) for a in corners]))
        assert model["holdout_g3_upper_with_budget"] < out["two_width_input_g3"][0]
        out[name] = model
    length, x = invert(central, mssm)
    scales = [invert(a,mssm)[1] for a in corners]
    ratio = math.exp(-x/2)
    r_range = bounds([math.exp(-xx/2) for xx in scales])
    c24 = [(1-r)/(3+2*r) for r in (ratio,-ratio)]
    assert all(math.isclose(abs((1-3*z)/(1+2*z)),ratio) for z in c24)
    weak_common = masses*ratio
    unified = direct_test(mssm,central,length,mz,masses,weak_common)
    # Keep both flavor mixing vertices as matched, freely vary the weak bare mass.
    controls = []
    for bare_weak in inp["weak_bare_mass_controls_GeV"]:
        weak = mat.copy()
        weak[1,0] *= ratio
        weak[0,1] *= ratio
        weak[1,1] = bare_weak
        spec = np.linalg.svd(weak,compute_uv=False)
        assert math.isclose(np.prod(spec),det*math.exp(-x),rel_tol=1e-12)
        if spec[-1]>mz:
            direct_test(mssm,central,length,mz,masses,spec)
        controls.append(dict(bare_weak_GeV=bare_weak,spectrum_GeV=spec.tolist(),
                             all_above_MZ=bool(spec[-1]>mz),
                             lighter_than_D=bool(spec[-1]<masses[-1])))
    central_log_step=math.log(2)
    dln=(invert(central,branch(True,soft*math.exp(central_log_step),mz))[1]-x)/central_log_step
    assert math.isclose(dln,19/18,abs_tol=1e-12)
    out["matched_MSSM_spectrum"] = dict(
        common_factor_central=ratio,common_factor_range=r_range,
        adjoint_over_singlet_real_solutions=c24,
        common_factor_weak_masses_GeV=weak_common.tolist(),
        common_factor_light_weak_range_GeV=[r*masses[-1] for r in r_range],
        matched_inverse_unified_coupling=unified,
        dsplit_dln_soft_scale=dln,
        soft_ceiling_from_threshold_domain_range_GeV=bounds([soft*math.exp((xmax-xx)/dln) for xx in scales]),
        unchanged_bare_weak_entry_light_mass_upper_central_GeV=det*math.exp(-x)/bare,
        unchanged_bare_weak_entry_light_mass_upper_over_two_widths_GeV=det*math.exp(-min(scales))/bare,
        bare_weak_spectrum_controls=controls)
    return out


def tidy(obj):
    if isinstance(obj,dict):
        return {key:tidy(value) for key,value in obj.items()}
    if isinstance(obj,list):
        return [tidy(value) for value in obj]
    return float(f"{obj:.11g}") if isinstance(obj,(float,np.floating)) else obj


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    verify = parser.parse_args().check
    result = tidy(compute())
    if verify:
        assert read("results.json") == result
        print("1696: physical determinants, gauge boxes, individual mass steps and competitors verified")
    else:
        (ROOT/"results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
