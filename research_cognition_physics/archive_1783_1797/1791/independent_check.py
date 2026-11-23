"""Two on-shell Higgs poles from a fermion scalar current at one inert-scalar loop."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def bubble(s, pole):
    """B0(s)-B0(0), and derivative, on the first real sheet below threshold."""
    a = math.sqrt(4 * pole**2 / s - 1)
    angle = math.atan(1 / a)
    return 2 - 2*a*angle, ((1+a*a)*angle/a - 1)/s


def triangle_real(s, h2, mass, nodes, weights):
    """Integrate exactly in x-y first; x+y is the remaining parameter."""
    a = mass**2 - s*nodes**2/4 - h2*nodes*(1-nodes)
    if np.min(a) <= 0:
        raise ValueError("first-sheet real triangle check needs energies below the cut")
    value = 2/np.sqrt(s*a)*np.arctan(nodes*np.sqrt(s/a)/2)
    return np.dot(weights, value)


def triangle_im_over_pi_beta(s, h2, mass, nodes, weights):
    """Explicit Feynman-parameter delta function for the s-channel cut."""
    c = s/4-h2
    root = 2*mass**2 / (math.sqrt(h2**2+4*c*mass**2) + h2)
    # Use the exact factored difference 1-root to avoid cancellation at threshold.
    one_minus_root = (s/4-mass**2)/(h2+c*(1+root))
    w = root + one_minus_root * nodes**2
    regular = 2*math.sqrt(one_minus_root/s) / np.sqrt(h2+c*(w+root))
    beta = math.sqrt((s-4*mass**2)/s)
    return float(np.dot(weights, regular)/beta)


def run(order):
    data = json.loads((HERE / "inputs.json").read_bytes())
    v, mh, m0 = data["vev_GeV"], data["higgs_mass_GeV"], data["charged_mass_GeV"]
    h2, g3 = mh**2, 3*mh**2/v
    masses = np.array([[m0, data["odd_neutral_mass_GeV"], m] for m in data["even_neutral_mass_GeV"]]).ravel()
    count = np.tile([2, 1, 1], 2)
    sign = np.repeat([1, -1], 3)
    bp = np.array([bubble(h2, m)[1] for m in masses])
    offset = 2*(masses**2-m0**2)/v**2
    axis = -np.dot(count*bp*sign, offset)/np.sum(count*bp)
    points = axis + np.array([-1,1])*data["portal_half_separation"]
    couplings = v*(offset[None,:] + points[:,None]*sign[None,:])
    wave = np.sum(count*couplings**2*bp, axis=1)/(32*math.pi**2)
    assert abs(wave[0]-wave[1]) < 2e-14
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = (nodes+1)/2, weights/2
    rows = []
    for energy in (mh, 250.0, 300.0, 350.0, 400.0, 430.0):
        s = energy**2
        triangle = np.array([triangle_real(s,h2,m,nodes,weights) for m in masses])
        bdiff = np.array([bubble(s,m)[0]-bubble(h2,m)[0] for m in masses])
        quotient = bdiff/(s-h2) if s != h2 else bp
        tr = np.sum(count*couplings**3*triangle, axis=1)/(16*math.pi**2*g3)
        seagull = -np.sum(count*couplings**2*bdiff, axis=1)/(32*math.pi**2*v*g3)
        propagator = -np.sum(count*couplings**2*quotient, axis=1)/(32*math.pi**2)
        physical = tr + seagull + propagator - wave
        rows.append({"sqrt_s_GeV": energy, "signed_delta_triangle": tr[0]-tr[1],
                     "signed_delta_bubble_vertex": seagull[0]-seagull[1],
                     "signed_delta_internal_line": propagator[0]-propagator[1],
                     "net_branch_difference_per_tree_s_channel": physical[0]-physical[1]})

    idx = 2  # 255 GeV neutral. Its branch point is distinct from the common 218,220 GeV pairs.
    m = masses[idx]
    cutrows = []
    for shift in (1e-3,1e-5,1e-7):
        s = 4*m**2*(1+shift)
        beta = math.sqrt((s-4*m**2)/s)
        k, q = s/2-h2, s/2*beta*math.sqrt(1-4*h2/s)
        angular = math.atanh(q/k)/q
        analytic_cut = triangle_im_over_pi_beta(s,h2,m,nodes,weights)
        assert abs(angular/analytic_cut-1) < 2e-9
        cutrows.append({"s_over_pair_threshold_minus_one": shift, "C0_Im_over_pi_beta": analytic_cut,
                        "tree_pair_angular_kernel": angular})
    g = couplings[:,idx]
    s = 4*m**2
    threshold_tree = -g/v - g*g3/(s-h2) + 2*g**2/(2*m**2-h2)
    singular_coeff = g*threshold_tree/(32*math.pi*g3)
    assert abs((2*g**3/(2*m**2-h2)-g**2/v-g**2*g3/(s-h2))[0]/
               (32*math.pi*g3)-singular_coeff[0]) < 1e-12
    # XENONnT 3.1 tonne-year ancillary SI csv, adjacent upper_limit rows (cm^2).
    anchors = [[215.44346900318843,8.184351415813917e-47],[225.7019719633921,8.562027553913661e-47]]
    target = data["odd_neutral_mass_GeV"]
    lp = math.log(target/anchors[0][0])/math.log(anchors[1][0]/anchors[0][0])
    recoil_limit = anchors[0][1]*(anchors[1][1]/anchors[0][1])**lp
    nuc, fn = .939, .308-2*.018
    factor = .3893793721e-27*fn**2*nuc**4/(4*math.pi*mh**4*(target+nuc)**2)
    cross = factor*(couplings[:,[1,4]]/v)**2
    recoil = {"data_url":"https://arxiv.org/src/2502.18005v2/anc/xenonnt_2025_si_wimp.csv",
              "anchors_mass_GeV_limit_cm2":anchors, "log_interpolated_90pct_cm2":recoil_limit,
              "working_lower_higgs_nucleon_element":fn, "nucleon_mass_GeV":nuc,
              "sigma_cm2_for_two_stable_states_in_each_branch":cross.tolist(),
              "necessary_total_local_fraction_upper_allowing_arbitrary_split":(recoil_limit/cross.min(axis=1)).tolist(),
              "sufficient_total_local_fraction_upper_any_split":(recoil_limit/cross.max(axis=1)).tolist()}
    return {"external_current": "fermion scalar current, two Higgs poles, leading inert loop, common tree Yukawa and other input match",
            "sign_convention": "branch x_minus minus branch x_plus", "portal_points": points.tolist(),
            "wave_W": wave.tolist(), "below_inert_threshold": rows,
            "isolated_neutral_pair_cut_energy_GeV": 2*m,
            "cut_onsets": cutrows, "tree_NN_to_hh_at_threshold": threshold_tree.tolist(),
            "scalar_current_Im_over_beta_per_tree_pole_channel": singular_coeff.tolist(),
            "signed_cut_coefficient_difference": singular_coeff[0]-singular_coeff[1],
            "same_stable_odd_particles_nuclear_recoil":recoil}


def portable(obj):
    if isinstance(obj, dict):
        return {k: portable(v) for k,v in obj.items()}
    if isinstance(obj, list):
        return [portable(v) for v in obj]
    if isinstance(obj, np.integer):
        return int(obj)
    return float(f"{obj:.9g}") if isinstance(obj, (float,np.floating)) else obj


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--check", action="store_true")
    args = cli.parse_args()
    answer = portable(run(80))
    assert answer == portable(run(120))
    text = json.dumps(answer, indent=2, allow_nan=False) + "\n"
    destination = HERE/"independent_results.json"
    if args.check:
        assert json.loads(destination.read_bytes()) == answer
    else:
        destination.write_text(text, encoding="utf-8")
    print("1791 scalar-current branch", hashlib.sha256(text.encode()).hexdigest()[:12])
