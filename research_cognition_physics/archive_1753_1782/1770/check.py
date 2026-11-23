"""Identify what three gauge readings fix in a canonical broken SU(5) completion."""
import argparse
import itertools
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())
LOW = runpy.run_path(str(HERE / CFG["gauge_model"]))
LIMITS = json.loads((HERE / CFG["flavour_search_limits"]).read_bytes())
Q = CFG["gauge_reference_GeV"]
TWO_PI = 2 * math.pi
CAP = math.sqrt(CFG["perturbative_cap_squared_over_pi"] * math.pi)
B_EFFECTIVE = LOW["SUSY"] + 4
B_SU5 = 1
TRIPLET = np.array([.4, 0., 1.])
ADJOINT = np.array([0., 2., 3.])
VECTOR = np.array([-10., -6., -4.])


def measured_at_reference(weak, signs=(0, 0, 0, 0), n=0, strong="lattice"):
    e1, e2, es, et = signs
    gs = LOW["electroweak"](n*e1, n*e2)
    mt, mt_err = LOW["CFG"]["top_MS_matching_GeV"]
    mt += n * et * mt_err
    alpha, sigma = LOW["CFG"]["strong_lattice_5f"]["ALPHA_2026"]
    a3_6 = 1/(alpha+n*es*sigma) + math.log(mt/LOW["MZ"])/(3*math.pi)
    if strong == "flavour":
        g3, dg = LOW["LOW"]["g3"]
        a3_6 = 4*math.pi/(g3+n*es*dg)**2
    trial = LOW["gauge"](weak, gs=gs, mt=mt)
    crossing = trial["crossing_GeV"]
    values = np.full(3, trial["inverse_alpha_U"])
    values[2] += a3_6 - trial["a3_SM"]
    return values - B_EFFECTIVE * math.log(Q/crossing) / TWO_PI


def invariants(a, ref=Q):
    log_h = TWO_PI * np.dot([-1, 3, -2], a) * 5/12
    log_k = TWO_PI * np.dot([5, -3, -2], a) / 36
    return ref*math.exp(log_h), ref*math.exp(log_k)


def heavy(a, h, k, x):
    sigma = k*(k/x)**2
    a5_ref = a[0] - (.4*math.log(Q/h)-10*math.log(Q/x))/TWO_PI
    a5_x = a5_ref - B_SU5*math.log(x/Q)/TWO_PI
    coupling = math.sqrt(4*math.pi/a5_x)
    length = x/coupling
    vev = length/5
    lambda_sigma = 2*sigma/length
    higgsino = LOW["SOFT"]
    lambda_h = (h-higgsino)/length
    mu_sigma = lambda_sigma*vev/4
    mu_h = higgsino + 3*lambda_h*vev
    for mu in (Q, x, Q*1.37):
        expected = a-B_EFFECTIVE*math.log(mu/Q)/TWO_PI
        a5 = a5_ref - B_SU5*math.log(mu/Q)/TWO_PI
        observed = a5+(TRIPLET*math.log(mu/h)+ADJOINT*math.log(mu/sigma)
                       +VECTOR*math.log(mu/x))/TWO_PI
        np.testing.assert_allclose(expected, observed, atol=4e-13, rtol=0)
    np.testing.assert_allclose([5*coupling*vev, 10*mu_sigma, mu_h+2*lambda_h*vev],
                               [x, sigma, h], rtol=3e-15)
    return dict(MX_GeV=x, M_sigma_GeV=sigma, g5_at_MX=coupling,
                bare_contact_scale_GeV=length, V_GeV=vev, lambda_sigma=lambda_sigma,
                lambda_H=lambda_h, mu_sigma_GeV=mu_sigma, mu_H_GeV=mu_h)


def root(fn, low, high):
    assert fn(low) < 0 < fn(high)
    for _ in range(70):
        middle = (low+high)/2
        if fn(middle) < 0:
            low = middle
        else:
            high = middle
    return (low+high)/2


def classify(scale):
    status = {}
    for name, row in LIMITS.items():
        lower, upper = row["lower_upper_bracket_on_required_scale_GeV"]
        if scale < lower:
            status[name] = "all orientations fail"
        elif scale < upper:
            status[name] = "some orientations pass"
        else:
            status[name] = "all orientations pass"
    return status


def solve(weak, strong):
    a = measured_at_reference(weak, strong=strong)
    h, k = invariants(a)
    np.testing.assert_allclose(invariants(a-B_EFFECTIVE*math.log(2)/TWO_PI, 2*Q),
                               (h, k), rtol=5e-14)
    def at(x):
        return heavy(a, h, k, x)
    def point_for_length(length):
        return at(root(lambda x: at(x)["bare_contact_scale_GeV"]-length, .05*k, 10*k))
    ends = {name: root(lambda x: CAP-at(x)[name], .05*k, 10*k)
            for name in ["lambda_sigma", "lambda_H"]}
    minimum = at(max(ends.values()))
    assert minimum["g5_at_MX"] < CAP
    length = minimum["bare_contact_scale_GeV"]
    largest = max(row["lower_upper_bracket_on_required_scale_GeV"][1]
                  for row in LIMITS.values())
    smallest = min(row["lower_upper_bracket_on_required_scale_GeV"][0]
                   for row in LIMITS.values())
    rate_bounds = {key: [(v/length)**4 for v in row["lower_upper_bracket_on_required_scale_GeV"]]
                   for key, row in LIMITS.items()}
    unit_samples = [at(x) for x in (max(ends.values()), 2*k, 3*k)]
    example_transfer = CFG["fixed_strength_illustration"]
    target = CFG["selected_effective_scale_GeV"]
    examples = []
    for bare in [math.sqrt(example_transfer)*target, 2*math.sqrt(example_transfer)*target]:
        pt = point_for_length(bare)
        pt["inside_all_three_caps"] = max(pt[key] for key in
                                           ["g5_at_MX", "lambda_H", "lambda_sigma"]) <= CAP
        pt["classification"] = classify(bare/math.sqrt(example_transfer))
        examples.append(pt)
    box=[]
    for signs in itertools.product((-1, 1), repeat=4):
        box.append(invariants(measured_at_reference(weak, signs, CFG["marginal_widths"], strong)))
    h_range=np.amin(box, axis=0)[0], np.amax(box, axis=0)[0]
    k_range=np.amin(box, axis=0)[1], np.amax(box, axis=0)[1]
    conservative=k_range[0]*(2/CAP**3)**(1/3)
    return dict(weak_masses_GeV=weak.tolist(), strong_input=strong, soft_mass_GeV=LOW["SOFT"],
                effective_inverse_gauge_couplings_at_Q=a.tolist(), M_HC_GeV=h,
                cubic_scale_K_GeV=k, cap=CAP, cap_thresholds_MX_GeV=ends,
                minimum_bare_contact=minimum, minimal_contact_rate_over_limit_at_unit_transfer=rate_bounds,
                transfer_up_to_which_all_points_and_orientations_pass=(length/largest)**2,
                transfer_above_which_all_orientations_can_fail_at_minimal_contact=(length/smallest)**2,
                unit_transfer_required_couplings_at_window_edges={
                    "smallest":point_for_length(smallest), "largest":point_for_length(largest)},
                same_running_different_absolute_scales=unit_samples,
                rate_ratio_first_to_last_at_equal_flavour_and_transfer=(
                    unit_samples[-1]["bare_contact_scale_GeV"]/length)**4,
                fixed_transfer_example=examples,
                two_width_box_H_GeV=list(h_range), two_width_box_K_GeV=list(k_range),
                two_width_universal_bare_scale_floor_GeV=conservative,
                two_width_floor_divided_by_central_nucleon_threshold_squared=(conservative/largest)**2)


def calculate():
    weak_mass=LOW["OLD"]["weak_fermions_GeV"][0]
    weak_pair=LOW["NEW_WEAK"]
    archived=LOW["floor_spectrum"](weak_mass,weak_pair)
    rows={"same_flavour_gauge_and_masses":solve(archived,"flavour"),
          "archived_masses_with_lattice_alpha":solve(archived,"lattice"),
          "direct_template_change_only_low_pair":solve(LOW["floor_spectrum"](224.,weak_pair),"lattice")}
    np.testing.assert_allclose(rows["archived_masses_with_lattice_alpha"]["cubic_scale_K_GeV"],
                               rows["direct_template_change_only_low_pair"]["cubic_scale_K_GeV"],
                               rtol=3e-14)
    matter=runpy.run_path(str(HERE/CFG["flavour_model"]))
    refit=json.loads((HERE/CFG["flavour_solutions"]).read_bytes())["refits"]["c1h"]
    weak_inputs=[refit["moduli"][k] for k in ("Vus","Vub","Vcb")]
    bg, _, _=matter["load_mother"](*weak_inputs, .25)
    start=bg["pt"]["flow"]["segments"](.25)[0]
    assert abs(start[0]-LOW["MZ"])<1e-8
    low_couplings=np.sqrt(start[3][0])
    np.testing.assert_allclose(low_couplings,
                               [LOW["LOW"][k][0] for k in ("g1_GUT","g2","g3")],rtol=1e-9)
    inverse_5f=4*math.pi/low_couplings[2]**2-math.log(LOW["CFG"]["top_MS_matching_GeV"][0]/LOW["MZ"])/(3*math.pi)
    lattice, width=LOW["CFG"]["strong_lattice_5f"]["ALPHA_2026"]
    rows["same_flavour_gauge_and_masses"].update(actual_mother_gauge_MZ=low_couplings.tolist(),
              matched_alpha_5f=1/inverse_5f, discrepancy_in_lattice_widths=(1/inverse_5f-lattice)/width)
    return rows


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args=parser.parse_args()
    result=json.loads(json.dumps(calculate()),parse_float=lambda x:float(f"{float(x):.11g}"))
    path=HERE/"results.json"
    if args.check:
        assert json.loads(path.read_bytes())==result
        print("1770 OK: exact multi-scale matching, superpotential, six searches and cap bounds")
    else:
        with path.open("w",encoding="utf8") as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write("\n")
