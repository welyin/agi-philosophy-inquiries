"""Inclusive R2 rates, pure-state comparison and the shared-mass normalization."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8-sig"))


def mixing(s12sq, s13sq, s23sq, delta):
    s, r, t = np.sqrt([s12sq, s13sq, s23sq])
    c, v, w = np.sqrt(1-np.array([s12sq, s13sq, s23sq]))
    z = np.exp(1j*delta)
    return np.array([[c*v,s*v,r*z.conjugate()],
                     [-s*w-c*t*r*z,c*w-s*t*r*z,t*v],
                     [s*t-c*w*r*z,-c*t-s*w*r*z,w*v]])


def phase_space(M, mq, ml):
    x, y = (mq/M)**2, (ml/M)**2
    return (1-x-y)*math.sqrt((1-x-y)**2-4*x*y)


def source_witness(prior, inputs, neutrinos):
    scenario = prior["rows"]["u"]["IO_joint_witness"]
    p = scenario["p_e_mu_tau"]
    example = inputs["comparison_vacuum"]
    U = mixing(example["sin2theta12"],scenario["sin2_theta13"],
               scenario["sin2_theta23"],math.radians(example["delta_deg"]))
    # The antineutrino weak mixing matrix is conjugated. The sign flip is physical
    # only when all other preparation/detection vertices are held fixed.
    mixing_bar = U.conj()
    n = mixing_bar[:,2].copy()
    other = n.copy()
    other[2] *= -1
    assert np.allclose(np.abs(other)**2,p) and np.allclose(np.abs(n)**2,p)
    osc = neutrinos["IO"]
    dm = [osc["atm_abs"][0]-osc["solar"][0],osc["atm_abs"][0],0.]
    phases = np.exp(-1j*example["phase_coefficient"]*example["L_km_per_E_GeV"]*np.array(dm))
    evolution = (mixing_bar*phases)@mixing_bar.conj().T
    assert np.linalg.norm(evolution @ n-n) < 1.e-12
    fake = np.abs(evolution @ other)**2
    assert np.linalg.norm(fake-p)>0.05

    # independent_check defines the reported singular value as ||E-kD||.
    delta_norm = scenario["rank1_correction_singular_value"]
    k = scenario["k"]
    bare_B = delta_norm/(4*k)
    unit_s_coeff = 2*math.sqrt(3)*k*bare_B
    pairs = []
    for sH in inputs["EW_doublet_projection_witness"]:
        yukawa = unit_s_coeff / sH
        mass_ratio = (inputs["masses_GeV"]["pole_benchmark"]/1000) / yukawa
        pairs.append({"sH":sH,"Y_R2_at_unit_transfer":yukawa,
                      "M_TeV_over_Y":mass_ratio})
    assert pairs[0]["M_TeV_over_Y"] < inputs["meg_IO_min_TeV_at_unit_norm"]
    assert pairs[1]["M_TeV_over_Y"] > inputs["meg_IO_sufficient_TeV_conservative"]
    return {"p_e_mu_tau":p, "alternative_massless_neutrino_fraction":abs(n.conj()@other)**2,
            "alternative_vacuum_detected_anti_e_mu_tau":fake.tolist(),
            "k":k,"B_norm":bare_B,"Higgs_projection_examples":pairs}


def compute():
    ins = read("inputs.json")
    prior = read("../1601/independent_results.json")
    main = read("../1601/results.json")["all_phases_and_closed_two_width_box"]
    charged = np.array(read(ins["mass_and_CKM"][0])["scales"]["1e16"])
    ckm = np.array(read(ins["mass_and_CKM"][1])["high_scale"]["angles_theta12_theta13_theta23_delta"])
    neutrinos = read(ins["oscillations"])
    low_atm = read(ins["atmospheric"])["theta23_sin2"]["IO"][1]
    e_low, e_high = neutrinos["IO"]["s13"][1:]
    cap = max(v["IO_sin2theta23_upper"] for v in prior["rows"].values())
    assert low_atm < cap < 0.5
    bounds = {}
    for label, result in prior["rows"].items():
        high_atm = result["IO_sin2theta23_upper"]
        mixed = [2*(1-e_high)**2*low_atm*(1-low_atm),
                 2*(1-e_low)**2*high_atm*(1-high_atm)]
        bounds[label] = {"p_mu": [(1-e_high)*low_atm,(1-e_low)*high_atm],
                         "p_tau": [(1-e_high)*(1-high_atm),(1-e_low)*(1-low_atm)],
                         "mixed_mu_tau_pair":mixed}
    nominal_low = min(v["mixed_mu_tau_pair"][0] for v in bounds.values())
    epsilon = ins["charged_flavor_state_relative_amplitude_error"]
    # Trace distance of two independent records is at most twice the single bound.
    stressed = nominal_low - 2*epsilon
    efficiency_ratio = 0.5
    selected = efficiency_ratio*stressed/(1-stressed+efficiency_ratio*stressed)
    trials = math.ceil(math.log(0.05)/math.log1p(-selected))

    lo, hi = charged[:,0]-2*charged[:,1], charged[:,0]+2*charged[:,1]
    kmin = main["general_IO_k_min_zero_transfer"]
    bottom_cap = main["bottom_probability_max_on_tag_all"]
    s2low = lo[2]*math.sqrt(1-bottom_cap)
    angle_lo, angle_hi = ckm[:,0]-2*ckm[:,1],ckm[:,0]+2*ckm[:,1]
    s12, _, s23 = np.sin(angle_lo[:3])
    c12, c13, c23 = np.cos(angle_hi[:3])
    _, s13max, s23max = np.sin(angle_hi[:3])
    vjd_min = min(c12*c13,s12*c23-s23max*s13max,s12*s23-s13max)
    Nmax = hi[0]/vjd_min
    tau_low = (1-e_high)*(1-cap)
    Bmin = (s2low*math.sqrt(tau_low)-Nmax)/4
    y_high_min = 2*math.sqrt(3)*kmin*Bmin
    amplitude_transfer = ins["tree_transfer_min"]
    y_low_min = y_high_min*amplitude_transfer
    m = ins["masses_GeV"]
    widths = {}
    for M in (m["pole_fermionic_lower_reference"],m["pole_benchmark"]):
        F = phase_space(M,m["top_upper"],m["tau_upper"])
        width_min = M*y_low_min**2*F/(16*math.pi)
        widths[str(M)] = {"phase_space_floor":F,"partial_width_min_GeV":width_min,
                           "proper_length_max_mm":ins["hbarc_GeV_mm"]/width_min}
    Bbar = ins["meg_IO_min_TeV_at_unit_norm"]
    M_TeV = m["pole_benchmark"]/1000
    keep = math.sqrt(1-ins["meg_same_chirality_relative_amplitude_budget"])
    return {"IO_electron": [e_low,e_high],"IO_flavor_bounds":bounds,
            "mixed_pair_after_amplitude_budget":stressed,
            "conditional_selected_mixed_probability":selected,
            "ideal_iid_background_free_number_for_zero_p005":trials,
            "s2low":s2low,"Nmax":Nmax,"B_singular_lower":Bmin,
            "sH_times_high_yukawa_lower":y_high_min,
            "sH_times_low_yukawa_lower":y_low_min,
            "one_TeV_sH_lower_MEG_without_cancellation":Bbar*y_low_min/M_TeV,
            "one_TeV_sH_lower_with_amplitude_budget":keep*Bbar*y_low_min/M_TeV,
            "conditional_width_lifetime":widths,
            "same_mass_competitors":source_witness(prior,ins,neutrinos)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    check = parser.parse_args().check
    rendered = json.dumps(compute(),indent=2,allow_nan=False)+"\n"
    path = ROOT/"results.json"
    if check:
        assert path.read_text(encoding="utf-8-sig")==rendered
        print("1604: charged pair, neutral preparation, width and MEG projection verified")
    else:
        path.write_text(rendered,encoding="utf-8")


if __name__ == "__main__":
    main()
