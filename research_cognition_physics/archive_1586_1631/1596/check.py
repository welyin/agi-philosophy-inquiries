"""From independent charged dipoles to the light-neutrino Majorana mass in xenon."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8"))


def br_cap(light, heavy, life, bound, q, dip):
    return math.sqrt(4*dip["hbar_GeV_s"]*bound /
                     (life*dip["alpha"]*heavy**3*(1-(light/heavy)**2)**3*q))


def least_component(row_sq, fixed, coeff_unknown, coeff_fixed):
    # row_sq <= unknown^2 + fixed^2 + (a*unknown+b*fixed)^2
    c = row_sq-(1+coeff_fixed**2)*fixed**2
    if c <= 0:
        return 0.
    linear = coeff_unknown*coeff_fixed*fixed
    return c / (math.sqrt(linear**2+(1+coeff_unknown**2)*c)+linear)


def calculate():
    inputs = read("inputs.json")
    neutrino = read(inputs["neutrino"])
    io = neutrino["IO"]
    dip, ref = read(inputs["electromagnetic"]), read(inputs["muon"])
    shared = read(inputs["common"])
    emu = read(inputs["recoil_results"])["recoil_ae_rectangles_one_two"]["Rb"]
    mu_mean, mu_width = ref["WP25_difference_as_published"]
    mu_upper = (mu_mean+2*mu_width)*1e-11
    me, mm = dip["me_GeV"], dip["mmu_GeV"]
    mt = shared["m_tau_GeV"][0]
    s12lo, s12hi = io["s12"][1:]
    s13lo, s13hi = io["s13"][1:]
    s23hi = read(inputs["direction"])["theta23_sin2"]["IO"][2]
    m2lo, m2hi = map(math.sqrt, io["atm_abs"][1:])
    solar_hi = io["solar"][2]
    m1ratio = math.sqrt(1-solar_hi/m2lo**2)
    beta_ratio_sq = (1-s13hi)*(1-(1-s12lo)*solar_hi/m2lo**2)
    generic = (1-s13hi)*m2lo*((1-s12hi)*m1ratio-s12hi)
    a = math.sqrt(s13hi/((1-s13hi)*(1-s23hi)))
    b = math.sqrt(s23hi/(1-s23hi))
    xenon = neutrino["kamland_zen"]
    four_nuclei = read(inputs["four_isotope_results"])["IBM2_plus_Mo"]["fixed"]["bound_eV"]
    results = {}
    for case in inputs["scenarios"]:
        q = case["rate_floor"]*(1-case["cancel"])**2
        meg = br_cap(me, mm, dip["mu_lifetime_s"], dip["meg_upper"], q, dip)
        tau_e = br_cap(me, mt, shared["tau_lifetime_lower_s"], inputs["tau_e_upper"], q, dip)
        tau_mu = br_cap(mm, mt, shared["tau_lifetime_lower_s"], shared["tau_mu_gamma_bound"], q, dip)
        gap = ((emu[0][0]-case["other_ae"])/(2*me**2) -
               (mu_upper+case["other_amu"])/(2*mm**2))
        assert gap > 0
        offdiag = meg*math.hypot(1/me, 1/mm) + (tau_e+tau_mu)/mt
        hmax = min(1., (offdiag+case["defect_per_m2_GeV_minus2"])/gap)
        ee_ratio = least_component(beta_ratio_sq, hmax, a, b)
        mass_floor = max(generic, m2lo*ee_ratio)
        phase_sine_sq = min(1., max(0., 1-(ee_ratio/(1-s13lo))**2) /
                                (4*s12lo*(1-s12lo)*m1ratio))
        lifetime_caps = [xenon["half_life_lower_years"]*(m/mass_floor)**2
                         for m in xenon["m_bb_upper_eV_endpoints_by_nme"]]
        hypothetical = inputs["comparison_mass_eV"]/m2lo
        required_h = least_component(beta_ratio_sq, hypothetical, b, a)
        results[case["label"]] = {
            "diagonal_response_gap_GeV_minus2":gap, "offdiag_allowance_GeV_minus2":offdiag,
            "emu_entry_ratio_to_largest_neutrino_mass_cap":hmax,
            "absolute_emu_entry_cap_eV":hmax*m2hi, "ee_entry_ratio_floor":ee_ratio,
            "effective_majorana_floor_eV":mass_floor,
            "relative_majorana_phase_abs_max_deg":2*math.degrees(math.asin(math.sqrt(phase_sine_sq))),
            "xenon_half_life_upper_years_by_NME_endpoint":lifetime_caps,
            "uniform_half_life_gain_to_test_even_with_one_extra_heavy_rate_IBM2":(four_nuclei/mass_floor)**2,
            "relative_antiphase_amplitude_min_to_reach_28meV":max(0.,1-0.028/mass_floor),
            "min_intertwining_defect_for_30meV_GeV_minus2":max(0.,required_h*gap-offdiag)}
    return {"generic_inverted_floor_eV":generic, "normalized_e_row_norm_sq_floor":beta_ratio_sq,
            "standard_inverted_upper_eV":(1-s13lo)*m2hi, "Rb_two_width_gap_can_vanish":emu[1][0] < 0,
            "cases":results}


def compare(saved, actual):
    if isinstance(actual, dict):
        assert saved.keys() == actual.keys()
        for key in actual:
            compare(saved[key], actual[key])
    elif isinstance(actual, list):
        assert len(saved) == len(actual)
        for left, right in zip(saved, actual):
            compare(left, right)
    else:
        assert math.isclose(saved, actual, rel_tol=4e-12, abs_tol=1.e-21)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--check", action="store_true")
    action.add_argument("--write-new", action="store_true")
    options = parser.parse_args()
    out = calculate()
    if options.write_new:
        with (HERE/"results.json").open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(out, indent=2, allow_nan=False)+"\n")
    elif options.check:
        compare(read("results.json"), out)
        print("1596: inverted Majorana mass from charged response gap and three LFV bounds verified")
    else:
        print(json.dumps(out, indent=2))
