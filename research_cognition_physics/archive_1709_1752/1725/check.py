"""Charged-proton null test on the unsplit tens of the split-matter benchmark."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
def read(file):
    return json.loads((HERE/file).read_bytes())


def gauge_transfer_below_soft(flow, new_masses):
    start, soft = flow["MZ"], flow["MS"]
    thresholds = list(flow["physics"]["thresholds"])
    thresholds.extend((mass, np.array([4/15,0,2/3])) for mass in new_masses)
    ends = sorted({start,soft,*[m for m,_ in thresholds if start < m < soft]})
    squared = np.array([flow["physics"]["weak"][v][0]**2 for v in ("g1_GUT","g2","g3")])
    anomalous = np.array([23/10,9/2,4.])
    log_gain = 0.
    for lo, hi in zip(ends,ends[1:]):
        beta = np.array([41/10,-19/6,-7.])
        for mass, weight in thresholds:
            if mass <= lo:
                beta += weight
        logspan = math.log(hi/lo)/(16*math.pi**2)
        drift = 2*beta*squared*logspan
        integral = squared*logspan*np.array([-math.log1p(-x)/x if x else 1. for x in drift])
        log_gain += np.dot(anomalous,integral)
        squared = squared/(1-drift)
    benchmark = read("../1721/results.json")["one_light_Higgs_below_30TeV"]["actual_active_columns"]
    np.testing.assert_allclose(np.sqrt(squared), benchmark["gauge"], rtol=2e-8)
    return math.exp(log_gain)


def compute():
    config = read("inputs.json")
    frame = runpy.run_path(str(HERE/config["flavour"]))
    flow, rep = frame["charged_background"](.2)
    down = frame["unpack"](rep["down_texture_10PeV"])
    electron = frame["unpack"](rep["charged_texture_10PeV"]).T
    up = np.array(rep["up_texture_diag_10PeV"])
    assert np.argmin(up) == 0
    d = np.linalg.svd(down)[0][:, -1]
    e = np.linalg.svd(electron)[0][:, -1]
    weighted = d * np.array([2,1,1])
    terms = e.conj() * weighted
    lengths = abs(terms)
    leader = np.argmax(lengths)
    amplitude_min = max(0., 2*lengths[leader] - sum(lengths))
    assert amplitude_min > 0
    flip = np.where(np.arange(3)==leader, 1., -1.)
    phases = flip*np.exp(-1j*np.angle(terms))
    np.testing.assert_allclose(abs(terms @ phases), amplitude_min)
    norm = np.linalg.norm(weighted)
    np.testing.assert_allclose(norm**2, 1+3*abs(d[0])**2)
    angle = math.asin(amplitude_min/norm)
    # With an additional projective rotation of the physical right-electron line
    # by eps, the worst possible amplitude has an exact spherical-angle minimum.
    def amplitude_after_rotation(degrees):
        return norm*math.sin(max(0, angle - math.radians(degrees)))

    gauge = read(config["gauge_point"])
    reference_mass = gauge["baseline_unification_GeV"]
    g2 = 4*math.pi/gauge["central_conditional_combined_point"]["common_inverse_coupling"]
    data = read(config["nucleon_data"])
    pion = data["pion_W0_GeV2"]
    W0 = abs(pion[0])/math.sqrt(2)
    Wlo = (abs(pion[0])-sum(pion[1:]))/math.sqrt(2)
    mp, mpi = config["proton_mass_GeV"], config["neutral_pion_mass_GeV"]
    width_limit = data["hbar_GeV_s"] / (config["electron_mode_partial_lifetime_lower_years"] * data["year_seconds"])
    prefactor = mp/(32*math.pi)*(1-(mpi/mp)**2)**2
    def limits(w, amplitude):
        if amplitude <= 0:
            raise ValueError("A vanishing amplitude gives no lower bound")
        coef = math.sqrt(width_limit/prefactor)/w
        minimum_mass = math.sqrt(g2*amplitude/coef)
        amplitude_ratio = (minimum_mass/reference_mass)**2
        return dict(low_energy_single_helicity_coefficient_upper_GeVm2=coef,
                    minimum_MX_GeV_at_unit_transfer=minimum_mass,
                    minimum_MX_over_crossing_at_unit_transfer=minimum_mass/reference_mass,
                    rate_over_experimental_upper_at_crossing_and_unit_transfer=amplitude_ratio**2,
                    maximum_abs_transfer_at_crossing=1/amplitude_ratio,
                    minimum_same_helicity_amplitude_cancellation_fraction_at_unit_transfer=max(0.,1-1/amplitude_ratio))
    central = limits(W0, amplitude_min)
    lower = limits(Wlo, amplitude_min)
    gauge_point = gauge["central_conditional_combined_point"]
    soft_gain = gauge_transfer_below_soft(flow,gauge_point["five_color_poles_GeV"])
    low_gain = 1.25 * soft_gain
    target_amp = central["low_energy_single_helicity_coefficient_upper_GeVm2"]*reference_mass**2/g2
    rotation_to_bound = max(0, angle - math.asin(min(1., target_amp/norm)))
    return dict(scale_of_family_tags_GeV=1e7,gauge_crossing_GeV=reference_mass,gauge_g_squared=g2,
                abs_e_10_components=abs(e).tolist(),abs_down_10_components=abs(d).tolist(),
                abs_e_terms=lengths.tolist(),amplitude_tagged=abs(sum(terms)),
                amplitude_minimum_diagonal_phases=amplitude_min,amplitude_maximum_diagonal_phases=sum(lengths),
                length_in_full_lepton_space=norm,minimum_additional_electron_line_rotation_to_zero_deg=math.degrees(angle),
                minimum_additional_electron_line_rotation_to_fit_at_unit_transfer_central_W_deg=math.degrees(rotation_to_bound),
                lattice_central_and_lower_GeV2=[W0,Wlo],central_lattice=central,
                one_width_each_lower_lattice=lower,
                gauge_only_transfer_MZ_to_soft_for_downward_Wilson=soft_gain,
                conditional_full_low_transfer_with_long_distance_1p25=low_gain,
                maximum_abs_residual_transfer_above_soft_central_and_lower_W=
                    [x["maximum_abs_transfer_at_crossing"]/low_gain for x in (central, lower)],
                central_10deg_unregistered_rotation=limits(W0, amplitude_after_rotation(10)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = json.loads(json.dumps(compute()), parse_float=lambda x: float(f"{float(x):.10g}"))
    file = HERE/"results.json"
    if args.check:
        assert json.loads(file.read_bytes()) == result
        print("1725 verified the charged-antilepton channel family coefficient, phase/rotation bounds and the mass/transport threshold.")
    else:
        with file.open("x", encoding="utf8") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
