"""Born barrel production and an inclusive, interference-safe counting test."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def poisson_upper(count, tail):
    """Invert the exact lower Poisson tail by recurrence and bisection."""
    def cdf(mean):
        term = math.exp(-mean)
        terms = [term]
        for j in range(1, count + 1):
            term *= mean / j
            terms.append(term)
        return math.fsum(terms)

    low, high = float(count), max(10.0, 2.0 * (count + 1))
    assert cdf(low) > tail > cdf(high)
    for _ in range(64):
        mid = (low + high) / 2
        if cdf(mid) > tail:
            low = mid
        else:
            high = mid
    return high


def barrel_pb(cfg, energy, mass, coupling, a):
    s = energy**2
    r = mass**2 / s
    integral = 2 * ((2 + 4*r/(1-r)**2) * math.atanh(a) - a)
    return (coupling**2 * (1-r) * integral / (4*s*cfg["alpha_inverse_born"])
            * cfg["GeV_inverse_squared_to_pb"])


def final_sample(cfg):
    data = cfg["final_sample"]
    a = data["abs_cos_edges"][-1]
    low_mass, high_mass = cfg["candidate_mass_GeV"]
    low_g, high_g = cfg["family_coupling"]
    rows = data["energy_low_high_GeV_and_luminosity_pb_inverse"]
    small, large = [], []
    for low_energy, high_energy, lum in rows:
        small.append(lum * barrel_pb(cfg, high_energy, low_mass, low_g, a))
        large.append(lum * barrel_pb(cfg, low_energy, high_mass, high_g, a))
        min_scaled_pt = (1 - (high_mass / low_energy)**2) * math.sqrt(1-a*a) / 2
        assert min_scaled_pt > data["high_energy_min_fraction_sqrt_s_photon_pt"]
    low, high = math.fsum(small), math.fsum(large)
    upper = poisson_upper(sum(data["observed"]), cfg["one_sided_alpha"])
    caps = []
    for background in data["coherent_caps"]:
        signal = (math.sqrt(upper) + math.sqrt(background))**2
        caps.append(dict(background_upper=background, selected_signal_upper=signal,
                         product_upper=signal/low, minimal_branch_transfer_upper=3*signal/low))
    return dict(luminosity_pb_inverse=math.fsum(r[2] for r in rows), nominal_yield=[low, high],
                upper_total_mean=upper, orthogonal_product_upper=upper/low, coherent=caps)


def calculate(cfg):
    energy = cfg["sqrt_s_GeV"]
    s = energy**2
    angle = math.radians(cfg["barrel_theta_min_degrees"])
    a = math.cos(angle)
    corners = []
    for mass in cfg["candidate_mass_GeV"]:
        r = mass**2 / s
        energy_photon = energy * (1 - r) / 2
        assert energy_photon * math.sin(angle) > cfg["single_photon_energy_cut_GeV"]
        for coupling in cfg["family_coupling"]:
            sigma = barrel_pb(cfg, energy, mass, coupling, a)
            corners.append(dict(mass_GeV=mass, g_F=coupling, photon_energy_GeV=energy_photon,
                                barrel_pb=sigma, nominal_events=sigma*cfg["integrated_luminosity_pb_inverse"]))

    smallest_yield = min(c["nominal_events"] for c in corners)
    upper = poisson_upper(cfg["observed_barrel_count"], cfg["one_sided_alpha"])
    baseline = cfg["nominal_neutrino_count"]
    scenarios = []
    for mult in cfg["background_cap_multipliers"]:
        cap = (math.sqrt(upper) + math.sqrt(mult * baseline))**2
        scenarios.append(dict(background_cap_multiplier=mult, selected_signal_cap=cap,
                              effective_invisible_product_cap=cap/smallest_yield,
                              minimum_transfer_for_entire_minimal_branch_exclusion=3*cap/smallest_yield))
    return dict(final_sample=final_sample(cfg), corners=corners, upper_total_mean=upper, orthogonal_product_cap=upper/smallest_yield,
                coherent=scenarios, minimal_active_branching_fraction=1/3,
                minimal_width_over_mass=[g**2/(16*math.pi) for g in cfg["family_coupling"]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = calculate(json.loads((HERE / "inputs.json").read_text(encoding="utf-8")))
    path = HERE / "results.json"
    if args.check:
        if result != json.loads(path.read_text(encoding="utf-8")):
            raise SystemExit("1627: saved results differ")
        print("1627 PASS: barrel cross section, recoil energy, count inversion and coherent norm bounds")
    else:
        path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
