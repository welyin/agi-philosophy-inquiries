"""Independent Koide roots, U(3) current contraction and two distinct experiments."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


OUT = Path(__file__).with_name("independent_results.json")
ME, MMU = 0.51099895069, 105.6583755
GF_TEV = 11.663787


def koide_roots():
    # In units of the physical muon mass, take sqrt of the third mass as unknown.
    v = math.sqrt(ME / MMU)
    coefficients = [1., -4 * (v + 1), 3 * (v*v + 1) - 2 * (v + 1)**2]
    masses = sorted(float(x*x*MMU) for x in np.roots(coefficients))
    for m in masses:
        a = np.sqrt([ME, MMU, m])
        assert abs((a @ a) / sum(a)**2 - 2 / 3) < 1e-13
    return masses


def currents():
    # Explicit matrix units for both real components of A12.
    e12 = np.array([[0., 1.], [0., 0.]])
    tl = [(e12 + e12.T)/2, (e12 - e12.T)/(2j)]
    tr = [-a.conj() for a in tl]
    # Coefficient of (mu_L e_L)(mu_R e_R) in -g^2 J^2/(2M^2).
    lr = -sum(x[1, 0] * y[1, 0] for x, y in zip(tl, tr))
    same_chirality = -sum(x[1, 0]**2 for x in tl) / 2
    vectorlike = -sum(x[1, 0]**2 for x in tl)
    assert abs(lr - .5) < 1e-13 and abs(same_chirality) < 1e-13 and abs(vectorlike) < 1e-13
    diag = np.diag([1., 0.]) / math.sqrt(2)
    return dict(coefficient_LR_times_M12_sq_over_g_sq=lr.real,
                axial_squared_times_M11_sq_over_g_sq=-diag[0, 0]**2/2,
                vectorlike_transfer=vectorlike.real)


def atomic_response():
    # SI lifetime to GeV^-1; use the reduced mass instead of m_e in the Bohr radius.
    alpha = 1 / 137.035999177
    reduced = ME * MMU / (ME + MMU) / 1000
    hbar = 6.582119569e-25
    lifetime_seconds = 2.1969811e-6
    gf_gev = GF_TEV / 1e6
    delta = 8 * gf_gev / math.sqrt(2) * (alpha * reduced)**3 / math.pi
    response = 1.5 * (delta * lifetime_seconds / hbar)**2
    # This leading-order atomic convention is a numerical cross-check of the published prefactor.
    return response


def compute():
    ratio = 2 * ME / (ME + MMU)
    # Enclose the D<=0.001 Belle II slope range in |s|<=0.0004.
    mu_over_e = MMU / ME
    rho_hi = 2 / (1 + mu_over_e**(1-.0004) * math.exp(-.001))
    rho_lo = 2 / (1 + mu_over_e**(1+.0004) * math.exp(.001))
    assert .0095961 < rho_lo < ratio < rho_hi < .0096562
    assert (rho_hi / ratio)**2 < 1.006255
    coefficient = 2 * math.pi / 14.0**2
    g_relative = coefficient * ratio / (math.sqrt(2) * GF_TEV)
    vacuum_factor, field_factor = 1.96e-5, .78
    # Published mixed-chirality spin transfer in a 0.1 tesla field.
    p = field_factor * vacuum_factor * g_relative**2
    reduced_response = atomic_response()
    assert 0.95 < reduced_response / vacuum_factor < 1.05
    g_psi = math.sqrt(8.3e-11 / (field_factor * vacuum_factor))
    gauge = 2 * math.sqrt(4 * math.pi / 137.035999177)
    m11 = gauge / (2 * math.sqrt(coefficient))
    m12 = m11 / math.sqrt(ratio)
    prop = (0.209 / m11)**2
    examples = {}
    for multiplier in (1, 10):
        m_e_family = 3. * multiplier
        c_aa = gauge**2 / (4 * m_e_family**2)
        g_lr = gauge**2 * math.sqrt(2) / (8 * (m_e_family**2 / ratio))
        assert abs(g_lr - c_aa * ratio / math.sqrt(2)) < 1e-14
        examples[str(multiplier)] = dict(
            M11_TeV=m_e_family, M12_TeV=m_e_family / math.sqrt(ratio),
            contact_TeV_minus_2=c_aa,
            p_at_0p1T=field_factor*vacuum_factor*(g_lr/GF_TEV)**2)
    return {
        "koide_roots_MeV": koide_roots(), "gauge_current_contraction": currents(),
        "M11_squared_over_M12_squared": ratio,
        "LEP_destructive_AA_only_95_coefficient_TeV_minus_2": coefficient,
        "M12_over_M11": 1 / math.sqrt(ratio),
        "G12_over_GF_from_Bhabha": g_relative, "G12_over_GF_from_PSI_only": g_psi,
        "P_under_historical_response": p, "P_under_reduced_mass_response":
            field_factor*reduced_response*g_relative**2,
        "atomic_free_probability_coefficient_from_constants": reduced_response,
        "PSI_direct_probability_over_prediction": 8.3e-11 / p,
        "amplitude_multiplier_needed_to_reach_PSI_limit": math.sqrt(8.3e-11/p),
        "stress_probability_with_coeff_1p20_ratio_1p10_amplitude_1p05":
            p * (1.20 * 1.10 * 1.05)**2,
        "conditional_mass_lower_if_gF_equals_2e0_TeV": [m11, m12],
        "largest_s_channel_inverse_propagator_error_at_that_mass": prop/(1-prop),
        "same_light_mass_ratios_and_different_flavor_exchange": examples,
    }


def clean(value):
    if isinstance(value, dict):
        return {key: clean(val) for key, val in value.items()}
    if isinstance(value, list):
        return [clean(val) for val in value]
    return float(f"{value:.11g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    opts = parser.parse_args()
    result = clean(compute())
    if opts.write:
        with OUT.open("x", encoding="utf-8") as f:
            f.write(json.dumps(result, indent=2) + "\n")
    elif result != json.loads(OUT.read_text(encoding="utf-8")):
        raise SystemExit("1619 independent results do not reproduce")
    print("1619 roots, left/conjugate-right currents, Bhabha and muonium passed")
