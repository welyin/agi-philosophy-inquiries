"""Read-only algebra checks for a count0, conditional soft-photon adoption."""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import sys


def calculate():
    c0, width = F(48, 25), F(1, 25)
    k3, k2, n = F(1), F(1), F(2)
    k1 = 36*n - 24*k3 - 18*k2
    kgamma = k2/2 + k1/36
    ratio = 2*kgamma/k3
    assert (k1, kgamma, ratio) == (F(30), F(4, 3), F(8, 3))
    assert F(2, 3)*k3 + kgamma == n

    gap = abs(ratio-c0)
    assert gap == F(56, 75)
    # For every n != 2, |r(n)-c0| >= 2-|r(2)-c0| > |r(2)-c0|.
    other_integer_lower_bound = 2-gap
    assert other_integer_lower_bound > gap
    lower_c, upper_c = c0-width, c0+width
    assert F(5, 3) < lower_c <= upper_c < F(11, 3)
    interval_gap = ratio-upper_c
    assert interval_gap == gap-width == F(53, 75) > 0
    # Demonstrate the conditional remainder bound, not a QCD error claim.
    remainder_budget = F(1, 10)
    conditional_gap = gap-width-remainder_budget
    assert conditional_gap == F(91, 150) > 0

    forbidden_r = F(2)
    forbidden_kgamma = forbidden_r*k3/2
    forbidden_lattice_value = F(2, 3)*k3+forbidden_kgamma
    assert forbidden_lattice_value == F(5, 3)
    assert forbidden_lattice_value.denominator != 1

    relaxed_k3, relaxed_kgamma = F(3), F(3)
    relaxed_r = 2*relaxed_kgamma/relaxed_k3
    assert F(2, 3)*relaxed_k3+relaxed_kgamma == F(5)
    assert relaxed_r-c0 == F(2, 25)
    assert (24*relaxed_k3+108)/36 == F(5)

    # Integer rephasing of a color fundamental of electric charge -1/3.
    # Reece's convention gives (delta K3, delta k_gamma) = (1, 1/3).
    # This is bare-coefficient transport, not a new physical K3 or E/N model.
    delta_k3, delta_kgamma = F(1), F(1, 3)
    shifted_k3, shifted_kgamma = k3+delta_k3, kgamma+delta_kgamma
    assert F(2, 3)*delta_k3+delta_kgamma == 1
    assert F(2, 3)*shifted_k3+shifted_kgamma == n+1
    shifted_bare_ratio = 2*shifted_kgamma/shifted_k3
    assert shifted_bare_ratio == F(5, 3) != ratio

    # Toy dimensionless normalization: ghat = (2*pi/alpha)*g.
    sqrt_chi = F(3)
    scaling = []
    for fa in (F(3), F(6)):
        mass = sqrt_chi/fa
        ghat = (ratio-c0)/fa
        assert mass*fa == sqrt_chi
        assert ghat/mass == gap/sqrt_chi
        scaling.append({"fa": str(fa), "mass": str(mass),
                        "normalized_g": str(ghat),
                        "normalized_g_over_mass": str(ghat/mass)})

    # LO two-flavor bookkeeping in Q_a=I/2 and mass-aligned Q_a bases.
    # This verifies transport of the already adopted matching, not QCD itself.
    z = F(12, 25)
    c_lo = F(2, 3)*(4+z)/(1+z)
    direct_symmetric = ratio-F(5, 3)
    meson_symmetric = -(1-z)/(1+z)
    qa_up = 1/(1+z)
    anomaly_mass_aligned = 6*(qa_up*F(4, 9)+(1-qa_up)*F(1, 9))
    assert anomaly_mass_aligned == c_lo
    assert direct_symmetric+meson_symmetric == ratio-c_lo
    assert c_lo != c0  # Do not conflate LO with the adopted NLO estimate.

    mu_gev2, alpha_display = 0.00570, 1/137
    central_floor = alpha_display*float(gap)/(2*math.pi*mu_gev2)
    assert 0.151 < central_floor < 0.153
    return {
        "scope": "count0 conditional algebra; no QCD calculation, fit, or empirical test",
        "scientific_increment": 0,
        "minimal_example": {"K3": str(k3), "K2": str(k2), "K1": str(k1),
                            "k_gamma": str(kgamma), "E_over_N": str(ratio)},
        "nearest_lattice_point": {
            "C0": str(c0), "gap": str(gap),
            "all_other_integers_lower_bound": str(other_integer_lower_bound),
            "adopted_interval_for_algebra_only": [str(lower_c), str(upper_c)],
            "interval_gap": str(interval_gap),
            "illustrative_remainder_budget": str(remainder_budget),
            "conditional_gap_after_both_budgets": str(conditional_gap)},
        "E_over_N_2_at_K3_1_rejected": {"lattice_value": str(forbidden_lattice_value)},
        "K3_3_boundary_example": {"K3": str(relaxed_k3), "k_gamma": str(relaxed_kgamma),
                                 "E_over_N": str(relaxed_r), "central_bracket": str(relaxed_r-c0)},
        "global_integer_rephasing_bookkeeping": {
            "delta_K3": str(delta_k3), "delta_k_gamma": str(delta_kgamma),
            "shifted_bare_K3": str(shifted_k3), "shifted_bare_k_gamma": str(shifted_kgamma),
            "shifted_bare_ratio_not_physical_matching_ratio": str(shifted_bare_ratio),
            "mass_phase_and_derivative_transport_required": True},
        "dimensionless_scaling_examples": scaling,
        "LO_basis_transport": {"z": str(z), "C_LO": str(c_lo),
                               "symmetric_direct": str(direct_symmetric),
                               "symmetric_meson": str(meson_symmetric),
                               "total_in_both_bases": str(ratio-c_lo)},
        "central_display_only": {"alpha": "1/137", "mass_times_fa_GeV2": mu_gev2,
                                 "minimum_abs_g_over_mass_GeV_minus2": round(central_floor, 9)},
        "all_assertions_passed": True}


if __name__ == "__main__":
    if sys.argv[1:]:
        raise SystemExit("No arguments: this checker is strictly read-only.")
    expected = json.loads(Path(__file__).with_name("results.json").read_text(encoding="utf-8"))
    assert calculate() == expected
    print("PASS: conditional soft-photon/mass overlap and exact lattice bookkeeping; count0.")
