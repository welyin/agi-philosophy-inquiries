"""Three Dirac widths at the Z pole and the associated-production radial bound."""
import argparse
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent


def read(relative):
    return json.loads((HERE / relative).read_text(encoding="utf-8-sig"))


def invis_fraction(x):
    return 0 if x >= .25 else (1 - x) * math.sqrt(1 - 4 * x)


def maximal_dirac_moment():
    data = read("../1422/inputs.json")
    n, err = data["Z_reanalysis"]
    floor = n - 2 * err
    max_transfer, unseen_cap = 1.01, .02
    loss = 3 - (floor - unseen_cap) / max_transfer
    assert 0 < loss < 1
    left, right = 0., .25
    for _ in range(70):
        mid = (left + right) / 2
        if 1 - invis_fraction(mid) < loss:
            left = mid
        else:
            right = mid
    assert right < (3 - math.sqrt(5)) / 4
    mz = data["MZ_GeV"]
    mass = mz * math.sqrt(right)
    return dict(input_floor=floor, max_transfer=max_transfer, extra_width_cap=unseen_cap,
                loss_cap=loss, largest_pole_mass=mass, x=right,
                maximal_fourth_moment_4sum=4 * mass**4)


def test_mass(floor, amp, kappa2, available):
    h4 = 2 * floor**4
    cap = amp / (h4 - available)
    extra = h4 - available - amp / kappa2
    scale = (h4 / (amp / kappa2 + available))**.25
    return dict(charged_pole_floor_GeV=floor, kappa2_upper_no_extra=cap,
                extra_negative_or_matching_GeV4=extra,
                minimum_single_Dirac_equivalent_GeV=(max(0, extra) / 4)**.25,
                minimum_charged_pole_over_GW=scale,
                largest_WH_ratio_transfer_separating=.55 * .84 / cap)


def certificate():
    neutral = maximal_dirac_moment()
    wh, low, _ = read("../1470/inputs.json")["associated_VHbb"]["WH"]
    mu_floor = wh - 2 * low
    # Distinct inclusive fractions bb, tau tau, cc, inclusive four fermions.
    ref_parts = [.577, .0632, .0291, .215, .0264]
    ref_floor = .84
    assert sum(ref_parts) > ref_floor
    transfer = 1.10
    kappa2 = mu_floor * ref_floor / transfer
    amp = 8 * math.pi**2 * 249**2 * 126**2
    # Color triplet Dirac top + 5 color triplet quarks + 3 charged leptons.
    # Deliberately allow all eight non-top species to have color multiplicity 3.
    ordinary = 12 * 175**4 + 8 * 12 * 5**4
    # Allow the running Dirac-neutrino masses to exceed the LEP pole by 10%.
    available = ordinary + 1.1**4 * neutral["maximal_fourth_moment_4sum"]
    return dict(neutrinos=neutral, SM_reference_fraction_sum=sum(ref_parts),
                normalized_partial_width_floor=ref_floor, WHbb_two_width_floor=mu_floor,
                WH_production_times_differential_width_transfer_cap=transfer,
                kappa2_empirical_lower=kappa2, kappa_empirical_lower=math.sqrt(kappa2),
                radial_factor_GeV4=amp, ordinary_negative_moment_cap=ordinary,
                registered_total_fermion_cap_GeV4=available,
                historical_type_II_2017=test_mass(570., amp, kappa2, available),
                updated_type_II_2026=test_mass(670., amp, kappa2, available))


def same(lhs, rhs):
    if isinstance(lhs, dict):
        return lhs.keys() == rhs.keys() and all(same(v, rhs[k]) for k, v in lhs.items())
    return math.isclose(lhs, rhs, rel_tol=2e-11, abs_tol=5e-13)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write", action="store_true")
    args = cli.parse_args()
    record = certificate()
    if args.write:
        with (HERE / "independent_results.json").open("x", encoding="utf-8") as f:
            json.dump(record, f, ensure_ascii=False, indent=2)
            f.write("\n")
    else:
        assert same(record, read("independent_results.json"))
        print("PASS: 1589 Dirac pole, radial WHbb and two charged-Higgs analyses.")
