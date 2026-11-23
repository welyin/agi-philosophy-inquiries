"""Lepton pole and flavour-changing Z tests for the complete charged spectrum."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def extremal_diagonal(rank, caps):
    remaining = float(rank)
    p = np.zeros(len(caps))
    for i in np.argsort(caps)[::-1]:
        p[i] = min(caps[i], remaining)
        remaining -= p[i]
    if remaining > 1e-14:
        return None
    return p


def rank_bounds(caps):
    bounds = {}
    for rank in (1, 2, 3):
        p = extremal_diagonal(rank, caps)
        if p is None:
            bounds[str(rank)] = {"feasible_under_caps": False}
            continue
        if rank == 1:
            vector = np.sqrt(p)
            proj = np.outer(vector, vector)
        else:
            vector = np.sqrt(1 - p)
            proj = np.eye(3) - np.outer(vector, vector)
        off = proj - np.diag(np.diag(proj))
        minimum = float(rank - p @ p)
        assert np.allclose(proj @ proj, proj) and np.allclose(np.diag(proj), p)
        assert abs(np.linalg.norm(off, "fro")**2 - minimum) < 1e-13
        bounds[str(rank)] = {"feasible_under_caps": True, "extremal_p": p.tolist(),
                             "ordered_offdiagonal_square_floor": minimum}
    return bounds


def run():
    data = json.loads((HERE / "inputs.json").read_bytes())
    prior = json.loads((HERE / data["left_charged_lepton_pole"]).read_bytes())
    bounds = prior["conventions"]
    smin, smax = bounds["weak_mixing_sin_squared_bounds"]
    kmin, kmax = bounds["common_Z_normalization_bounds"]
    nuisance = bounds["lepton_specific_absolute_budget"]
    widths = data["pole_marginal_widths"]
    right = np.array(data["gR_absolute_e_mu_tau"])
    left = np.array([prior["lep_sld"]["absolute_gL"][f] for f in ("electron", "muon", "tau")])
    left_effective_min = (left[:, 0] - widths * left[:, 1] - nuisance) / kmax
    assert np.all(left_effective_min > smax)
    left_trace_min = (2 * (smin + left_effective_min)).sum()
    caps = 2 * (smax + (right[:, 0] + widths * right[:, 1] + nuisance) / kmin)
    effective_min = (right[:, 0] - widths * right[:, 1] - nuisance) / kmax
    low_branch_max = 2 * (smax - effective_min)
    high_branch_min = 2 * (smin + effective_min)
    trace_windows = []
    for signs in itertools.product((False, True), repeat=3):
        lower = sum(high_branch_min[i] for i in range(3) if signs[i])
        upper = sum(caps[i] if signs[i] else low_branch_max[i] for i in range(3))
        trace_windows.append((sum(signs), float(lower), float(upper)))
    two_high_max = max(hi for n, lo, hi in trace_windows if n <= 2)
    three_high_min = min(lo for n, lo, hi in trace_windows if n == 3)
    assert two_high_max < 2 < three_high_min and max(caps) < 1
    pole_solutions = [(l, r, ell) for l, r, ell in itertools.product(range(4), range(4), (0, 1))
                      if l + r == 3 and l - r == 3 - ell]
    assert pole_solutions == [(3, 0, 0)]
    ranks = rank_bounds(caps)
    floor = min(item["ordered_offdiagonal_square_floor"] for item in ranks.values()
                if item["feasible_under_caps"])
    electron = data["electron_width_GeV"]
    total = data["Z_width_GeV"]
    ebranch_floor = (electron[0] - widths * electron[1]) / (total[0] + widths * total[1])
    e_coupling_square = sum((column[0, 0] + widths * column[0, 1])**2 for column in (left, right))
    coefficient = ebranch_floor * kmin**2 / (4 * e_coupling_square)
    predicted = coefficient * floor
    obs_upper = sum(data["lfv_branch_upper_e_mu_e_tau_mu_tau"])
    with_triplet = [sum(lam * (lam - max(caps)) for lam in eigenvalues)
                    for eigenvalues in itertools.combinations_with_replacement((0, 1, 2), 3)
                    if 2 in eigenvalues]
    triplet_floor = min(with_triplet)
    assert triplet_floor > 2 and triplet_floor > floor
    witness_p = np.array([0.03, 0.94, 0.03])
    witness_right = kmax * abs(smax - witness_p / 2)
    witness_left = kmax * abs(smax - 0.5)
    assert np.all(abs(witness_right - right[:, 0]) < nuisance)
    assert np.all(abs(witness_left - left[:, 0]) < nuisance)
    witness = {"right_rank": 1, "diagonal_probabilities": witness_p.tolist(),
               "diagonal_right_moduli": witness_right.tolist(), "left_modulus": witness_left,
               "ordered_offdiagonal_square": float(1 - witness_p @ witness_p),
               "quark_doublet_multiplicities_L_R": [2, 0]}
    pairs = [(l, r) for l, r in itertools.product(range(4), repeat=2) if l-r == 3]
    assert 2 < left_trace_min < 3 and pairs == [(3, 0)]
    m = data["electromagnetic_example_mass_cap_GeV"]
    energy = data["electromagnetic_example_cm_GeV"]
    t = (m / energy)**2
    assert t < 0.25
    kernel = math.sqrt(1-4*t) * (1+2*t)
    return {"left_trace_floor": float(left_trace_min), "right_diagonal_caps": caps.tolist(),
            "projector_rank_bounds": ranks, "electron_branch_lower": ebranch_floor,
            "coefficient_before_transfer": coefficient, "lfv_branch_sum_floor_at_unit_transfer": predicted,
            "lfv_marginal_limit_sum": obs_upper, "necessary_rate_transfer_for_nonzero_rank_at_most": obs_upper / predicted,
            "quark_chiral_doublets_from_anomaly": pairs, "diagonal_only_competitor": witness,
            "triplet_ordered_offdiag_square_floor": triplet_floor,
            "triplet_branch_sum_floor_at_unit_transfer": coefficient * triplet_floor,
            "right_trace_ceiling_at_most_two_high_signs": two_high_max,
            "right_trace_floor_three_high_signs": three_high_min,
            "pole_width_only_solution_qL_qR_leptonR": pole_solutions,
            "known_massless_electromagnetic_weight": 8,
            "finite_mass_example_min_extra_photon_R": kernel / 3}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    path = HERE / "results.json"
    if args.check:
        assert json.loads(path.read_bytes()) == json.loads(output)
    else:
        path.write_text(output, encoding="utf-8")
    print(output)
