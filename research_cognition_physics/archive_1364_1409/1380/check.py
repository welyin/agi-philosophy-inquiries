"""Round 1380: direct CP rate, arbitrary-order closure and aligned heavy-family bound."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def asymmetry(minus, plus):
    return (minus - plus) / (minus + plus)


def minimum_path_ratio(asym):
    return abs(asym) / (1 + math.sqrt(1 - asym**2))


def ratio_of_rates(asym):
    return (1 + asym) / (1 - asym)


def witness_matrix(parameters):
    s12, s23, s13 = (parameters[key] for key in ("s12", "s23", "s13"))
    c12, c23, c13 = (math.sqrt(1 - s**2) for s in (s12, s23, s13))
    phase = np.exp(1j * parameters["delta_rad"])
    r12 = np.array([[c12, s12, 0], [-s12, c12, 0], [0, 0, 1]])
    r13 = np.array([[c13, 0, s13 / phase], [0, 1, 0], [-s13*phase, 0, c13]])
    r23 = np.array([[1, 0, 0], [0, c23, s23], [0, -s23, c23]])
    return r23 @ r13 @ r12


def polygon_certificate(matrix, b, d):
    z = matrix[1] * matrix[0].conj()
    quartet = (z[b] * z[d].conj()).imag
    other_sides = sum(abs(value) for j, value in enumerate(z) if j not in (b, d))
    bound = min(abs(z[b]), abs(z[d])) * other_sides
    assert abs(z.sum()) < 5e-14 and abs(quartet) <= bound + 5e-14
    return float(quartet), float(bound)


def charge_witness(matrix, target):
    d, s, b = 0, 1, 2
    weak = np.array([
        abs(matrix[0, s])**2 * matrix[1, b] * matrix[1, d].conj(),
        abs(matrix[1, s])**2 * matrix[0, b] * matrix[0, d].conj(),
    ])
    both_charges = np.array([weak, weak.conj()])
    hadronic = np.linalg.solve(both_charges, np.sqrt([1+target, 1-target]))
    minus, plus = abs(both_charges @ hadronic)**2
    quartet, bound = polygon_certificate(matrix, b, d)
    strong_product = hadronic[0] * hadronic[1].conj()
    weak_product = weak[0] * weak[1].conj()
    assert abs(weak_product.imag - quartet * abs(matrix[0,s]*matrix[1,s])**2) < 1e-14
    assert abs(minus - plus + 4*weak_product.imag*strong_product.imag) < 1e-12
    four = np.eye(4, dtype=complex)
    four[:3, :3] = matrix
    assert polygon_certificate(four, b, d) == (quartet, bound)
    return {"quartet_im": quartet, "polygon_upper_bound": bound,
            "charge_rates_in_common_arbitrary_units": [float(minus), float(plus)],
            "charge_asymmetry": float(asymmetry(minus, plus)),
            "fourth_decoupled_row_same_quartet": True}


def arbitrary_phase(target, phi):
    coeffs = np.array([[1, np.exp(1j*phi)], [1, np.exp(-1j*phi)]])
    strong = np.linalg.solve(coeffs, np.sqrt([1+target, 1-target]))
    rates = abs(coeffs @ strong)**2
    return {"unit_rad_weak_phase": phi, "strong_magnitudes": abs(strong).tolist(),
            "asymmetry": float(asymmetry(*rates))}


def higgs_results(inputs):
    previous = read_json(ROOT / inputs["previous_higgs_results"])
    data = read_json(ROOT / inputs["previous_higgs_inputs"])["atlas_four_leptons"]
    norm = data["central_normalization"]
    eta = norm["inclusive_bbh_pb"] * norm["br_zz"] / norm["stage_gg_plus_bb_times_br_pb"]
    gg, vv = data["gg_plus_bb"], data["vbf"]
    k = inputs["sensitivity_multiplier"]
    upper = (gg["mu"]+k*math.hypot(*gg["width_plus"])) / (vv["mu"]-k*math.hypot(*vv["width_minus"]))
    amp = complex(*previous["SM_quark_loop_real_imag"])
    floor = {}
    for n in inputs["heavy_extra_pairs"]:
        rgg = abs(amp + 8*n/3)**2 / abs(amp)**2
        floor[str(n)] = (1-eta)*rgg + eta
    assert abs(floor["1"] - previous["min_gg_plus_bb_over_vbf_ratio"]) < 1e-12
    return {"atlas_central": gg["mu"]/vv["mu"], "atlas_rectangular_upper": upper,
            "family_floor": floor, "smallest_four_family_floor_over_upper": floor["1"]/upper}


def evaluate():
    inputs = read_json(ROOT / "inputs.json")
    obs = inputs["lhcb"]
    a = obs["asymmetry"]
    offset = inputs["sensitivity_multiplier"] * (obs["stat"]+obs["syst"])
    interval = [a-offset, a+offset]
    r = minimum_path_ratio(a)
    assert abs(asymmetry((1-r)**2, (1+r)**2) - a) < 1e-14
    rng = np.random.default_rng(1380)
    for n in (2, 3, 4, 6):
        for _ in range(24):
            matrix, _ = np.linalg.qr(rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n)))
            poly, _ = polygon_certificate(matrix, 0, 1)
            assert n != 2 or abs(poly) < 5e-14
    matrix = witness_matrix(inputs["three_by_three_existence_witness"])
    assert np.linalg.norm(matrix @ matrix.conj().T - np.eye(3)) < 1e-14
    return {"charge_asymmetry_sensitivity_interval": interval,
            "minus_over_plus_central": ratio_of_rates(a),
            "minus_over_plus_box": list(map(ratio_of_rates, interval)),
            "minimum_weaker_path_fraction_of_stronger": {"central": r, "whole_box": minimum_path_ratio(interval[1])},
            "reported_likelihood_significance_sigma": obs["reported_likelihood_sigma"],
            "additional_cp_source_r_exact_for_central": r,
            "three_family_existence": charge_witness(matrix, a),
            "arbitrarily_small_phase_witnesses": [arbitrary_phase(a, phi) for phi in inputs["small_weak_phases_rad"]],
            "aligned_higgs": higgs_results(inputs)}


def close_tree(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(close_tree(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(close_tree(x,y) for x,y in zip(a,b))
    if isinstance(a, (float, int)) and not isinstance(a, bool):
        return math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-13)
    return a == b


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    outcome = evaluate()
    result = ROOT / "results.json"
    if args.write_new:
        with result.open("x", encoding="utf-8") as stream:
            json.dump(outcome, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        assert close_tree(outcome, read_json(result)), "Stored and recomputed results differ."
    print("1380 OK")
