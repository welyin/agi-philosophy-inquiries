"""Check adopted BBN predictive summaries; no nuclear-network computation.

Default: read inputs/results and verify without writing.
--write-results: create results.json once, refusing to overwrite an existing file.
Only the Python standard library is needed. Scientific increment is zero.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist


HERE = Path(__file__).resolve().parent


def calculate(raw: bytes) -> dict:
    data = json.loads(raw.decode("utf-8"))
    assert data["scientific_increment"] == 0
    assert data["new_research_round"] is None
    calibration = data["calibration"]
    assert not calibration["uses_held_out_DH_data"]
    assert calibration["uses_BBN_helium_relation"]
    observed = data["held_out_observation"]
    assert not observed["used_in_prediction_calibration"]
    assert math.isfinite(observed["mean"]) and observed["sigma"] > 0
    contract = data["statistical_contract"]
    assert contract["prediction_summary_approximated_as_gaussian"]
    assert contract["observation_summary_approximated_as_gaussian"]
    assert contract["prediction_observation_errors_assumed_independent"]
    assert contract["do_not_add_omega_b_error_again"]
    assert not contract["method_averaging_allowed"]
    assert not contract["parameter_refitting_allowed"]
    coverage = contract["residual_interval_coverage"]
    assert 0 < coverage < 1
    quantile = NormalDist().inv_cdf((1 + coverage) / 2)

    rows = []
    for prediction in data["predictions"]:
        assert prediction["calibration_id"] == calibration["id"]
        assert prediction["omega_b_uncertainty_already_included"]
        assert prediction["nuclear_rate_and_neutron_lifetime_uncertainty_included"]
        assert not prediction["uses_held_out_DH_in_rate_fit_or_calibration"]
        assert math.isfinite(prediction["mean"])
        assert math.isfinite(prediction["sigma"]) and prediction["sigma"] > 0
        residual = observed["mean"] - prediction["mean"]
        sigma = math.hypot(prediction["sigma"], observed["sigma"])
        z = residual / sigma
        half_width = quantile * sigma
        lower, upper = residual - half_width, residual + half_width
        rows.append({
            "prediction_id": prediction["id"],
            "prediction_mean_yD": prediction["mean"],
            "prediction_sigma_yD": prediction["sigma"],
            "observation_mean_yD": observed["mean"],
            "observation_sigma_yD": observed["sigma"],
            "residual_observation_minus_prediction_yD": residual,
            "combined_sigma_yD": sigma,
            "signed_standardized_residual": z,
            "two_sided_gaussian_p": math.erfc(abs(z) / math.sqrt(2)),
            "residual_interval_95_yD": [lower, upper],
            "interval_contains_zero": lower <= 0 <= upper,
        })
    assert len({row["prediction_id"] for row in rows}) == len(rows)
    return {
        "schema": "cosmology_public_summary_check_results_v1",
        "inputs_sha256": hashlib.sha256(raw).hexdigest(),
        "scientific_increment": 0,
        "new_research_round": None,
        "scope": "Arithmetic verification of adopted predictive and observation summaries only",
        "normal_quantile_97_5_percent": quantile,
        "residual_definition": data["observable"]["residual_definition"],
        "rows": rows,
        "averaged_prediction": None,
        "physical_prediction_independently_recomputed": False,
        "joint_B993_cosmological_history_verified": False,
    }


def compare(saved, actual, location="results") -> None:
    if isinstance(actual, dict):
        assert isinstance(saved, dict) and saved.keys() == actual.keys(), location
        for key in actual:
            compare(saved[key], actual[key], f"{location}.{key}")
    elif isinstance(actual, list):
        assert isinstance(saved, list) and len(saved) == len(actual), location
        for index, (old, new) in enumerate(zip(saved, actual)):
            compare(old, new, f"{location}[{index}]")
    elif isinstance(actual, float):
        assert type(saved) in (int, float) and math.isfinite(saved), location
        assert math.isclose(saved, actual, rel_tol=1e-13, abs_tol=1e-15), location
    else:
        assert type(saved) is type(actual) and saved == actual, location


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = calculate((HERE / "inputs.json").read_bytes())
    destination = HERE / "results.json"
    if args.write_results:
        with destination.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(result, handle, indent=2, ensure_ascii=False, allow_nan=False)
            handle.write("\n")
    else:
        compare(json.loads(destination.read_text(encoding="utf-8")), result)
    print(json.dumps({
        "status": "created" if args.write_results else "verified_read_only",
        "comparisons": len(result["rows"]),
        "scientific_increment": 0,
        "GP_combined_sigma": result["rows"][0]["combined_sigma_yD"],
        "GP_z": result["rows"][0]["signed_standardized_residual"],
        "GP_two_sided_p": result["rows"][0]["two_sided_gaussian_p"],
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
