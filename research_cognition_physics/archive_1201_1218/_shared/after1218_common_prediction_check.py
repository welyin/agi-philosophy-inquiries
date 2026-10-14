"""Recombine accepted 1205/1217 certificates; no new model or channel claim.
All numerical decisions use Fraction. --write creates the result once; default is read-only.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = STAGE.parent
OUT = HERE / "after1218_common_prediction_results.json"

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))

def floor_strict(q: F, denominator: int = 10**9) -> F:
    p = F((q.numerator * denominator) // q.denominator, denominator)
    if p == q:
        p -= F(1, denominator)
    assert p < q
    return p

def build():
    sources = {}
    for number in (1205, 1217):
        folder = STAGE / str(number)
        acceptance = load(folder / "acceptance.json")
        assert acceptance["status"] == "ACCEPTED_WITH_SCOPE"
        for filename in ("proof.md", "results.json"):
            path = folder / filename
            actual = sha(path)
            assert actual == acceptance["files_sha256"][filename]
            sources[str(path.relative_to(ROOT)).replace("\\", "/")] = actual
        sources[str((folder / "acceptance.json").relative_to(ROOT)).replace("\\", "/")] = sha(folder / "acceptance.json")
    r5 = load(STAGE / "1205/results.json")
    r17 = load(STAGE / "1217/results.json")
    for rel, expected in r5["source_sha256"].items():
        path = STAGE / rel
        assert sha(path) == expected
        sources[str(path.relative_to(ROOT)).replace("\\", "/")] = expected
    r2 = load(STAGE / "1202/results.json")
    r4 = load(STAGE / "1204/results.json")
    cold = [F(s) for s in r2["one_marginal_rational"]]
    hot = [F(s) for s in r4["new_reference"]["one_marginal_rational"]]
    assert all(p > 0 for p in cold + hot)
    assert sum(cold) == sum(hot) == 1
    joint = sorted(a * b for a in cold for b in hot)
    tail = sum(joint[:12])
    assert tail == F(r17["groups"]["actual_A_and_K_spectral_bounds"]["exact_unital_tail"])
    # Distinct ideal source pools are independent conditional on the conserved
    # original-label success event. Correlations with retained inventory/R remain.
    # This is not an assertion that the real post-preparation environment is reset.
    failures = []
    dimensions = []
    for prep in (r2["preparation"], r4["new_reference"]):
        dimension = int(prep["dimension"])
        bits = int(prep["seed_bits_each"])
        trials = int(prep["seed_registers"])
        assert trials == 10 and 0 < dimension < 2**bits
        fail = (1 - F(dimension, 2**bits))**trials
        failures.append(fail)
        dimensions.append(str(dimension))
    assert dimensions == r5["fixed_program"]["dimensions"]
    p_success = (1 - failures[0]) * (1 - failures[1])
    p_star = F(r5["fixed_program"]["worst_joint_success_probability_lower"])
    assert p_success > p_star == F(998, 1000)
    clock = F(r5["normal_clock"]["epsilon_clock_contract"])
    assert 0 < F(r5["normal_clock"]["epsilon_clock_upper"]) < clock
    assert clock == F(4, 10**15)
    assert r5["normal_clock"]["same_conserved_success_event_at_both_contact_sections"]
    assert r5["normal_clock"]["conserved_seed_event_is_not_a_physical_oracle_gate"]
    assert r5["normal_clock"]["success_conditioned_reference_clock_not_reset"]
    assert r5["normal_clock"]["shape_error_not_used_in_clock_entropy"]
    rows = []
    for v, key in ((1, "t600_lambda1"), (2, "t300_lambda1")):
        shape = F(r5["native_flow_transport"][key]["rectangular_shape_error"])
        assert shape == F(1, 10000 * v)
        # Ideal rectangle/copy on AB is unital; both ideal pulled-back
        # A-excited and K-nonzero effects have trace 12 on dimension 16.
        # The two REAL events are transported separately to the common
        # scalar comparator. They need not agree trial by trial.
        error = clock + shape
        conditional_lower = tail - error
        raw_lower = p_star * conditional_lower
        assert conditional_lower > 0 and raw_lower > 0
        rows.append({
            "velocity_label": v,
            "QQ_duration": 600 // v,
            "factory_parent": "1205 fixed H, M=10^112, common T=60000000000606",
            "events_separately_certified": ["A_excited", "K_nonzero"],
            "clock_isometry_to_same_smooth_program_upper": str(clock),
            "smooth_to_rectangle_isometry_upper": str(shape),
            "sum_probability_transport_upper": str(error),
            "S_conditional_probability_lower_exact": str(conditional_lower),
            "S_conditional_probability_strict_display_lower": str(floor_strict(conditional_lower)),
            "all_raw_trials_probability_lower_using_success_bound": str(raw_lower),
            "all_raw_trials_probability_strict_display_lower": str(floor_strict(raw_lower)),
            "initial_and_final_success_projector": "same conserved original-label S; not computed P_ok",
            "does_not_postselect_actual_trials": True,
        })
    sources[str(Path(__file__).resolve().relative_to(ROOT)).replace("\\", "/")] = sha(Path(__file__))
    return {
        "schema": "after1218_existing_certificate_common_prediction_v1",
        "status": "SYNTHESIS_CHECK_PASS_NOT_NEW_ROUND_OR_ACCEPTANCE",
        "all_checks_passed": True,
        "check_groups": 4,
        "source_sha256": dict(sorted(sources.items())),
        "groups": {
            "accepted_source_hashes": {"pass": True},
            "source_and_ideal_event": {
                "pass": True, "dimension": 16, "ideal_event_rank": 12,
                "exact_joint_bottom_12_spectral_mass": str(tail),
                "ideal_source": "independent cold and hot output marginals conditional on conserved S",
                "factory_inventory_and_source_R_retained": True,
                "cold_and_hot_port_exchange_has_same_joint_spectrum": True,
                "not_exact_Gibbs_and_not_all_factory_faithful": True,
            },
            "factory_success": {
                "pass": True,
                "pool_dimensions": dimensions,
                "seed_bits_per_candidate": [173, 195],
                "candidates_each_pool": 10,
                "joint_success_strict_lower": str(p_star),
                "joint_success_additional_strict_display_lower": str(floor_strict(p_success, 10**12)),
                "computed_flag_not_used_as_conserved_event": True,
            },
            "native_factory_event_transport": {
                "pass": True, "rows": rows,
                "actual_1205_post_preparation_channel_unital_claimed": False,
                "1203_parent_or_error_substituted_for_1205": False,
                "finite_observation_confidence_claimed": False,
                "S_conditional_heat_bounds_promoted_to_unconditional_heat": False,
                "universal_finite_time_third_law_claimed": False,
                "new_hardware_or_permissions": False,
                "additional_preparation_or_readout_noise_included": False,
            },
        },
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    value = build()
    if args.write:
        with OUT.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
    else:
        assert value == load(OUT), "saved results mismatch"
    print(json.dumps({
        "status": "PASS", "all_checks_passed": True,
        "results_match": True, "check_groups": 4,
        "mode": "initial_write" if args.write else "read_only_replay",
    }, sort_keys=True))

if __name__ == "__main__":
    main()

