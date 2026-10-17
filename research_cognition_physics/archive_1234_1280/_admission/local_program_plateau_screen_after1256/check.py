"""Unnumbered conditional local-program plateau screen.
Read old outward interval endpoints only. No ODE/PDE solver is imported.
Default compares without writes; --write creates results.json exclusively once.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent.parent
SOURCES = {
    "1247/results.json": "dbcd590e7816bc423c49a4e1bf83bab6fa50a84dd5578249cc493a15cd72bf06",
    "1247/proof.md": "f3a79f8f43aa235d2fa13cc87a232e2b64643fc1bade375bfee96c204ab4bf3d",
    "1246/proof.md": "9a0247d1b9b43df6d14f8c211a68a36f7251c34a79e544102ea5d85137b40c6e",
}


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def interval(pair):
    result = [F(x) for x in pair]
    assert len(result) == 2 and result[0] <= result[1]
    return result


def subtract(a, b):
    return [a[0]-b[1], a[1]-b[0]]


def encode(obj):
    if isinstance(obj, F):
        return str(obj)
    if isinstance(obj, dict):
        return {key: encode(value) for key, value in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [encode(value) for value in obj]
    return obj


def calculate():
    source_records = []
    for relative, expected in SOURCES.items():
        actual = sha(STAGE / relative)
        assert actual == expected, (relative, "frozen source changed")
        source_records.append({"path_from_stage": relative, "sha256": actual})
    saved = json.loads((STAGE / "1247/results.json").read_text(encoding="utf-8"))
    assert saved["round"] == 1247 and saved["status"] == "PASS"
    fm = saved["group_finite_memory"]
    assert fm["status"] == "PASS"
    hold = F(1,4)
    T = F(9,4)
    R, Rs, period = F(5,2), F(3), F(10)
    target_light_bound, receiver_label_radius = F(5,4), F(5,4)
    assert R < Rs < period/2
    branches, Ks, rods = [], [], []
    for key in ("alpha0", "alpha1"):
        checkpoint = fm[key]["checkpoints"]["2"]
        state = checkpoint["state"]
        u2, L2 = interval(state["u"]), interval(state["L"])
        assert u2[0] > 0 and u2[1] < F(3,8)
        assert L2[0] > 0
        L_upper = L2[1] + hold*u2[1]
        assert L_upper < target_light_bound
        margin = R-receiver_label_radius-L_upper
        assert margin > 0
        K2, rod2 = interval(state["K"]), interval(checkpoint["Brod"])
        Ks.append(K2)
        rods.append(rod2)
        # B_rod = u^2/10; old outward interval must enclose endpoint arithmetic.
        assert rod2[0] <= u2[0]**2/10 and rod2[1] >= u2[1]**2/10
        branches.append({"branch": key, "u2_interval": u2,
                         "L2_interval": L2, "L_T_upper": L_upper,
                         "margin_to_five_quarters": target_light_bound-L_upper,
                         "center_platform_margin": R-L_upper,
                         "receiver_radius": receiver_label_radius,
                         "receiver_platform_margin": margin,
                         "K2_interval": K2, "Brod2_interval": rod2})
    delta_K = subtract(Ks[1], Ks[0])
    assert delta_K == interval(fm["delta_K_final"])
    eta = F(fm["read_error_per_history_K"])
    assert eta == F(1,10**10)
    margin = [x-2*eta for x in delta_K]
    assert margin == interval(fm["margin_after_two_errors"])
    assert margin[0] > F("0.00000006485193048082292") > 0
    delta_rod = subtract(rods[1], rods[0])
    assert delta_rod[0] == F(11204728464402934542638046227,
                           15625000000000000000000000000000000)
    assert delta_rod[0] > F("0.000000717102621721787") > 0
    source_group = saved["group_source_and_retention"]
    assert F(source_group["initial_rho_both"]) == F(111,100)
    assert F(source_group["initial_Mc_both"]) == F(187,200)
    return encode({
        "status": "PASS_WITH_UNPROVED_DEVELOPMENT_PRECONDITIONS",
        "kind": "unnumbered_local_program_plateau_inheritance",
        "arithmetic_groups": 3,
        "new_scientific_groups": 0,
        "new_round": None,
        "adopted_axioms": 0,
        "provenance": {"source_files": source_records,
                       "check_sha256": sha(Path(__file__)),
                       "readme_sha256": sha(HERE/"README.md"),
                       "no_old_integrator_imported": True,
                       "no_new_ODE_or_PDE_integrated": True},
        "platform": {"inner_radius": R, "outer_radius": Rs,
                     "torus_period": period, "injection_radius": period/2,
                     "hold_duration": hold, "task_terminal_clock_age": T,
                     "common_strict_light_upper": target_light_bound,
                     "receiver_label_radius": receiver_label_radius},
        "initial_source_reused": {"M_c": F(187,200), "rho": F(111,100),
                       "local_B_initial_source_identity_is_analytic_reuse": True,
                       "preparation_resources_not_proved_free": True},
        "branches": branches,
        "conditional_record_difference": {
            "K1_minus_K0_interval": delta_K,
            "read_error_per_history": eta,
            "lower_bound_after_two_read_errors": margin[0],
            "old_margin_interval_reused": margin,
            "read_error_contract_requires_explicit_local_task_extension": True},
        "conditional_local_geometry": {
            "quantity": "Brod = k^2*u^2 = u^2/10 for one fixed rod",
            "clock_age": F(2), "difference_interval": delta_rod,
            "actual_geometry_instrument_not_added": True},
        "unproved_preconditions": [
            "LWP-DD for complete Einstein/two-forced-flow/internal-C2-pulse system",
            "same complete initial-data development contains task causal domain",
            "Lorentz/positive-F/positive-mass/no-caustic domain through task",
            "material identity and actual preparation/terminal instrument contract"],
        "scope": {
            "all_nonzero_development_predictions_conditional": True,
            "transition_annulus_not_integrated_or_certified": True,
            "no_direct_patching_of_two_homogeneous_solutions": True,
            "not_independent_remote_sender_task": True,
            "small_amplitude_signal_not_inferred_from_endpoint_gap": True,
            "no_quantum_CC_or_full_GR_completion": True}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="exclusive first save")
    args = parser.parse_args()
    destination = HERE/"results.json"
    if args.write and destination.exists():
        raise FileExistsError("Refusing to replace existing results.json")
    if not args.write and not destination.exists():
        raise FileNotFoundError("No saved result; first save needs --write")
    result = calculate()
    if args.write:
        with destination.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        action = "CREATED_ONCE"
    else:
        assert result == json.loads(destination.read_text(encoding="utf-8")), "saved result mismatch"
        action = "READ_ONLY_MATCH"
    print(json.dumps({"status": result["status"], "action": action,
                      "arithmetic_groups": result["arithmetic_groups"],
                      "light_upper_both": "<5/4",
                      "K_reported_difference_lower": result["conditional_record_difference"]["lower_bound_after_two_read_errors"],
                      "new_scientific_groups": 0}, ensure_ascii=False))


if __name__ == "__main__":
    main()
