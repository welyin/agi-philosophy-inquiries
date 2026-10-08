"""Limited round 1019 delivery verification; no live navigation is frozen."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
from urllib.parse import unquote

import numpy as np
import spin_statistics_locality as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent / "research_note_1019.md"
OUT = HERE / "research_round_1019_checks.json"
OWN = ["selection_audit.md", "spin_statistics_locality.py",
       "spin_statistics_locality_results.json", "review.md", "verify_round1019.py",
       "input_dependency_update_v0_8.md", "NEXT.md"]
HISTORICAL = {
    "archive_1009_/research_note_1018.md",
    "archive_1009_/1018/gravity_ideal_selection_results.json",
    "archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md",
    "archive_1009_/1009/input_dependency_ledger_v0_1.md",
}
LINK = re.compile(r"(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)|^\[[^\]\n]+\]:\s*(\S+)\s*$", re.M)
EXCLUDED = re.compile(r"\\\[.*?\\\]|\$\$.*?\$\$|^```.*?^```\s*$", re.S | re.M)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a, b):
    assert math.isclose(a, b, rel_tol=2e-11, abs_tol=2e-14), (a, b)


def local_links(text):
    masked = EXCLUDED.sub(lambda match: " "*len(match[0]), text)
    for match in LINK.finditer(masked):
        value = match[1] if match[1] is not None else match[2]
        stripped = value.strip()
        if stripped.startswith("<") and ">" in stripped:
            target = stripped[1:stripped.index(">")]
        else:
            target = re.split(r"\s+[\"']", stripped, 1)[0]
        local = unquote(target.split("#", 1)[0])
        if not local or re.match(r"^[a-zA-Z]+:", local) or local.startswith("//"):
            continue
        yield target, local


def verify(prospective=False):
    fresh = science.run()
    science.compare(fresh, json.loads((HERE / "spin_statistics_locality_results.json").read_text("utf8")))
    assert fresh["round"] == 1019 and fresh["all_scientific_calibrations_passed"]
    assert fresh["new_calibration_groups"] == 1 and fresh["cumulative_test_groups"] == 3797
    assert fresh["new_cognitive_axioms"] == 0
    for name in ("microcausality_generated_from_FUCP", "spin_statistics_complete_theorem_numerically_proved",
                 "QED_local_gauge_invariance_claimed", "physical_device_implementation_claimed",
                 "failure_only_at_infinite_energy_claimed", "norm_cutoff_equals_physical_minimum_scale",
                 "quadrature_errors_certified"):
        assert fresh[name] is False, name
    for name in ("no_signaling_failure_conditional_on_declared_local_operation_menu",
                 "tail_bounds_rigorous_analytic", "whole_state_energy_moments_include_subspace_leakage",
                 "preparation_of_entangled_witness_is_explicit_extra_menu_assumption"):
        assert fresh[name] is True, name
    smearing = fresh["smearing"]
    assert smearing["a"] == 1 and smearing["tau"] == .25
    assert smearing["spatial_support_radius"] == 2 and smearing["total_integral"] == 1
    assert smearing["distribution_test_smoothing_limit_requires_operator_norm_continuity"]
    examples = fresh["separated_support_examples"]
    assert [row["separation"] for row in examples] == [6, 8]
    max_energy_second = 0.
    for row in examples:
        L = row["separation"]
        assert row["strict_spacelike_support_gap"] == L-4-.5 > 0
        assert row["quadrature_refinement_is_not_rigorous_error_certificate"]
        qs = row["quadrature"]
        assert [q["cutoff"] for q in qs] == [64., 128., 256.]
        assert [q["gauss_order"] for q in qs] == [16, 32, 32]
        for q in qs:
            assert q["quadrature_error_certified"] is False
            assert len(q["moments"]) == len(q["cross_moments"]) == 3
            assert all(x > 0 for x in q["moments"])
            for k, tail in enumerate(q["analytic_absolute_tail_bounds"]):
                close(tail, 324/(math.pi**2*(6-k)*q["cutoff"]**(6-k)))
        assert row["quadrature_refinement_differences"][-1]["maximum_moment_difference"] < 1e-9

        b = row["analytical_bounds"]
        # Independently reconstruct all decisive bounds from rational input,
        # with no use of a quadrature result or fitted overlap.
        nl = float(F(9,10)**4*F(383,384)**4)/(8*math.pi**2)
        nu = 1/math.pi**2
        sl = 1/(4*math.pi**2*(L+4)**2)
        su = 1/(4*math.pi**2*((L-4)**2-F(1,4)))
        rl, ru = sl/nu, su/nl
        assert 0 < rl < ru < 1
        dl = rl*math.sqrt(1-ru**2)/math.sqrt(2)
        theta = dl/4
        lower = dl**2/8
        for key, value in dict(n_lower=nl, n_upper=nu, s_lower=sl, s_upper=su,
                               r_lower=rl, r_upper=ru, d_lower=dl, fixed_theta=theta,
                               exact_continuum_probability_lower_bound=lower).items():
            close(b[key], value)
        assert b["inequalities_derived_without_numerical_quadrature"]
        cutoff = b["cutoff_control"]
        assert cutoff["cutoff"] == 2000
        tail = 54/(math.pi**2*2000**6)
        lost = tail/nl
        error = 2*(1+abs(theta))*math.sqrt(lost)
        close(cutoff["normalization_tail_bound"], tail)
        close(cutoff["normalized_lost_norm_squared_upper_bound"], lost)
        close(cutoff["projection_operator_norm_error_upper_bound"], math.sqrt(lost))
        close(cutoff["same_state_delta_probability_error_upper_bound"], error)
        assert error < lower/2 and cutoff["error_less_than_half_signal"]
        assert cutoff["direct_quadrature_to_cutoff_performed"] is False
        assert cutoff["cutoff_operations_strictly_local"] is False
        assert cutoff["all_states_assumed_bandlimited"] is False

        car = row["CAR_calibration"]
        assert car["matrix_dimension"] == 16 and car["number_of_CAR_modes"] == 4
        for key in ("maximum_CAR_error", "projector_error", "charge_zero_commutator_error", "even_parity_commutator_error"):
            assert car[key] < 3e-14
        r = car["overlap"]
        assert rl < r < ru
        d = r*math.sqrt(1-r*r)/math.sqrt(2)
        close(car["vacuum_commutator_norm_squared"], d*d)
        close(car["commutator_operator_norm"], math.sqrt(2)*d)
        close(car["computed_d"], d)
        for label, expected_theta in (("adapted_theta_witness", d/4),
                                      ("fixed_analytic_theta_witness", theta)):
            witness = car[label]
            close(witness["theta"], expected_theta)
            assert witness["postselection_used"] is False
            assert 0 <= witness["initial_probability"] <= 1
            assert 0 <= witness["final_probability"] <= 1
            close(witness["probability_difference"], witness["final_probability"]-witness["initial_probability"])
            assert witness["measured_Taylor_remainder"] <= 2*expected_theta**2+3e-15
            assert witness["absolute_difference"] >= d*expected_theta/2-3e-15
        assert car["fixed_analytic_theta_witness"]["absolute_difference"] >= lower
        energy = car["finite_energy"]
        h1 = np.array(energy["one_particle_H_compression"])
        h2 = np.array(energy["one_particle_H_squared_compression"])
        leakage = np.array(energy["leakage_matrix"])
        assert np.linalg.norm(h2-h1@h1-leakage) < 2e-13
        assert min(np.linalg.eigvalsh(h1)) > 0 and min(np.linalg.eigvalsh(leakage)) > 0
        assert energy["two_wavepacket_space_energy_invariant"] is False
        assert energy["local_Klein_Gordon_energy_density_assumed_positive"] is False
        assert energy["finite_moments_are_analytic_tail_consequences"]
        assert energy["exact_energy_values_are_quadrature_calibrations"]
        assert len(energy["states"]) == 8
        for state in energy["states"]:
            assert state["mean_energy"] >= 0 and math.isfinite(state["energy_second_moment"])
            assert state["energy_second_moment"] >= state["mean_energy"]**2-2e-12
            max_energy_second = max(max_energy_second, state["energy_second_moment"])

    cutoff_samples = fresh["cutoff_projection_algebra"]["samples"]
    assert len(cutoff_samples) == 4
    for sample in cutoff_samples:
        close(sample["projection_difference_norm"], math.sqrt(sample["lost_norm_squared"]))
    assert set(fresh["historical_source_sha256"]) == HISTORICAL
    for name, digest in fresh["historical_source_sha256"].items():
        assert sha(BASE/name) == digest, name

    assets = [NOTE]+[HERE/name for name in OWN]
    assert len(assets) == 8
    assert all(path.name not in {"README.md", "RESEARCH_STATE.md", "research_direction.md"} for path in assets)
    count = 0
    for path in assets:
        assert path.is_file(), path
        if path.suffix == ".md":
            text = path.read_text("utf-8-sig")
            assert text.count("$$") % 2 == 0, path
            for target, local in local_links(text):
                resolved = (path.parent/local.replace("\\", "/")).resolve()
                assert resolved.exists() or (prospective and resolved == OUT), (path, target)
                count += 1
    note = NOTE.read_text("utf8")
    assert all(f"## {i}." in note for i in range(1,11))
    assert "独立科学签审通过" in (HERE/"review.md").read_text("utf8")
    return dict(round=1019, date="2026-10-08", all_delivery_checks_passed=True,
        scientific_result_reproduced=True, new_calibration_groups=1,
        cumulative_research_groups=3797, local_links_checked=count,
        frozen_current_files=len(assets), historical_input_files=len(HISTORICAL),
        live_navigation_frozen=False, neighboring_round_frozen=False,
        spacelike_support_examples=2, CAR_matrix_dimension=16,
        positive_analytic_signal_bounds_independent_of_quadrature=True,
        finite_cutoff_error_below_half_signal=True,
        quadrature_errors_certified=False, strict_tails_analytic=True,
        cutoff_operations_strictly_local=False, finite_energy_states_checked=16,
        maximum_calibrated_energy_second_moment=max_energy_second,
        energy_subspace_leakage_preserved=True,
        local_operation_menu_explicit_input=True, local_QED_gauge_claimed=False,
        physical_device_implementation_claimed=False,
        complete_spin_statistics_theorem_numerically_proved=False,
        new_cognitive_axiom=False, goal_completed=False, visual_checks_performed=False,
        analytic_scope=fresh["scope"],
        historical_source_sha256=fresh["historical_source_sha256"],
        source_sha256={str(path.relative_to(ROOT)).replace("\\", "/"):sha(path) for path in assets})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--write", action="store_true")
    args = parser.parse_args(); result = verify(prospective=args.write)
    if args.write:
        with OUT.open("x", encoding="utf8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False); stream.write("\n")
    else:
        assert result == json.loads(OUT.read_text("utf8"))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith("sha256")}, ensure_ascii=False, indent=2))
