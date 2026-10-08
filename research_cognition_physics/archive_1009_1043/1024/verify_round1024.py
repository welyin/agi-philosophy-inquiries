"""Freeze only the 1024 audit delivery and fixed historical inputs, never live navigation."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
from urllib.parse import unquote

import numpy as np
import parent_bridge_audit as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent / "research_note_1024.md"
OUT = HERE / "research_round_1024_checks.json"
OWN = ["parent_bridge_audit.py", "parent_bridge_audit_results.json", "review.md",
       "bridge_ledger.json", "input_dependency_update_v0_13.md", "NEXT.md", "verify_round1024.py"]
HISTORICAL = {f"archive_1009_/research_note_{n}.md" for n in range(1009, 1024)} | {
    "archive_956_989/981/drafts/common_parent_contract_v1.md",
    "archive_990_1008/993/common_candidate_v1.md",
    "archive_1009_/1009/input_dependency_ledger_v0_1.md"}
LINK = re.compile(r"(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)|^\[[^\]\n]+\]:\s*(\S+)\s*$", re.M)
EXCLUDED = re.compile(r"\\\[.*?\\\]|\$\$.*?\$\$|^```.*?^```\s*$", re.S | re.M)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank_mod_prime(matrix, prime=101):
    """An independent exact lower bound for rational matrix rank."""
    a = [[F(x).numerator * pow(F(x).denominator, -1, prime) % prime for x in row] for row in matrix]
    pivot = 0
    for col in range(len(a[0])):
        found = next((i for i in range(pivot, len(a)) if a[i][col]), None)
        if found is None:
            continue
        a[pivot], a[found] = a[found], a[pivot]
        inverse = pow(a[pivot][col], -1, prime)
        a[pivot] = [x * inverse % prime for x in a[pivot]]
        for i in range(pivot + 1, len(a)):
            factor = a[i][col]
            a[i] = [(x - factor * y) % prime for x, y in zip(a[i], a[pivot])]
        pivot += 1
    return pivot


def closed_hessian_coefficients(kappa):
    # Construct from the independently derived closed block Hessian, not differentiation.
    result = {}
    zero = (0,) * 5
    result[zero] = np.diag(list(map(F, [-1, -1, -1, -1, 1])))
    for i in range(4):
        e = [0] * 5; e[i] = 2
        values = [F(1)] * 4 + [F(kappa, 2)]; values[i] += 2
        result[tuple(e)] = np.diag(values)
        for j in range(i + 1, 4):
            e = [0] * 5; e[i] = e[j] = 1
            m = np.full((5, 5), F(0), dtype=object); m[i, j] = m[j, i] = F(2)
            result[tuple(e)] = m
        if kappa:
            e = [0] * 5; e[i] = e[4] = 1
            m = np.full((5, 5), F(0), dtype=object); m[i, 4] = m[4, i] = F(kappa)
            result[tuple(e)] = m
    e = [0] * 5; e[4] = 2
    result[tuple(e)] = np.diag([F(kappa, 2)] * 4 + [F(3)])
    return result


def local_links(text):
    for match in LINK.finditer(EXCLUDED.sub(lambda match: " " * len(match[0]), text)):
        value = match[1] if match[1] is not None else match[2]
        stripped = value.strip()
        target = stripped[1:stripped.index("> ")] if stripped.startswith("<") and "> " in stripped else None
        if target is None:
            target = stripped[1:stripped.index(">")] if stripped.startswith("<") and ">" in stripped else re.split(r"\s+[\"']", stripped, 1)[0]
        local = unquote(target.split("#", 1)[0])
        if local and not re.match(r"^[a-zA-Z]+:", local) and not local.startswith("//"):
            yield target, local


def verify(prospective=False):
    fresh = science.run()
    science.compare(fresh, json.loads(science.OUT.read_text("utf8")))
    assert fresh["round"] == 1024 and fresh["all_audit_certificates_passed"]
    assert fresh["new_cognitive_axioms"] == fresh["new_physical_calibration_groups"] == 0
    assert fresh["new_audit_certificate_groups"] == 1 and fresh["cumulative_research_groups"] == 3801
    assert not fresh["parent_hypotheses_derived_from_cognition"]
    assert not fresh["full_parent_process_certified"]
    assert not fresh["machine_checks_prove_human_applicability_classifications"]

    fields = fresh["chiral_mass_interface"]
    y = np.array([1] * 6 + [-4] * 3 + [2] * 3 + [-3] * 2 + [6])
    em = np.array([2] * 3 + [-1] * 3 + [-2] * 3 + [1] * 3 + [0, -3, 3])
    assert np.array_equal(fields["three_generation_6Y"], np.tile(y, 3))
    assert fields["number_of_left_handed_Weyl_fields"] == 45
    assert sum(y) == sum(y**3) == 0
    assert F(fields["linear_anomaly_exact"]) == F(fields["cubic_anomaly_exact"]) == 0
    assert 0 not in y and not (set(y) & set(-y))
    assert fields["fully_chiral_under_hypercharge"]
    components = fields["dimension_four_graph_components"]
    assert len(components) == 2
    assert [(r["charge_levels"], r["partition_sizes"], r["symmetric_mass_rank_upper_bound"]) for r in components] == [([-4, 1, 2], [18, 18], 36), ([-3, 6], [3, 6], 6)]
    # Independently verify every charge-level support edge and the bipartition bound.
    levels = sorted(set(map(int, y)))
    edges = {(a, b) for a in levels for b in levels if abs(a + b) == 3}
    assert edges == {(-4, 1), (1, -4), (1, 2), (2, 1), (-3, 6), (6, -3)}
    Mh = np.array(fields["H_mass_matrix"], dtype=int)
    Mc = np.array(fields["H_conjugate_mass_matrix"], dtype=int)
    M2 = np.array(fields["HH_Weinberg_mass_matrix"], dtype=int)
    assert Mh.shape == Mc.shape == M2.shape == (45, 45)
    assert np.array_equal(Mh, Mh.T) and np.array_equal(Mc, Mc.T) and np.array_equal(M2, M2.T)
    charged, full = Mh + Mc, Mh + Mc + M2
    assert rank_mod_prime(charged.tolist()) == fields["charged_matrix_exact_rank"] == 42
    assert rank_mod_prime(full.tolist()) == fields["Weinberg_extended_exact_rank"] == 45
    assert fields["dimension_four_maximum_rank"] == 42
    assert np.count_nonzero(M2) == 3 and list(np.flatnonzero(np.diag(M2))) == [12, 27, 42]
    assert len(fields["charged_pairs"]) == 21 and np.count_nonzero(charged) == 42
    for g in range(3):
        o = 15 * g
        assert np.array_equal(charged[o:o+3, o+6:o+9], (3*g+1)*np.eye(3, dtype=int))
        assert np.array_equal(charged[o+3:o+6, o+9:o+12], (3*g+2)*np.eye(3, dtype=int))
        assert charged[o+13, o+14] == charged[o+14, o+13] == 3*g+3
    assert fields["each_quark_mass_shared_across_three_colors"]
    assert fields["color_ward_generator_count"] == len(fields["color_ward_identity_errors"]) == 8
    assert max(fields["color_ward_identity_errors"]) == 0
    for M, rhs in ((Mh, -3), (Mc, 3), (M2, -6)):
        charges = np.tile(y, 3)
        assert np.array_equal((charges[:, None] + charges[None, :]) * M, rhs * M)
    assert np.array_equal(fields["per_generation_3Qem"], em)
    assert fields["electric_zero_charge_multiplicity"] == 3 and fields["electric_Higgs_vev_charge"] == 0
    assert fields["electric_opposite_charge_pairs"] == [1, 2, 3] and not fields["fully_chiral_under_electric_charge"]
    charges = np.tile(em, 3)
    assert np.all((charges[:, None] + charges[None, :]) * full == 0)
    assert max(fields["polynomial_charge_identity_errors"].values()) == 0
    assert max(fields["finite_phase_covariance_errors"]) < 1e-11
    assert fields["round1022_full_rank_dimension_four_hypothesis_failed"]
    assert fields["Weinberg_extension_changes_operator_dimension"] and fields["two_distinct_U1_generators_not_interchangeable"]

    hessians = fresh["Higgs_portal_commutant"]
    assert len(hessians) == 2
    ranks = []
    for row in hessians:
        kappa = int(F(row["kappa_exact"]))
        expected = closed_hessian_coefficients(kappa)
        assert sorted(expected) == list(map(tuple, row["monomial_exponents"]))
        assert len(expected) == row["coefficient_monomial_count"] == (16 if kappa else 12)
        for e, values in zip(row["monomial_exponents"], row["coefficient_matrices"]):
            assert np.array_equal(expected[tuple(e)], np.array([[F(x) for x in line] for line in values], dtype=object))
        elementary = []
        for i in range(5):
            for j in range(i, 5):
                q = np.full((5, 5), F(0), dtype=object); q[i, j] = q[j, i] = F(1)
                elementary.append(q)
        system = []
        for h in expected.values():
            matrices = [q @ h - h @ q for q in elementary]
            for i in range(5):
                for j in range(i + 1, 5):
                    system.append([m[i, j] for m in matrices])
        rank = rank_mod_prime(system)
        assert rank == row["exact_constraint_rank"] == (14 if kappa else 13)
        nulls = [np.eye(5, dtype=int)] if kappa else [np.diag([1, 1, 1, 1, 0]), np.diag([0, 0, 0, 0, 1])]
        # Known independent rational null vectors give the matching rank upper bound.
        assert all(np.array_equal(q @ h, h @ q) for q in nulls for h in expected.values())
        assert 15 - rank == row["exact_commutant_dimension"] == len(nulls)
        saved_nulls = [np.array([[F(x) for x in line] for line in q], dtype=object) for q in row["nullspace_matrices"]]
        assert rank_mod_prime([q.ravel().tolist() for q in saved_nulls]) == len(nulls)
        assert all(np.array_equal(q @ h, h @ q) for q in saved_nulls for h in expected.values())
        assert max(row["exact_closed_formula_evaluation_errors"]) == 0
        assert row["all_polynomial_coefficients_tested"] and row["potential_algebra_only_not_full_gauge_dynamics"]
        ranks.append(rank)

    trunc = fresh["gauge_truncation"]
    assert trunc["signature"] == "-+++" and F(trunc["hypercharge_y_exact"]) == F(1, 2)
    assert F(trunc["exact_source_coefficient_over_g"]) == 1
    assert len(trunc["time_samples"]) == 4
    for row in trunc["time_samples"]:
        assert row["scalar_equation_residual"] < 1e-12 and row["sigma_equation_residual"] == 0
        assert abs(row["hypercharge_action_source_over_g"] - 1) < 1e-12
        assert row["source_imaginary_error"] == row["Maxwell_left_hand_divergence_at_A_zero"] == 0
        assert abs(row["Maxwell_Euler_residual_over_g"] - 1) < 1e-12
    assert trunc["scalar_only_equations_satisfied"] and not trunc["full_gauge_equations_satisfied"]
    assert trunc["no_failure_of_full_parent_dynamics_inferred"] and trunc["no_universal_gravity_coupling_conclusion_inferred"]

    ledger = json.loads((HERE / "bridge_ledger.json").read_text("utf-8-sig"))
    rows = ledger["rows"]
    assert ledger["round"] == 1024 and len(rows) == len({r["id"] for r in rows}) == 17
    statuses = dict(Counter(r["status"] for r in rows))
    assert statuses == {"adopted_identity": 3, "mapped_interface": 2, "conditional_tool": 8, "missing_bridge": 4}
    assert {n for row in rows for n in row["source_rounds"]} == set(range(1009, 1024))
    assert ledger["new_cognitive_axioms"] == ledger["new_physical_calibration_groups"] == 0
    assert ledger["new_audit_certificate_groups"] == 1 and ledger["cumulative_research_groups"] == 3801
    ledger_sources = set(ledger["parent_specs"])
    for row in rows:
        assert row["cognitive_origin_closed"] is False and row["full_parent_process_certified"] is False
        assert all(isinstance(row[key], str) and row[key] for key in ("object", "claim", "requires", "not_proved"))
        ledger_sources.update(row["source_paths"])
    assert ledger_sources <= HISTORICAL
    assert set(fresh["historical_source_sha256"]) == HISTORICAL
    for name, digest in fresh["historical_source_sha256"].items():
        assert sha(BASE / name) == digest, name
    assets = [NOTE] + [HERE / name for name in OWN]
    assert len(assets) == 8 and all(p.name not in {"README.md", "RESEARCH_STATE.md", "research_direction.md"} for p in assets)
    count = 0
    for path in assets:
        assert path.is_file(), path
        if path.suffix == ".md":
            text = path.read_text("utf-8-sig"); assert text.count("$$") % 2 == 0, path
            for target, local in local_links(text):
                resolved = (path.parent / local.replace("\\", "/")).resolve()
                assert resolved.exists() or (prospective and resolved == OUT), (path, target)
                count += 1
    assert all(f"## {i}." in NOTE.read_text("utf8") for i in range(1, 11))
    assert "独立科学签审通过" in (HERE / "review.md").read_text("utf8")
    return dict(round=1024, date="2026-10-08", all_delivery_checks_passed=True,
                audit_result_reproduced=True, new_physical_calibration_groups=0,
                new_audit_certificate_groups=1, cumulative_research_groups=3801,
                local_links_checked=count, frozen_current_files=8, historical_input_files=len(HISTORICAL),
                live_navigation_frozen=False, neighboring_round_frozen=False,
                ledger_rows=17, ledger_status_counts=statuses, ledger_historical_sources=len(ledger_sources),
                human_semantic_classifications_not_machine_proved=True,
                chiral_field_count=45, maximum_dimension_four_rank=42, Weinberg_extended_rank=45,
                color_Ward_generator_checks=8, exact_color_delta_blocks_checked=True,
                Hessian_constraint_ranks=ranks, Hessian_commutant_dimensions=[1, 2],
                rational_rank_lower_bounds_independently_checked_modulo_prime=101,
                scalar_solution_nonzero_gauge_source_exact="1", all_parent_process_flags_false=True,
                new_physical_prediction=False, parent_hypotheses_derived_from_cognition=False,
                goal_completed=False, visual_checks_performed=False,
                maximum_mass_covariance_error=max(fields["finite_phase_covariance_errors"]),
                historical_source_sha256=fresh["historical_source_sha256"],
                source_sha256={str(path.relative_to(ROOT)).replace("\\", "/"): sha(path) for path in assets})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--write", action="store_true")
    args = parser.parse_args(); result = verify(prospective=args.write)
    if args.write:
        with OUT.open("x", encoding="utf8") as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    else:
        assert result == json.loads(OUT.read_text("utf8"))
    print(json.dumps({key: value for key, value in result.items() if not key.endswith("sha256")}, ensure_ascii=False, indent=2))
