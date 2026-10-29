"""Count-zero stage auxiliary verification; not a new scientific round.

Checks frozen 1323--1326 receipts and hashes, paper links/math conventions,
the unchanged goal, and analytic parameter transport from frozen results.
No scientific matrix is rescanned or old implementation imported. A PASS is
not a declaration of goal completion or acceptance of new cognitive axioms.
Default/--check are read-only; --write-new exclusively creates the result.
"""
from __future__ import annotations
import argparse
import cmath
import hashlib
import json
import math
import re
from pathlib import Path
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
RESEARCH = STAGE.parent
ROOT = RESEARCH.parent
PAPER = RESEARCH / "16_未知参考交接的作用选择与动态来源_阶段论文.md"
EXPECTED_GOAL = "95bbeb961464c39cdd4e3867d6f00faec868f38fd5cb2cf43c59d4393b7e3eb9"
RESULT = HERE / "stage_joint_results.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path):
    return path.relative_to(ROOT).as_posix()


def same_result(a, b, name="result"):
    if isinstance(a, dict):
        assert isinstance(b, dict) and a.keys() == b.keys(), name
        for key in a:
            same_result(a[key], b[key], name + "." + key)
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), name
        for i, (x, y) in enumerate(zip(a, b)):
            same_result(x, y, f"{name}[{i}]")
    elif isinstance(a, float):
        assert isinstance(b, (int, float)) and not isinstance(b, bool), name
        assert math.isfinite(b) and math.isclose(a, b, rel_tol=2e-12, abs_tol=2e-12), name
    else:
        assert type(a) is type(b) and a == b, name


def markdown_checks(path):
    text = path.read_text(encoding="utf-8")
    assert not any(ord(c) < 32 and c not in "\r\n\t" for c in text), "control character"
    # Exclude fenced code examples from rendering-convention checks.
    prose = re.sub(r"(?ms)^\s*(```|~~~).*?^\s*\1\s*$", "", text)
    display_lines = [i + 1 for i, line in enumerate(prose.splitlines()) if line.strip() == "$$"]
    assert len(display_lines) % 2 == 0 and display_lines, "unpaired/missing display fences"
    assert prose.count("$$") == len(display_lines), "display $$ must stand on its own line"
    assert not any(old in prose for old in (r"\[", r"\]", r"\(", r"\)")), "legacy delimiters"
    without_display = re.sub(r"(?s)\$\$.*?\$\$", "", prose)
    assert len(re.findall(r"(?<!\\)\$", without_display)) % 2 == 0, "unpaired inline math"
    local_links = []
    external_count = 0
    # All links in the manuscript use ordinary inline Markdown, optionally
    # enclosing a space-bearing target in angle brackets.
    for match in re.finditer(r"!?\[[^\]\n]*\]\((<[^>]+>|[^)\n]+)\)", text):
        target = match.group(1).strip()
        if target.startswith("<") and target.endswith(">"):
            target = target[1:-1]
        else:
            target = re.sub(r'\s+"[^"\n]*"$', "", target)
        if re.match(r"^(https?|mailto|app|codex):", target, flags=re.I):
            external_count += 1
            continue
        target_path, _, anchor = target.partition("#")
        if not target_path:
            continue
        target_path = unquote(re.sub(r":\d+$", "", target_path))
        if re.match(r"^/[A-Za-z]:[/\\]", target_path):
            target_path = target_path[1:]
        linked = Path(target_path)
        if not linked.is_absolute():
            linked = path.parent / linked
        linked = linked.resolve()
        # The sole forward target allowed is this script's own exclusive
        # output; main verifies it exists after publication and on replay.
        assert linked.exists() or linked == RESULT.resolve(), ("missing paper link", target)
        local_links.append(dict(target=target, resolved=rel(linked) if linked.is_relative_to(ROOT)
                                else str(linked), anchor_present=bool(anchor),
                                own_result_publication_target=linked == RESULT.resolve()))
    assert local_links, "no local evidence links found"
    return dict(paper=rel(path), sha256=digest(path), local_links=local_links,
                external_links_not_refetched=external_count,
                display_equations=len(display_lines)//2,
                inline_dollar_count=len(re.findall(r"(?<!\\)\$", without_display)),
                scope="syntax, control characters and local target existence; not image rendering")


def compute():
    bindings = {}
    checks = []

    def bind(path, expected=None):
        actual = digest(path)
        if expected is not None:
            assert actual == expected, ("hash mismatch", str(path), actual, expected)
        key = rel(path)
        if key in bindings:
            assert bindings[key] == actual
        bindings[key] = actual
        return actual

    def close(name, a, b, tolerance=2e-11):
        residual = abs(a-b)
        assert residual <= tolerance, (name, residual)
        checks.append(dict(name=name, residual=float(residual), tolerance=tolerance))

    bind(HERE / "goal.md", EXPECTED_GOAL)
    receipt = read_json(HERE / "goal_start_receipt.json")
    assert receipt["goal_document_sha256"] == EXPECTED_GOAL
    bind(HERE / "goal_start_receipt.json")
    bind(HERE / "goal_scope_review.md", receipt["scope_review_sha256"])
    rounds = []
    for round_id in range(1323, 1327):
        folder = STAGE / str(round_id)
        acceptance_path = folder / "acceptance.json"
        acceptance = read_json(acceptance_path)
        bind(acceptance_path)
        assert acceptance["round"] == round_id
        assert acceptance["status"] == "ACCEPTED_CONDITIONAL_SCOPE"
        assert acceptance["new_axioms"] == 0
        assert acceptance["scientific_evidence_files"] == 10
        expected_paths = {rel(folder / name) for name in
                          ("proof.md", "check.py", "results.json", "NEXT.md",
                           "independent_check.py", "independent_results.json",
                           "independent_review.md", "scope_review.md")}
        expected_paths |= {rel(STAGE / f"research_note_{round_id}.md"),
                           rel(STAGE / f"_admission/round{round_id}_admission.md")}
        assert set(acceptance["sha256"]) == expected_paths
        for path, wanted in acceptance["sha256"].items():
            bind(ROOT / path, wanted)
        author = read_json(folder / "results.json")
        independent = read_json(folder / "independent_results.json")
        assert author["round"] == round_id and author["status"] == "PASS"
        assert author["count"] == len(author["checks"])
        assert independent["status"] == "PASS_INDEPENDENT_CHECKS_ONLY"
        assert independent["author_code_imported"] is False
        independent_count = independent["checks"]
        assert isinstance(independent_count, int) and independent_count > 0
        for name, wanted in author["source_sha256"].items():
            bind(ROOT / name, wanted)
        if "source_bindings" in independent:
            for name, wanted in independent["source_bindings"].items():
                bind(STAGE / name, wanted)
        else:
            for name, wanted in independent["source_sha256"].items():
                bind((folder / name).resolve(), wanted)
        executed = acceptance["executed_checks"]
        assert len(executed) == 2
        for item in executed:
            assert item["exit_code"] == 0
            bind(ROOT / item["script"], item["script_sha256"])
            expected_count = independent_count if Path(item["script"]).name == "independent_check.py" else author["count"]
            assert item["checks"] == expected_count
        rounds.append(dict(round=round_id, acceptance_sha256=digest(acceptance_path),
                           scientific_files=10, author_checks=author["count"],
                           independent_checks=independent_count,
                           status=acceptance["status"]))

    paper_checks = markdown_checks(PAPER)
    bind(PAPER)
    bind(Path(__file__).resolve())
    result1325 = read_json(STAGE / "1325/results.json")
    result1326 = read_json(STAGE / "1326/results.json")

    # Same fixed scalar portals as 1326: no coupling re-fit. This is a
    # parameter-setting bridge, not a fresh dynamical matrix experiment.
    gd, gu = result1326["model"]["normalized_couplings"]
    unit = .001
    ratios = [gd/unit, gu/unit]
    close("scalar_D_odd_integer", ratios[0], 13.)
    close("scalar_U_odd_integer", ratios[1], 11.)
    transfer_time = math.pi/(2*unit)
    close("scalar_common_clock", transfer_time, 500*math.pi)
    scalar_probabilities = [math.sin(g*transfer_time)**2 for g in (gd, gu)]
    for i, value in enumerate(scalar_probabilities):
        close(f"scalar_port_{i}_full_transfer", value, 1.)
    cutoff = result1326["model"]["boson_cutoff"]
    free_phase_error = max(abs(cmath.exp(-1j*K*transfer_time)-1) for K in range(cutoff+1))
    close("scalar_complete_safe_K_free_phase", free_phase_error, 0., 3e-11)
    scalar_bridge = dict(portals=[gd, gu], odd_multipliers=ratios,
        common_time=transfer_time, transfer_probabilities=scalar_probabilities,
        safe_K_max=cutoff, free_phase_error=free_phase_error,
        transfer_phases=[[float((-1j*math.sin(g*transfer_time)).real),
                          float((-1j*math.sin(g*transfer_time)).imag)] for g in (gd, gu)],
        old_setting="f=p=0; E_D=0,E_L=1,E_U=2, omega=1; H0=K",
        driven_setting="E_D=E_L=E_U=0; p=.005; f=1 or 2; same portals and reference carrier",
        switching_work="Tr rho [f S_f + D(p) - n_L - 2 n_U]; changed preparation separately registered",
        scope="analytic transport of the old full-input theorem, not a new scientific round")

    budget = result1325["reference_budget"]
    scale = .1
    scaled_beta = scale*budget["beta"]
    scaled_tau = budget["tau"]/scale
    scaled_time = budget["common_full_transfer_time"]/scale
    close("three_flavor_scaled_clock", scaled_time, transfer_time)
    close("scaled_V_tau_coefficient", scale*scaled_tau, budget["tau"])
    close("scaled_V_T_coefficient", scale*scaled_time, budget["common_full_transfer_time"])
    original_spectra = {"D": [.13, .03, .01], "U": [.11, .05, .01]}
    scaled_spectra = {side: [scale*x for x in values] for side, values in original_spectra.items()}
    for side, values in scaled_spectra.items():
        for i, value in enumerate(values):
            close(f"three_flavor_{side}_{i}_odd_transfer", math.sin(value*scaled_time)**2, 1.)
    rows = []
    for case in budget["cases"]:
        eta = scale*case["eta"]
        curvature = scale*scale*case["curvature"]
        rate = scale*case["finite_specified_rate"]
        allowance = scale*case["taylor_rate_allowance"]
        close(case["label"]+".scaled_curvature", curvature, eta*eta)
        close(case["label"]+".same_probability_from_scaled_rate",
              (rate*scaled_tau)**2, case["reference_exit_probability"])
        passed = eta <= scaled_beta
        assert passed == case["passes_curvature_budget"]
        if passed:
            assert eta+allowance < scaled_beta
        else:
            assert eta-allowance > scaled_beta
        rows.append(dict(label=case["label"], eta=eta, curvature=curvature,
            budget_passed=passed, finite_rate=rate, taylor_rate_allowance=allowance,
            guaranteed_rate_interval=[eta-allowance, eta+allowance],
            transported_reference_exit_probability=case["reference_exit_probability"],
            transportation="exact e^(-i (scale V)(t/scale)) = e^(-i V t); no matrix rescanned"))
    # In the physical picture the short-time free phase is not unchanged.
    short_free_phase_difference = abs(cmath.exp(-1j*scaled_tau)-cmath.exp(-1j*budget["tau"]))
    assert short_free_phase_difference > .1
    long_original = max(abs(cmath.exp(-1j*K*budget["common_full_transfer_time"])-1)
                        for K in range(cutoff+1))
    close("original_long_clock_free_phase", long_original, 0., 3e-11)
    scaling = dict(scale=scale, original_spectra=original_spectra, scaled_spectra=scaled_spectra,
        original_beta=budget["beta"], scaled_beta=scaled_beta,
        original_tau=budget["tau"], scaled_tau=scaled_tau,
        original_transfer_time=budget["common_full_transfer_time"], scaled_transfer_time=scaled_time,
        cases=rows, short_time_physical_K1_free_phase_difference=short_free_phase_difference,
        output_identity_scope="complete interaction-picture reference/matter/partner output; physical-picture readout co-rotates with H0",
        long_time_scope="both old and scaled long clocks have e^(-i K T)=1 for all safe integer K",
        independent_parameter_points="the scaled three-flavor budget matrices are not the single-flavor 1326 low-band table",
        unchanged_initial_data="same reference, phase, unknown input and partner, record/readout and full failed branches in the interaction picture",
        source_scope="V and its interaction-energy contribution scale; H0 remains registered; no assertion that all physical source numbers stay identical")
    return dict(status="PASS_STAGE_AUXILIARY_CHECKS_ONLY", new_scientific_rounds=0,
        new_axioms=0, stage_completion_declared=False,
        goal_completion_not_authorized_by_this_script=True,
        goal_sha256=EXPECTED_GOAL, current_application_goal_status="not queried",
        round_evidence=rounds, scientific_files_verified=sum(r["scientific_files"] for r in rounds),
        previously_recorded_author_checks=sum(r["author_checks"] for r in rounds),
        previously_recorded_independent_checks=sum(r["independent_checks"] for r in rounds),
        checks_are_hash_and_metadata_verification_not_reruns=True,
        paper_checks=paper_checks, scalar_same_portal_bridge=scalar_bridge,
        three_flavor_scaling_transport=scaling,
        analytic_arithmetic_checks=checks,
        max_transport_residual=max(x["residual"] for x in checks),
        source_sha256=dict(sorted(bindings.items())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = compute()
    if args.write_new:
        with RESULT.open("x", encoding="utf-8") as output:
            json.dump(result, output, ensure_ascii=False, indent=2)
            output.write("\n")
    else:
        same_result(result, read_json(RESULT))
    for link in result["paper_checks"]["local_links"]:
        linked = Path(link["resolved"])
        if not linked.is_absolute():
            linked = ROOT/linked
        assert linked.exists(), ("unfulfilled local link after result publication", str(linked))
    print(json.dumps(dict(status=result["status"], mode="write-new" if args.write_new else "read-only",
        scientific_files_verified=result["scientific_files_verified"],
        paper_local_links=len(result["paper_checks"]["local_links"]),
        analytic_arithmetic_checks=len(result["analytic_arithmetic_checks"]),
        max_transport_residual=result["max_transport_residual"],
        stage_completion_declared=False)))


if __name__ == "__main__":
    main()
