"""Verify final evidence without rewriting historical progress or navigation.

Default: --manifest final_evidence_manifest.json (relative to this directory).
Default execution prints a fresh result and writes nothing. --write-results
exclusively creates stage_completion_checks.json; an existing receipt is never
overwritten. This tool never changes an application goal or replays the old
progress verifier against the now-current management files.

Manifest contract (all embedded paths are project-relative):
  expected_sha256: {path: sha256}; paper_path; review_paths (nonempty list);
  audit_json_path; required_audit_status='PASS_IN_DECLARED_SCOPE';
  latest_progress_path; latest_progress_sha256;
  round_range: [1219,1233] or {first:1219,last:1233}; expected_groups=63;
  expected_status='completed_in_declared_scope'.
Optional: final_goal_receipt_path and final_goal_receipt_sha256 (together),
  required_notice_fragments (list of exact strings, if desired).

The hash map must include the final paper, reviews, audit, common-domain draft,
this verifier, all acceptance/scientific bindings, and current management files.
It must not include the manifest itself or the output receipt (no hash cycle).
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
BASE = ROOT / "research_cognition_physics"
SELF = Path(__file__).resolve()
OUTPUT = HERE / "stage_completion_checks.json"
PUBLISHER = HERE / "publish_round.py"
LEDGER = HERE / "dependency_ledger.json"
NOTICE = HERE / "navigation_notice.md"
LAYER = BASE / "_migration/thermal_geometry_active_1219/layer.json"
RUNNER = ROOT / "scripts/run_research_thermal_geometry_active.py"
UPDATER = HERE / "update_navigation.py"
COMMON = HERE / "common_domain_and_remaining_freedoms_draft.md"
NAV_NAMES = ("README.md", "research_direction.md", "RESEARCH_STATE.md", "ROADMAP.md")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rel(path):
    return path.resolve().relative_to(ROOT).as_posix()


def project_path(value):
    require(isinstance(value, str) and value, "Empty/non-string project path")
    path = Path(value)
    require(not path.is_absolute() and not path.drive, "Expected project-relative path: " + value)
    resolved = (ROOT / path).resolve()
    require(resolved.is_relative_to(ROOT), "Path escapes project: " + value)
    require(resolved.is_file(), "Missing evidence: " + value)
    return resolved


def checked_map(mapping):
    require(isinstance(mapping, dict) and bool(mapping), "Empty SHA map")
    canonical = {}
    for name, digest in mapping.items():
        path = project_path(name)
        require(isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
                "Invalid SHA256: " + name)
        key = rel(path)
        require(key not in canonical, "Duplicate canonical evidence path: " + key)
        require(sha(path) == digest, "SHA256 mismatch: " + key)
        canonical[key] = digest
    return canonical


def require_bound(path, bindings):
    key = rel(path)
    require(key in bindings, "Final manifest omits required evidence: " + key)
    require(sha(path) == bindings[key], "Changed evidence: " + key)
    return key


def run_json(script, *arguments):
    result = subprocess.run(
        [sys.executable, "-B", "-X", "utf8", str(script), *arguments],
        cwd=ROOT, capture_output=True, text=True, encoding="utf-8", check=False)
    require(result.returncode == 0,
            "Read-only check failed: " + str(script) + "\n" + result.stdout + result.stderr)
    return json.loads(result.stdout)


def load_publisher():
    spec = importlib.util.spec_from_file_location("thermal_stage_publication_for_completion", PUBLISHER)
    require(spec is not None and spec.loader is not None, "Cannot load frozen publisher")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify(manifest_path):
    require(sys.flags.optimize == 0, "Run without -O: frozen publisher uses assertion checks")
    manifest_digest = sha(manifest_path)
    manifest = read(manifest_path)
    require(manifest.get("expected_status") == "completed_in_declared_scope", "Wrong stage status contract")
    require(manifest.get("required_audit_status") == "PASS_IN_DECLARED_SCOPE", "Wrong audit status contract")
    require(manifest.get("expected_groups") == 63, "Expected exactly 63 author groups")
    span = manifest.get("round_range")
    if isinstance(span, dict):
        span = [span.get("first"), span.get("last")]
    require(span == [1219, 1233], "Expected rounds 1219 through 1233")
    rounds = list(range(span[0], span[1]+1))
    bindings = checked_map(manifest.get("expected_sha256"))
    require(rel(manifest_path) not in bindings and rel(OUTPUT) not in bindings,
            "Manifest/output hash cycle is not allowed")

    paper = project_path(manifest["paper_path"])
    reviews_raw = manifest.get("review_paths")
    require(isinstance(reviews_raw, list) and bool(reviews_raw), "Stage reviews must be explicit")
    reviews = [project_path(name) for name in reviews_raw]
    require(len(set(reviews)) == len(reviews), "Duplicate stage review")
    audit_path = project_path(manifest["audit_json_path"])
    require(paper != audit_path and paper not in reviews and audit_path not in reviews,
            "Paper, audit JSON, and reviews must be distinct")
    progress_path = project_path(manifest["latest_progress_path"])
    require(progress_path == (HERE / "progress_checks_1233.json").resolve(), "Wrong latest historical progress")
    require(sha(progress_path) == manifest["latest_progress_sha256"], "Historical progress hash changed")
    mandatory = [paper, audit_path, *reviews, COMMON, SELF, PUBLISHER, UPDATER, RUNNER,
                 progress_path, LEDGER, NOTICE, LAYER, STAGE/"README.md", STAGE/"文件索引.md",
                 *(BASE/name for name in NAV_NAMES)]
    for path in mandatory:
        require_bound(path, bindings)
    audit = read(audit_path)
    require(audit.get("status") == manifest["required_audit_status"], "Stage audit has not passed in declared scope")

    # Freeze the historical scientific/goal part, not its old mutable management.
    progress = read(progress_path)
    require(progress.get("status") == "PASS" and progress.get("latest_round") == 1233,
            "Invalid historical progress receipt")
    require(progress.get("accepted_stage_rounds") == rounds and progress.get("stage_check_groups") == 63,
            "Historical scientific totals changed")
    require(progress.get("stage_complete") is False, "Original progress completion claim changed")
    require(progress.get("application_goal_status") == "active", "Original progress goal status changed")
    require(sha(PUBLISHER) == progress["publisher_sha256"], "Historical publisher changed")
    science = checked_map(progress["science_files_sha256"])
    historical_audits = checked_map(progress["audit_files_sha256"])
    frozen_startup = checked_map(progress["startup_immutable_bindings_sha256"])
    for collection in (science, historical_audits, frozen_startup):
        for name, digest in collection.items():
            require(bindings.get(name) == digest, "Final map omits/differs from frozen evidence: " + name)
    startup_path = project_path(progress["startup_receipt_path"])
    require(sha(startup_path) == progress["startup_receipt_sha256"], "Historical startup receipt changed")
    require_bound(startup_path, bindings)
    startup = read(startup_path)
    require(startup.get("accepted_scientific_rounds") == [] and startup.get("new_scientific_rounds") == 0,
            "Historical startup was rewritten")
    prompt = (HERE/"goal_prompt.txt").read_text(encoding="utf-8").strip()
    objective_sha = hashlib.sha256(prompt.encode("utf-8")).hexdigest()
    start_goal = read(HERE/"goal_start_receipt.json")["goal"]
    require(start_goal["objective"] == prompt and start_goal["status"] == "active", "Original goal changed")
    require(objective_sha == progress["goal_objective_sha256"] == startup["goal_objective_sha256"],
            "Original objective hash differs")

    publisher = load_publisher()
    accepted = []
    require(set(progress["acceptance_sha256"]) == {str(number) for number in rounds},
            "Historical acceptance set differs")
    for number in rounds:
        path = (STAGE/str(number)/"acceptance.json").resolve()
        require_bound(path, bindings)
        require(sha(path) == progress["acceptance_sha256"][str(number)], "Acceptance hash changed")
        acceptance, files = publisher.acceptance(path)
        require(acceptance["round"] == number and acceptance.get("root_executed_both_codes") is True,
                "Missing accepted root recomputation")
        for name, digest in files.items():
            require(science.get(name) == digest and bindings.get(name) == digest,
                    "Acceptance/science/final map differ: " + name)
        # Count the explicitly accepted author groups once. Independent groups
        # are evidence in the bindings but are deliberately never added here.
        accepted.append(dict(round=number, author_check_groups=acceptance["check_groups"],
                             acceptance_sha256=sha(path), evidence_files=len(files)))
    require(sum(item["author_check_groups"] for item in accepted) == 63, "Author group sum differs from 63")

    ledger = read(LEDGER)
    require(ledger.get("status") == manifest["expected_status"], "Current ledger is not completed in declared scope")
    require(ledger.get("first_round") == 1219 and ledger.get("latest_accepted_round") == 1233,
            "Current ledger round endpoints differ")
    require(ledger.get("accepted_stage_rounds") == rounds and ledger.get("stage_check_groups") == 63,
            "Current ledger scientific totals differ")
    require(ledger.get("goal_objective_sha256") == objective_sha, "Current ledger objective changed")
    require(ledger.get("stage_requirement_status", {}).get("stage_deliverables_all_completed") is True,
            "Ledger deliverables are not marked complete")
    require(ledger.get("old_GR_goal_completed") is False, "Old GR scope was incorrectly promoted")
    layer = read(LAYER)
    require(layer.get("notice_path") == rel(NOTICE) and layer.get("notice_sha256") == sha(NOTICE),
            "Layer does not bind current notice")
    fragments = manifest.get("required_notice_fragments", [])
    require(isinstance(fragments, list) and all(isinstance(item,str) for item in fragments),
            "Invalid notice fragments")
    notice = NOTICE.read_text(encoding="utf-8-sig")
    require(all(item in notice for item in fragments), "Final notice omits a required scope fragment")

    navigation_check = run_json(UPDATER, "--check")
    require(navigation_check.get("status") == "CHECKED" and navigation_check.get("changed_files") == 0,
            "Current root navigation does not match current notice")
    navigation_verify = run_json(RUNNER, "--verify-navigation")
    require(navigation_verify.get("status") == "PASS" and navigation_verify.get("manifest_sha256") == sha(LAYER),
            "Current navigation layer verification failed")
    plan_path = project_path(progress["publication_plan_path"])
    require(plan_path == (HERE/"publication_1233.json").resolve(), "Wrong historical publication plan")
    require(sha(plan_path) == progress["publication_plan_sha256"], "Historical plan changed")
    require_bound(plan_path, bindings)
    # Do not call publisher.verify: its old management hashes are historical.
    closure = publisher.closure_check(read(plan_path), science)
    require(closure.get("status") == "EXTERNAL_EDIT" and closure.get("whole_closure_passed") is False,
            "Known historical external edit was concealed or silently rebound")
    require(closure.get("historical_complete_replay_status") == "NOT_PASSED_EXTERNAL_DOCUMENT_EDIT",
            "Historical replay scope was incorrectly promoted")
    require(progress.get("historical_closure_status") == "EXTERNAL_EDIT", "Historical external-edit status changed")

    final_goal = None
    receipt_path = manifest.get("final_goal_receipt_path")
    receipt_sha = manifest.get("final_goal_receipt_sha256")
    require((receipt_path is None) == (receipt_sha is None), "Final goal receipt needs both path and hash")
    if receipt_path is not None:
        receipt = project_path(receipt_path)
        require(sha(receipt) == receipt_sha, "Final goal receipt hash differs")
        require_bound(receipt, bindings)
        goal = read(receipt)["goal"]
        require(goal["objective"] == prompt and isinstance(goal.get("status"),str), "Final goal receipt changed objective")
        final_goal = dict(path=rel(receipt),sha256=receipt_sha,status=goal["status"],
                          objective_sha256=objective_sha)
    # The paper's elementary record-entropy corollary reuses the frozen 1230
    # marginal and mutual-information bounds; it is not another research round.
    marginal_error = Fraction(121, 10**22)
    marginal_entropy_guard = 8*marginal_error**2/(1-4*marginal_error**2)
    record_entropy_gap = Fraction(4, 10**22)-marginal_entropy_guard
    require(marginal_entropy_guard < Fraction(12, 10**40), "Record entropy guard failed")
    require(record_entropy_gap > Fraction(3999, 10**25), "Record entropy gap failed")
    # Detect concurrent edits during the independent subprocesses/imported check.
    require(sha(manifest_path) == manifest_digest, "Final manifest changed during verification")
    require(checked_map(manifest["expected_sha256"]) == bindings, "Evidence changed during verification")
    return dict(schema="thermal_geometry_stage_completion_checks_v1",
                status="PASS_IN_DECLARED_SCOPE", stage_status=manifest["expected_status"],
                scope="Final declared stage deliverables and unchanged accepted science; not unrestricted historical replay",
                manifest_path=rel(manifest_path), manifest_sha256=manifest_digest,
                verifier_path=rel(SELF), verifier_sha256=sha(SELF),
                audit_path=rel(audit_path), audit_status=audit["status"],
                paper_path=rel(paper), stage_review_paths=[rel(path) for path in reviews],
                common_domain_path=rel(COMMON), final_bindings_verified=len(bindings),
                accepted_rounds=rounds, author_check_groups=63,
                independent_groups_counted_again=False, acceptance_checks=accepted,
                record_entropy_corollary=dict(marginal_entropy_guard_exact=str(marginal_entropy_guard),
                    joint_record_entropy_drop_lower_exact=str(record_entropy_gap),
                    drop_strictly_above="3.999e-22 nats", additional_scientific_round=False),
                historical_progress_path=rel(progress_path), historical_progress_sha256=sha(progress_path),
                historical_progress_science_verified=True, historical_management_not_revalidated_as_current=True,
                original_goal_objective_sha256=objective_sha, original_goal_unchanged=True,
                current_ledger_sha256=sha(LEDGER), current_notice_sha256=sha(NOTICE),
                current_navigation_layer_sha256=sha(LAYER),
                navigation_check=navigation_check, navigation_verify=navigation_verify,
                historical_closure_status="EXTERNAL_EDIT", whole_historical_closure_passed=False,
                frozen_closure_verify=closure, final_goal_receipt=final_goal,
                application_goal_changed_by_verifier=False, old_GR_goal_completed=False,
                old_scientific_files_modified=False, navigation_modified=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default="final_evidence_manifest.json")
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    manifest_path = (HERE / args.manifest).resolve()
    require(manifest_path.is_relative_to(ROOT) and manifest_path.is_file(), "Manifest is missing or outside project")
    if args.write_results:
        require(not OUTPUT.exists(), "Completion receipt already exists; it will not be overwritten")
    result = verify(manifest_path)
    result["read_only"] = not args.write_results
    if args.write_results:
        with OUTPUT.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write("\n")
    print(json.dumps(result, ensure_ascii=False, allow_nan=False))


if __name__ == "__main__":
    main()
