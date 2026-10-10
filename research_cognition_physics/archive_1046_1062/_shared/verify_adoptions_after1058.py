"""Read-only verification of the count-zero interfaces after round 1058.

--freeze creates a receipt exclusively. Navigation may evolve afterwards;
frozen scientific assets and original roadmap requirements may not.
"""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import re
import runpy
import subprocess
import sys
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
PHASE = HERE.parent
ROOT = PHASE.parent
RECEIPT = HERE / "adoptions_after1058_checks.json"
BASE = runpy.run_path(str(HERE / "verify_adoptions_after1055.py"))
ROADMAP_SHA = "0572fd4469c75ffa86d4e08162ca11784dfdbaa983762c520076a5e152f10fa4"
OWN = [
    "gravitational_memory_source_adoption.md",
    "gravitational_memory_source_adoption_review.md",
    "gravitational_memory_source_adoption/check.py",
    "gravitational_memory_source_adoption/results.json",
    "standard_siren_amplitude_adoption.md",
    "standard_siren_amplitude_adoption_review.md",
    "../_admission/after1058_memory_source/selection.md",
    "../_admission/after1058_common_prediction/selection.md",
    "../_admission/after1058_cognitive_selection/selection.md",
    "../_admission/after1058_cognitive_selection/review.md",
    "integration_adoptions_after1058.md",
    "verify_adoptions_after1058.py",
]
HISTORY = [
    "archive_301_341/research_note_322.md",
    "archive_301_341/research_note_326.md",
    "archive_531_553/research_note_549.md",
    "archive_742_763/research_note_763.md",
    "archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md",
    "archive_956_989/981/drafts/common_parent_contract_v1.md",
    "archive_1009_1043/research_note_1033.md",
    "archive_1046_/research_note_1047.md",
    "archive_1046_/research_note_1056.md",
    "archive_1046_/1056/proof.md",
    "archive_1046_/1056/results.json",
    "archive_1046_/1056/mainline_acceptance.json",
    "archive_1046_/1057/mainline_acceptance.json",
    "archive_1046_/1058/mainline_acceptance.json",
    "archive_1046_/_shared/progress_1057_1058_checks.json",
    "archive_1046_/_shared/goal_start.json",
    "archive_1046_/_shared/multimessenger_propagation_adoption.md",
    "archive_1046_/_admission/shared_interfaces_after1056/classical_geometry_selection.md",
]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def local_links(path):
    body = path.read_text(encoding="utf-8")
    if path.parent == ROOT and path.name in BASE["TAILS"]:
        body = body.split(BASE["TAILS"][path.name][0], 1)[0]
    # Display mathematics can contain [g](x); it is not a file link.
    body = re.sub(r"\\\[[\s\S]*?\\\]", "", body)
    body = re.sub(r"\$\$[\s\S]*?\$\$", "", body)
    count = 0
    for target in re.findall(r"\[[^\]\n]*\]\(([^)\n]+)\)", body):
        target = target.strip().strip("<>")
        if re.match(r"[a-zA-Z][\w+.-]*://", target) or target.startswith("#"):
            continue
        assert (path.parent / unquote(target.split("#", 1)[0])).is_file(), (path, target)
        count += 1
    return count


def run(freeze=False):
    algebra = subprocess.run(
        [sys.executable, "-B", "-X", "utf8",
         str(HERE / "gravitational_memory_source_adoption/check.py")],
        capture_output=True, text=True, encoding="utf-8", check=True,
    )
    result = json.loads(algebra.stdout)
    assert result["all_checks_passed"] and result["adoption_count"] == 0
    assert not result["finite_radius_or_time_error_certified"]
    assert not result["actual_gravity_instrument_certified"]
    for rel, expected in result["historical_sha256"].items():
        assert digest(PHASE / rel) == expected, rel

    # Both adoption reviews must sign the actual, unchanged author snapshot.
    for name in ("gravitational_memory_source", "standard_siren_amplitude"):
        author = HERE / f"{name}_adoption.md"
        review = (HERE / f"{name}_adoption_review.md").read_text(encoding="utf-8")
        assert "通过" in review and digest(author) in review.lower(), name
    audit = HERE / "../_admission/after1058_cognitive_selection/selection.md"
    review = audit.with_name("review.md").read_text(encoding="utf-8")
    assert "通过" in review and digest(audit) in review
    independent = json.loads((HERE / "gravitational_memory_source_adoption/review_checks.json").read_text(encoding="utf-8"))
    assert independent["all_checks_passed"] and independent["new_science_groups"] == 0
    for rel, expected in independent["author_sha256"].items():
        assert digest(HERE / rel) == expected, rel
    integration_review = (HERE / "integration_adoptions_after1058_review.md").read_text(encoding="utf-8")
    assert "通过" in integration_review
    assert digest(HERE / "integration_adoptions_after1058.md") in integration_review

    for name, (marker, expected) in BASE["TAILS"].items():
        data = (ROOT / name).read_bytes()
        assert sha256(data[data.index(marker.encode("utf-8")):]).hexdigest() == expected
    data = (ROOT / "ROADMAP.md").read_bytes()
    marker = "## M1 检验认知要求的独立选择力".encode("utf-8")
    assert sha256(data[data.index(marker):]).hexdigest() == ROADMAP_SHA

    owned = list(OWN)
    for rel in ("gravitational_memory_source_adoption/review_checks.json",
                "integration_adoptions_after1058_review.md"):
        if (HERE / rel).exists():
            owned.append(rel)
    assets = {rel: digest(HERE / rel) for rel in owned}
    links = sum(local_links(HERE / rel) for rel in owned if rel.endswith(".md"))
    navlinks = sum(local_links(p) for p in BASE["NAV"])
    stable = {
        "date": "2026-10-09", "kind": "mature_interfaces_and_selection_audit",
        "latest_accepted_round_at_batch": 1058, "cumulative_scientific_groups": 3834,
        "stage_scientific_groups": 13, "new_scientific_groups": 0,
        "new_empirical_groups": 0, "new_cognitive_axioms": 0,
        "assets_sha256": assets,
        "historical_sha256": {rel: digest(ROOT / rel) for rel in HISTORY},
        "asset_local_links_checked": links,
        "original_roadmap_requirements_sha256": ROADMAP_SHA,
        "independent_adoption_reviews_passed": True,
        "lightweight_algebra_readonly_check_passed": True,
        "newly_excluded_accepted_physical_branches": [],
        "finite_memory_detector_certified": False,
        "whole_roadmap_complete": False,
        "goal_tool_status_observed": "active",
        "goal_objective_or_status_mutated": False,
        "all_checks_passed": True,
    }
    if freeze:
        record = dict(stable)
        record["navigation_snapshot"] = {
            "sha256": {str(p.relative_to(ROOT)): digest(p) for p in BASE["NAV"]},
            "local_links_checked": navlinks,
        }
        with RECEIPT.open("x", encoding="utf-8") as stream:
            json.dump(record, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        record = json.loads(RECEIPT.read_text(encoding="utf-8"))
        record.pop("navigation_snapshot")
        assert record == stable, "Frozen batch changed; preserve history and create a new batch."
    return stable, navlinks


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    record, navlinks = run(args.freeze)
    print(json.dumps({"all_checks_passed": True,
                      "assets": len(record["assets_sha256"]),
                      "history": len(record["historical_sha256"]),
                      "current_navigation_links": navlinks,
                      "new_science_groups": 0,
                      "latest_round_at_batch": 1058,
                      "whole_roadmap_complete": False}, ensure_ascii=False))
