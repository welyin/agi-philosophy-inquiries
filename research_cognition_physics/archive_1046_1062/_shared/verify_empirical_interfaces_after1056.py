"""Freeze/read-only checks of the empirical-interface adoption batch.

This signs bounded adoptions and history integrity, not new physical proofs.
The default does not write. --freeze exclusively creates the initial receipt.
"""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import re
import runpy
import subprocess
import sys

HERE = Path(__file__).resolve().parent
BASE = runpy.run_path(str(HERE / "verify_adoptions_after1055.py"))
ROOT, PHASE = BASE["ROOT"], BASE["PHASE"]
RECEIPT = HERE / "empirical_interfaces_after1056_checks.json"
OWN = [
    "atom_interferometry_source_adoption.md",
    "atom_interferometry_source_adoption_review.md",
    "multimessenger_propagation_adoption.md",
    "multimessenger_propagation_adoption_review.md",
    "multimessenger_propagation_adoption/check.py",
    "multimessenger_propagation_adoption/results.json",
    "../_admission/shared_interfaces_after1056/field_instrument_selection.md",
    "../_admission/shared_interfaces_after1056/spin_statistics_selection.md",
    "integration_empirical_interfaces_after1056.md",
    "verify_empirical_interfaces_after1056.py",
    "_history/verify_empirical_interfaces_after1056_initial.py",
    "_history/empirical_interfaces_after1056_initial_checks.json",
]
HISTORY = [
    "archive_956_989/981/drafts/common_parent_contract_v1.md",
    "archive_935_955/research_note_941.md",
    "archive_956_989/research_note_974.md", "archive_956_989/research_note_979.md",
    "archive_1009_1043/research_note_1019.md", "archive_1009_1043/research_note_1024.md",
    "archive_1009_1043/research_note_1035.md",
    "archive_1009_1043/research_note_1040.md", "archive_1009_1043/research_note_1041.md",
    "archive_1046_/1056/mainline_acceptance.json",
    "archive_1046_/_shared/integration_1056.md",
    "archive_1046_/_shared/navigation_1056_checks.json",
    "archive_1046_/_shared/operational_adoptions_after1055_checks.json",
    "archive_1046_/_shared/atomic_clock_redshift_adoption.md",
    "archive_1046_/_shared/equivalence_principle_adoption.md",
    "archive_1046_/_shared/optomechanical_readout_adoption.md",
]
ROADMAP_TAIL_SHA = "0572fd4469c75ffa86d4e08162ca11784dfdbaa983762c520076a5e152f10fa4"


def checked_links(path, pending=False):
    content = path.read_text(encoding="utf-8")
    if path.parent == ROOT and path.name in BASE["TAILS"]:
        content = content.split(BASE["TAILS"][path.name][0], 1)[0]
    # Exclude display math: e.g. [T+C[g](x),A] is not a Markdown link.
    content = re.sub(r"\\\[.*?\\\]", "", content, flags=re.S)
    pending_links = 0
    if pending:
        content, pending_links = re.subn(r"\]\((?:archive_1046_/_shared/|_shared/)?empirical_interfaces_after1056_checks\.json\)",
                                        "](#pending-receipt)", content)
    # A to-be-created receipt is still one local link. Keep the asset metadata
    # identical between the initial freeze and subsequent read-only runs.
    return BASE["local_links"](path, False, content) + pending_links


def check_script(path):
    result = subprocess.run([sys.executable, "-B", "-X", "utf8", str(path)],
                            check=True, capture_output=True, text=True, encoding="utf-8")
    return json.loads(result.stdout)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    if args.freeze:
        assert not RECEIPT.exists(), "Preserve frozen receipts."
    numerical = check_script(HERE / "multimessenger_propagation_adoption/check.py")
    assert numerical["all_checks_passed"] and numerical["new_science_groups"] == 0
    assert not numerical["raw_data_fit"] and not numerical["new_confidence_interval"]
    for topic in ("atom_interferometry_source", "multimessenger_propagation"):
        report = HERE / f"{topic}_adoption.md"
        review = (HERE / f"{topic}_adoption_review.md").read_text(encoding="utf-8")
        assert "通过" in review and BASE["digest"](report) in review.lower()
    assets = []
    for rel in OWN:
        path = (HERE / rel).resolve()
        links = checked_links(path, args.freeze) if path.suffix == ".md" else 0
        assets.append({"file": rel, "sha256": BASE["digest"](path), "local_links": links})
    history = {rel: BASE["digest"](ROOT / rel) for rel in HISTORY}
    for name, (marker, expected) in BASE["TAILS"].items():
        raw = (ROOT / name).read_bytes()
        assert sha256(raw[raw.index(marker.encode("utf-8")):]).hexdigest() == expected
    raw = (ROOT / "ROADMAP.md").read_bytes()
    marker = "## M1 检验认知要求的独立选择力".encode("utf-8")
    assert sha256(raw[raw.index(marker):]).hexdigest() == ROADMAP_TAIL_SHA
    navigation_links = sum(checked_links(p, args.freeze) for p in BASE["NAV"])
    stable = {
        "date": "2026-10-09", "kind": "bounded_empirical_adoption_and_history_audit",
        "assets": assets, "history_sha256": history,
        "original_roadmap_requirements_sha256": ROADMAP_TAIL_SHA,
        "latest_accepted_round": 1056, "cumulative_scientific_groups": 3832,
        "stage_scientific_groups": 11, "new_scientific_groups": 0,
        "new_cognitive_axioms": 0, "independent_reviews_passed": True,
        "propagation_dictionary_readonly_check_passed": True,
        "new_experimental_fit": False, "same_global_process_certified": False,
        "whole_roadmap_complete": False, "goal_status_at_adoption": "active",
        "initial_receipt_link_count_bug_preserved_in_history": True,
        "all_passed": True,
    }
    if args.freeze:
        output = dict(stable)
        output["navigation_snapshot"] = {
            "sha256": {str(p.relative_to(ROOT)): BASE["digest"](p) for p in BASE["NAV"]},
            "local_link_occurrences": navigation_links,
        }
        with RECEIPT.open("x", encoding="utf-8") as out:
            json.dump(output, out, ensure_ascii=False, indent=2)
            out.write("\n")
    else:
        saved = json.loads(RECEIPT.read_text(encoding="utf-8"))
        saved.pop("navigation_snapshot")
        assert stable == saved, "Frozen asset/history mismatch."
    print(json.dumps({"all_passed": True, "assets": len(assets), "history": len(history),
                      "current_navigation_local_links": navigation_links,
                      "latest_accepted_round": 1056, "new_science_groups": 0}))


if __name__ == "__main__":
    main()
