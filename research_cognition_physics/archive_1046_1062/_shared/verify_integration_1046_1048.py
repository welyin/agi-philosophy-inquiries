"""Verify the accepted evidence manifest without repeating scientific runs."""
from pathlib import Path
import argparse
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
PHASE = HERE.parent
ROOT = PHASE.parent
OUT = HERE / "acceptance_1046_1048.json"
REVIEWERS = {
    1046: "/root/round1037_dynamics_selection",
    1047: "/root/next_structural_selection",
    1048: "/root/phase1044_goal_scope_review",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    rounds = []
    all_history = {}
    for n, reviewer in REVIEWERS.items():
        base = PHASE / str(n)
        receipt = base / f"research_round_{n}_checks.json"
        checked = json.loads(receipt.read_text(encoding="utf-8"))
        assert checked.get("all_checks_passed",
                           checked.get("all_scientific_checks_passed", False))
        owned = checked["owned_sha256"]
        for rel, expected in owned.items():
            assert sha(base / rel) == expected, (n, rel)
        history = checked["historical_sha256"]
        for rel, expected in history.items():
            assert sha(ROOT / rel) == expected, (n, rel)
            assert all_history.get(rel, expected) == expected, rel
            all_history[rel] = expected
        review = base / "independent_review.md"
        assert review.is_file(), review
        # Approval is a mainline judgement after reading the review, not inferred
        # by scanning for a word. Freeze the actual reviewed evidence here.
        rounds.append({
            "round": n, "accepted_in_declared_scope": True,
            "independent_agent_reviewer": reviewer,
            "review_sha256": sha(review),
            "author_receipt_sha256": sha(receipt),
            "owned_sha256": owned,
            "new_scientific_groups": 1,
            "new_cognitive_axioms": 0,
        })
    integration = HERE / "integration_1046_1048.md"
    link_count = 0
    # The stage index is live navigation and must remain free to gain rounds.
    # Only this batch's immutable integration belongs in the frozen receipt.
    for p in (integration,):
        for target in re.findall(r"\]\(([^)]+)\)", p.read_text(encoding="utf-8")):
            if target.startswith(("https://", "http://", "#")):
                continue
            assert (p.parent / target.split("#", 1)[0]).resolve().exists(), target
            link_count += 1
    return {
        "schema": "roadmap_increment_acceptance_v2",
        "date": "2026-10-08",
        "mainline_acceptance": "passed_with_explicit_scopes",
        "rounds": rounds,
        "latest_completed_round": 1048,
        "baseline_completed_round": 1045,
        "baseline_scientific_groups": 3821,
        "new_scientific_groups": 3,
        "cumulative_scientific_groups": 3824,
        "new_cognitive_axioms": 0,
        "historical_sha256": all_history,
        "integration_sha256": sha(integration),
        "verifier_sha256": sha(Path(__file__)),
        "frozen_integration_links_checked": link_count,
        "live_navigation_excluded_from_freeze": True,
        "initial_receipt_preserved": "_history/acceptance_1046_1048_initial.json",
        "goal_start_sha256": sha(HERE / "goal_start.json"),
        "goal_status": "active",
        "roadmap_complete": False,
        "single_parent_unification_certified": False,
        "new_empirical_validation": False,
        "checks": {
            "author_science_replays": "Passed and independently reviewed; see each round evidence",
            "this_run": "Hashes, evidence receipts and integration links only",
            "additional_scientific_groups_from_audit": 0,
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = build()
    if args.write:
        with OUT.open("x", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write("\n")
    else:
        assert result == json.loads(OUT.read_text(encoding="utf-8"))
    print(json.dumps({
        "rounds": [r["round"] for r in result["rounds"]],
        "accepted": True,
        "historical_assets": len(result["historical_sha256"]),
        "cumulative_scientific_groups": result["cumulative_scientific_groups"],
        "goal_status": "active", "roadmap_complete": False,
        "mode": "exclusive_write" if args.write else "read_only",
    }))


if __name__ == "__main__":
    main()
