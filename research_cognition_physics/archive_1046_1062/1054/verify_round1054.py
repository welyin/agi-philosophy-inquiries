"""Read-only science/history/link replay; --freeze exclusively creates the receipt."""
from pathlib import Path
import argparse
import json
import re
import optical_probability_coordinates as science

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RECEIPT = HERE / "research_round_1054_checks.json"
OWN = ["../research_note_1054.md", "proof.md", "optical_probability_coordinates.py",
       "results.json", "dependency_update.md", "NEXT.md", "review.md", "sources.md",
       "verify_round1054.py"]


def verify(freeze=False):
    current = science.run()
    saved = json.loads(science.RESULT.read_text(encoding="utf-8"))
    science.compare(current, saved)
    for rel, value in saved["historical_sha256"].items():
        assert science.digest(ROOT/rel) == value, rel
    link_count = 0
    for rel in OWN:
        path = HERE/rel
        assert path.is_file(), rel
        if path.suffix != ".md":
            continue
        for target in re.findall(r"\]\(([^)]+)\)", path.read_text(encoding="utf-8")):
            target = target.strip().strip("<>")
            if target.startswith(("http://", "https://", "#")):
                continue
            destination = (path.parent/target.split("#", 1)[0]).resolve()
            if not (freeze and destination == RECEIPT):
                assert destination.exists(), (rel, target)
            link_count += 1
    return {
        "round": 1054, "date": "2026-10-08", "author_checks_passed": True,
        "independent_review_completed_at_author_freeze": False,
        "new_project_conditional_physical_dictionary_groups": 1,
        "new_cognitive_axioms": 0, "M3A_or_roadmap_completed_here": False,
        "owned_sha256": {rel: science.digest(HERE/rel) for rel in OWN},
        "historical_sha256": saved["historical_sha256"],
        "local_links_checked": link_count,
        "exact_slope_lower": saved["rational_certificate"]["exact"]["exact_slope_lower"],
        "inverse_lipschitz_upper": "200000000",
        "scope": saved["scope"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    result = verify(args.freeze)
    if args.freeze:
        with RECEIPT.open("x", encoding="utf-8") as output:
            json.dump(result, output, ensure_ascii=False, indent=2)
            output.write("\n")
    else:
        assert result == json.loads(RECEIPT.read_text(encoding="utf-8"))
    print(json.dumps({"round": 1054, "author_checks_passed": True,
                      "mode": "exclusive_freeze" if args.freeze else "read_only",
                      "owned": len(OWN), "history": len(result["historical_sha256"]),
                      "local_links": result["local_links_checked"]}))


if __name__ == "__main__":
    main()
