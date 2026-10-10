"""Read-only scientific replay, history/link check, and author-asset receipt."""
from pathlib import Path
import argparse
import json
import re
import coulomb_newton_source as science

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RECEIPT = HERE / "research_round_1053_checks.json"
OWN = ["../research_note_1053.md", "proof.md", "coulomb_newton_source.py",
       "results.json", "dependency_update.md", "NEXT.md", "review.md",
       "sources.md", "verify_round1053.py"]


def verify(freeze=False):
    current = science.run()
    saved = json.loads(science.RESULT.read_text(encoding="utf-8"))
    science.compare(current, saved)
    for rel, sha in saved["historical_sha256"].items():
        assert science.digest(ROOT/rel) == sha, rel
    links = 0
    for rel in OWN:
        path = HERE/rel
        assert path.is_file(), rel
        if path.suffix != ".md":
            continue
        for target in re.findall(r"\]\(([^)]+)\)", path.read_text(encoding="utf-8")):
            target = target.strip().strip("<>")
            if target.startswith(("https://", "http://", "#")):
                continue
            dest = (path.parent/target.split("#", 1)[0]).resolve()
            if not (freeze and dest == RECEIPT):
                assert dest.exists(), (rel, target)
            links += 1
    return {
        "round": 1053, "date": "2026-10-08", "author_checks_passed": True,
        "independent_review_completed_at_author_freeze": False,
        "new_conditional_physical_source_calibration_groups": 1,
        "new_cognitive_axioms": 0,
        "M4_or_full_roadmap_complete": False,
        "owned_sha256": {rel: science.digest(HERE/rel) for rel in OWN},
        "historical_sha256": saved["historical_sha256"],
        "local_links_checked": links,
        "source_dictionary_bound_exact": saved["rational_certificate"]["exact"]["driven_source_error_strict_upper"],
        "source_scope": saved["scope"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    output = verify(args.freeze)
    if args.freeze:
        with RECEIPT.open("x", encoding="utf-8") as stream:
            json.dump(output, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        assert output == json.loads(RECEIPT.read_text(encoding="utf-8"))
    print(json.dumps({"round": 1053, "author_checks_passed": True,
                      "mode": "exclusive_freeze" if args.freeze else "read_only",
                      "owned": len(OWN), "history": len(output["historical_sha256"]),
                      "local_links": output["local_links_checked"],
                      "independent_review_at_freeze": False}))


if __name__ == "__main__":
    main()
