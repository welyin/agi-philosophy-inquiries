"""Read-only science/history/link replay; --freeze exclusively creates receipt."""
from pathlib import Path
import argparse
import json
import re
import weak_decay_source as science

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RECEIPT = HERE / "research_round_1056_checks.json"
OWN = ["../research_note_1056.md", "proof.md", "sources.md", "dependency_update.md",
       "weak_decay_source.py", "results.json", "NEXT.md", "verify_round1056.py"]


def run(freeze):
    current = science.run()
    saved = json.loads(science.RESULT.read_text(encoding="utf-8"))
    science.compare(current, saved)
    for rel, digest in saved["historical_sha256"].items():
        assert science.digest(ROOT/rel) == digest, rel
    link_count = 0
    for name in OWN:
        path = HERE/name
        assert path.is_file(), name
        if path.suffix != ".md":
            continue
        body = path.read_text(encoding="utf-8")
        body = re.sub(r"\\\[[\s\S]*?\\\]", "", body)
        body = re.sub(r"\$\$[\s\S]*?\$\$", "", body)
        for target in re.findall(r"\[[^\]\n]*\]\(([^)\n]+)\)", body):
            target = target.strip().strip("<>")
            if target.startswith(("https://", "http://", "#")):
                continue
            destination = (path.parent/target.split("#", 1)[0]).resolve()
            assert destination.is_file(), (name, target)
            link_count += 1
    record = {
        "round": 1056, "date": "2026-10-09", "author_checks_passed": True,
        "owned_sha256": {name: science.digest(HERE/name) for name in OWN},
        "historical_sha256": saved["historical_sha256"],
        "local_links_checked": link_count,
        "rational_bound_over_m": "-1/4500",
        "record_only_recovery_error_strictly_greater_than": "1/4500",
        "scientific_addition": "one conditional joint weak-decay record/source dictionary",
        "new_project_groups": 1, "new_cognitive_axioms": 0,
        "whole_roadmap_complete": False,
        "independent_final_review_complete_at_author_freeze": False,
        "scope": saved["scope"],
    }
    if freeze:
        with RECEIPT.open("x", encoding="utf-8", newline="\n") as output:
            json.dump(record, output, ensure_ascii=False, indent=2)
            output.write("\n")
    else:
        assert record == json.loads(RECEIPT.read_text(encoding="utf-8"))
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    result = run(args.freeze)
    print(json.dumps({"round": 1056, "all_passed": True, "assets": len(OWN),
                      "history": len(result["historical_sha256"]),
                      "links": result["local_links_checked"],
                      "mode": "exclusive_freeze" if args.freeze else "read_only"}))
