"""Read-only replay, retained history, and local evidence for round 1049."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import top_threshold as science

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RECEIPT = HERE / "research_round_1049_checks.json"
OWN = [
    "../research_note_1049.md", "proof.md", "top_threshold.py",
    "top_threshold_results.json", "NEXT.md", "verify_round1049.py",
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    result = science.run()
    saved = json.loads(science.RESULT.read_text(encoding="utf-8"))
    science.compare(result, saved)
    history = saved["historical_sha256"]
    for rel, expected in history.items():
        assert digest(ROOT / rel) == expected, rel
    link_count = 0
    for rel in OWN:
        path = HERE / rel
        if path.suffix != ".md":
            continue
        for target in re.findall(r"\]\(([^)]+)\)", path.read_text(encoding="utf-8")):
            if target.startswith(("http://", "https://", "#")):
                continue
            dest = target.split("#", 1)[0].strip("<>")
            assert (path.parent / dest).resolve().exists(), (rel, dest)
            link_count += 1
    return {
        "round": 1049,
        "all_scientific_checks_passed": True,
        "new_scientific_groups": 1,
        "new_cognitive_axioms": 0,
        "owned_sha256": {rel: digest(HERE / rel) for rel in OWN},
        "historical_sha256": history,
        "local_links_checked": link_count,
        "full_roadmap_completed": False,
        "scope": saved["scope"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    receipt = verify()
    if args.freeze:
        with RECEIPT.open("x", encoding="utf-8") as out:
            json.dump(receipt, out, ensure_ascii=False, indent=2)
            out.write("\n")
    else:
        assert receipt == json.loads(RECEIPT.read_text(encoding="utf-8"))
    print(json.dumps({
        "round": 1049, "all_passed": True,
        "mode": "exclusive_freeze" if args.freeze else "read_only",
        "owned": len(receipt["owned_sha256"]),
        "history": len(receipt["historical_sha256"]),
        "local_links": receipt["local_links_checked"],
    }))


if __name__ == "__main__":
    main()
