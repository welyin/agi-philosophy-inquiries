"""Read-only diagnosis of the legacy round-1009 publication-scope mismatch.

The original failing replay is preserved.  --record only exclusively creates
this diagnosis JSON; default mode compares invariant scientific conditions,
not the growing stage's number of links/files or current navigation hashes.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
PROJECT = STAGE.parent.parent
OLD = STAGE / "1009"
OUT = HERE / "legacy_1009_scope_results.json"
NAVIGATION = {
    "research_cognition_physics/archive_1009_/README.md",
    "research_cognition_physics/archive_1009_/文件索引.md",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evidence() -> dict:
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    sys.dont_write_bytecode = True
    verifier = OLD / "verify_round1009.py"
    receipt = OLD / "research_round_1009_checks.json"
    saved = json.loads(receipt.read_text("utf8"))
    old_sources = saved["source_sha256"]
    assert len(old_sources) == 9
    assert NAVIGATION <= old_sources.keys()
    frozen = {name: digest for name, digest in old_sources.items()
              if name not in NAVIGATION}
    before = {name: sha(PROJECT / name) for name in frozen}
    assert before == frozen
    receipt_before, verifier_before = sha(receipt), sha(verifier)
    sys.path.insert(0, str(OLD))
    spec = importlib.util.spec_from_file_location("legacy_1009_read_only_scope", verifier)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Does not invoke __main__ and passes no write mode.  Every original
    # scientific comparison, historical source check and current link assertion
    # in run() is executed; the publication-time aggregate equality is diagnosed.
    fresh = module.run(False)
    differing = sorted(k for k in saved.keys() | fresh.keys()
                       if saved.get(k) != fresh.get(k))
    allowed = {"source_sha256", "new_stage_local_links_checked"}
    assert set(differing) <= allowed, differing
    new_sources = fresh["source_sha256"]
    missing = sorted(old_sources.keys() - new_sources.keys())
    assert not missing
    changed = sorted(name for name in old_sources
                     if new_sources[name] != old_sources[name])
    assert set(changed) <= NAVIGATION, changed
    assert {name: new_sources[name] for name in frozen} == frozen
    after = {name: sha(PROJECT / name) for name in frozen}
    assert before == after
    assert receipt_before == sha(receipt)
    assert verifier_before == sha(verifier)
    for key in ("scientific_result_reproduced", "archived_scientific_inputs_unchanged",
                "previous_migration_manifest_unchanged", "all_delivery_checks_passed"):
        assert fresh[key] is True
    raw_path = HERE / "stage_replay_results.json"
    raw = json.loads(raw_path.read_text("utf8"))
    row = next(x for x in raw["results"] if x["round"] == 1009)
    assert row["exit_code"] == 1 and row["regression_passed"] is False
    assert "line 69" in row["stderr"] and "AssertionError" in row["stderr"]
    assert raw["failed_rounds"] == [1009]
    assert row["receipt_sha256_before"] == receipt_before
    assert row["verifier_sha256_before"] == verifier_before
    return {
        "round": 1043, "kind": "legacy_publication_scope_diagnosis",
        "new_scientific_calibration_groups": 0,
        "original_1009_default_verifier_passed": False,
        "original_failure_preserved": True,
        "run_false_full_scientific_and_source_checks_passed": True,
        "only_allowed_publication_scope_fields_differ": True,
        "all_seven_frozen_original_assets_unchanged": True,
        "missing_original_sources": missing,
        "permitted_mutable_original_sources": sorted(NAVIGATION),
        "actually_changed_original_sources": changed,
        "differing_top_level_keys_at_this_run": differing,
        "link_counts_at_this_run": {"publication": saved["new_stage_local_links_checked"],
                                    "current": fresh["new_stage_local_links_checked"]},
        "source_counts_at_this_run": {"publication": len(old_sources),
                                      "current": len(new_sources),
                                      "added": len(new_sources.keys() - old_sources.keys())},
        "dynamic_counts_and_navigation_hashes_are_not_scientific_invariants": True,
        "frozen_original_source_sha256": frozen,
        "original_verifier_sha256": verifier_before,
        "original_receipt_sha256": receipt_before,
        "preserved_raw_stage_replay_sha256": sha(raw_path),
        "diagnostic_script_sha256": sha(Path(__file__)),
        "writes_to_original_assets": False,
        "semantic_goal_completion_certified": False,
    }


def stable(value: dict) -> dict:
    excluded = {"actually_changed_original_sources", "differing_top_level_keys_at_this_run",
                "link_counts_at_this_run", "source_counts_at_this_run"}
    return {k: v for k, v in value.items() if k not in excluded}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    args = parser.parse_args()
    if args.record and OUT.exists():
        raise FileExistsError(OUT)
    result = evidence()
    if args.record:
        with OUT.open("x", encoding="utf8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    elif OUT.exists():
        assert stable(result) == stable(json.loads(OUT.read_text("utf8")))
    print(json.dumps({k: v for k, v in result.items()
                      if not k.endswith("sha256")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
