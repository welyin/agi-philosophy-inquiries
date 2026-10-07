"""Verify saved evidence and the priority correction; no new science experiment."""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1]
RESEARCH = STAGE.parent
ROOT = RESEARCH.parent
TARGET = HERE / "priority_reaudit_checks.json"

def read(p):
    return json.loads(p.read_text("utf-8-sig"))

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def run():
    files = {}
    def check_manifest(manifest):
        for rel, digest in manifest.items():
            assert sha(ROOT / rel) == digest, rel
            assert rel not in files or files[rel] == digest, rel
            files[rel] = digest
    # Freeze preservation is checked against the original delivery receipts.
    for n in range(947, 953):
        receipt = read(STAGE / f"{n}/research_round_{n}_checks.json")
        assert receipt["all_delivery_checks_passed"]
        check_manifest(receipt["frozen_inputs"])
        check_manifest(receipt["new_scientific_and_entry_files"])
    sources = {
        947: "protocol_field_transport",
        948: "neutral_matter_bridge",
        950: "material_vacuum_matching",
        951: "body_sector_bridge",
        952: "portal_lapse_source",
    }
    results = {}
    for n, stem in sources.items():
        result = read(STAGE / f"{n}/{stem}_results.json")
        assert result["all_scientific_checks_passed"]
        check_manifest(result["source_hashes"])
        results[n] = result
    q = results[947]["finite_joint_control"]["full_cq_instrument_distance_bound"]
    assert q == results[951]["inherited_by_exact_sector_identity"]["full_cq_instrument_distance_bound"]
    local = results[950]["local_matching_counterexample"]["common_relative_response_change"]
    assert abs(local - 2 / 43) < 1e-14
    finite = results[952]["finite_lapse"]
    assert .0005165 < finite["wrong_minus_correct_record_W"] < .0005167
    assert not results[952]["scope"]["new_947_finite_time_probability_or_preparation_bound"]
    assert not results[951]["scope"]["full_SM_Einstein_parent_matching_certified"]
    nav = [RESEARCH / n for n in ("README.md", "research_direction.md", "RESEARCH_STATE.md")]
    nav += [STAGE / n for n in ("README.md", "文件索引.md", "阶段成果总览.md",
        "跨阶段主题索引.md", "_shared/notes/unified_physics_condition_ledger_current.md")]
    docs = [HERE / "priority_reaudit.md"] + nav
    links = 0
    for p in docs:
        text = p.read_text("utf-8-sig")
        if p in nav:
            assert "953工作期执行更正：先判共同接合，暂缓动态细化" in text
            assert "正式952／累计3737" in text and "priority_reaudit.md" in text
            assert not re.search(r"接\[953\]\([^)]+STATUS\.md\)核", text)
        for link in re.findall(r"\]\(([^)]+)\)", text):
            if re.match(r"^[a-zA-Z]+://", link) or link.startswith("#"):
                continue
            target = (p.parent / link.split("#")[0].strip("<>")).resolve()
            assert target.exists() or target == TARGET.resolve(), (str(p), link)
            links += 1
    assert not (STAGE / "research_note_953.md").exists(), "This audit is not a completed scientific round"
    ledger = nav[-1].read_text("utf-8-sig").split("## 六条共同协议：全局缺口对应与检验优先级（截至952）", 1)[1].split("### 当前取舍", 1)[0]
    assert re.findall(r"^\|(C\d\d) ", ledger, re.M) == [f"C{i:02d}" for i in range(1, 28)]
    created = [HERE / n for n in ("priority_reaudit.md", "verify_priority_reaudit.py", "update_priority_navigation.py")]
    for p in created:
        if p.suffix == ".py":
            ast.parse(p.read_text("utf-8"))
    return {
        "date": "2026-10-07", "all_audit_checks_passed": True,
        "formal_round": 952, "cumulative_test_groups": 3737, "new_scientific_test_groups": 0,
        "round953_scientifically_completed": False, "full_goal_completed": False,
        "original_dynamic_plan_deferred": True, "goal_or_automation_changed": False,
        "old_evidence_values": {"instrument_bound": q, "local_matching_relative_change": local,
            "lapse_kernel_difference": finite["wrong_minus_correct_record_W"]},
        "frozen_evidence_hashes": files, "local_links_checked": links,
        "audit_file_hashes": {str(p.relative_to(ROOT)): sha(p) for p in created},
        "navigation_snapshot_hashes": {str(p.relative_to(ROOT)): sha(p) for p in nav},
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.write:
        assert not TARGET.exists(), "Never overwrite the audit receipt"
    out = run()
    if args.write:
        with TARGET.open("x", encoding="utf-8") as f:
            json.dump(out, f, ensure_ascii=False, indent=2)
            f.write("\n")
    else:
        before = read(TARGET)
        for key in ("frozen_evidence_hashes", "audit_file_hashes", "old_evidence_values"):
            assert out[key] == before[key], key
    print(json.dumps({k: v for k, v in out.items() if not k.endswith("hashes")}, ensure_ascii=False, indent=2))
