"""Delivery checks for 962. These are not a proof of physical unification."""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import re
import sys

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/"scripts"))
from research_layout import Layout
TARGET = HERE/"research_round_962_checks.json"


def read(p):
    return json.loads(p.read_text("utf-8-sig"))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(writing=False):
    frozen = {}
    def add(d):
        for k, v in d.items():
            assert k not in frozen or frozen[k] == v, k
            frozen[k] = v
    for n in range(776, 962):
        old = read(STAGE/f"{n}/research_round_{n}_checks.json")
        for k in ("frozen_inputs", "new_scientific_and_entry_files"):
            add(old[k])
    for path in ("875/drafts/clock_transport_working_checks.json",
                 "878/drafts/working_checks.json", "884/drafts/working_checks.json",
                 "887/drafts/working_checks.json", "894/drafts/working_checks.json",
                 "896/drafts/working_checks.json", "897/drafts/working_checks.json",
                 "899/drafts/working_checks.json", "900/drafts/working_checks.json",
                 "904/drafts/working_checks.json", "907/drafts/working_checks.json",
                 "909/drafts/working_checks.json", "912/drafts/working_checks.json",
                 "913/drafts/working_checks.json", "915/drafts/working_checks.json",
                 "917/drafts/working_checks.json"):
        add(read(STAGE/path)["files"])
    extra = {
        "909/drafts/scope_reaudit_checks.json": ("preserved_files", "evidence_and_audit_hashes"),
        "911/drafts/finite_scope_checkpoint_checks.json": ("evidence_and_audit_hashes",),
        "912/drafts/effective_scope_reaudit_checks.json": ("evidence_and_audit_hashes",),
        "914/drafts/finite_scope_decision_checks.json": ("evidence_hashes", "new_document_and_verifier_hashes"),
        "919/drafts/effective_scope_after_918_checks.json": ("frozen_inputs_verified", "new_document_and_verifier_hashes"),
        "953/drafts/priority_reaudit_checks.json": ("frozen_evidence_hashes", "audit_file_hashes")}
    for path, keys in extra.items():
        value = read(STAGE/path)
        for key in keys:
            add(value[key])
    for rel, digest in frozen.items():
        assert sha(ROOT/rel) == digest, rel
    layout = Layout().verify()
    result = read(HERE/"clock_relation_screen_results.json")
    assert result["round"] == 962 and result["all_scientific_checks_passed"]
    for rel, digest in result["source_hashes"].items():
        assert sha(ROOT/rel) == digest, rel
    assert result["raw_and_encoded_unitary_max_error"] < 1e-12
    assert result["analytic_nonzero_clock_source_certificate"]["lower_float"] > 1e-10
    last = result["cases"][-1]
    assert last["graph_trace_distance"] > .08 and .99 < last["clock_visibility"] < .999
    assert last["total_current_balance_error"] < 1e-12
    assert result["scope"]["stop_this_candidate_optimization"]
    assert not result["scope"]["full_goal_completed"]
    ack = read(HERE/"drafts/app_goal_confirmation.json")
    goal = (HERE/"drafts/revised_goal_20261007.txt").read_text("utf-8")
    assert ack["status"] == "active" and ack["objective"].strip() == goal.strip()
    assert ack["edited_by_user"] and not ack["modified_by_assistant"]
    note = STAGE/"research_note_962.md"
    prose = note.read_text("utf-8")
    assert prose.count("$$") == 16
    assert re.findall(r"\\tag\{(\d+)\}", prose) == [str(n) for n in range(1, 9)]
    newfiles = [note, HERE/"clock_relation_screen.py", HERE/"clock_relation_screen_results.json",
                Path(__file__), STAGE/"963/drafts/STATUS.md"] + [
        HERE/"drafts"/n for n in ("revised_goal_20261007.txt", "previous_goal_section.md",
                                 "update_goal_documents.py", "app_goal_confirmation.json",
                                 "update_navigation962.py")]
    for p in newfiles:
        if p.suffix == ".py":
            ast.parse(p.read_text("utf-8"))
    nav = [STAGE.parent/n for n in ("README.md", "research_direction.md", "RESEARCH_STATE.md")] + [
        STAGE/n for n in ("README.md", "文件索引.md", "阶段成果总览.md", "跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping = nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至962）", 1)[1].split("### 当前取舍", 1)[0]
    assert re.findall(r"^\|(C\d\d) ", mapping, re.M) == [f"C{i:02d}" for i in range(1, 28)]
    links = 0
    for doc in [p for p in newfiles if p.suffix == ".md"] + nav:
        content = re.sub(r"\$\$.*?\$\$", "", doc.read_text("utf-8-sig"), flags=re.S)
        for link in re.findall(r"\]\(([^)]+)\)", content):
            if re.match(r"^[a-zA-Z]+://", link) or link.startswith("#"):
                continue
            target = (doc.parent/link.split("#")[0].strip("<>")).resolve()
            assert target.exists() or (writing and target == TARGET.resolve()), (doc, link)
            links += 1
    nums = []
    for p in STAGE.parent.rglob("research_note_*.md"):
        m = re.fullmatch(r"research_note_(\d+).md", p.name)
        if m and p.parent.name.startswith("archive_"):
            nums.append(int(m[1]))
    assert sorted(n for n in nums if n <= 962) == list(range(1, 963))
    assert "001—962轮共962份" in nav[0].read_text("utf-8-sig")
    assert "231—962的732份" in nav[5].read_text("utf-8-sig")
    for p in nav:
        assert "962：钟保护、关系反馈与稳定性的量词" in p.read_text("utf-8-sig")
    oldfiles = [STAGE/"961/research_round_961_checks.json", HERE/"drafts/STATUS.md"]
    prev = read(oldfiles[0])
    out = dict(round=962, date="2026-10-07", all_delivery_checks_passed=True,
        fresh_test_groups=1, formal_reports=962,
        cumulative_numbered_test_groups_from_961=prev["cumulative_numbered_test_groups_from_960"]+1,
        historical_unique_files_verified=len(frozen), historical_manifest_evidence=layout,
        local_links_checked=links, app_goal_user_updated_and_active=True,
        old_results_reused_and_not_claimed_as_new=True,
        M2_precision_scope_clarified=True, geometry_and_unification_not_derived=True,
        full_goal_completed=False, visual_checks_performed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
    if not writing:
        before = read(TARGET)
        for k in ("frozen_inputs", "new_scientific_and_entry_files"):
            assert before[k] == out[k], k
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.write:
        assert not TARGET.exists()
    out = run(args.write)
    if args.write:
        with TARGET.open("x", encoding="utf-8") as dest:
            json.dump(out, dest, ensure_ascii=False, indent=2)
            dest.write("\n")
    print(json.dumps({k:v for k,v in out.items() if k not in
                    ("frozen_inputs", "new_scientific_and_entry_files")},
                     ensure_ascii=False, indent=2))
