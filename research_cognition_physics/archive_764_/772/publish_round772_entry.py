"""Publish a verified working entry; do not increment scientific round counters."""
import ast
import importlib.util
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
TARGET = HERE/"round772_drafts/first_response_entry_checks.json"
assert not TARGET.exists()
history = dict(core.read(HERE/"round584_drafts/historical_evidence_manifest.json")["evidence_hashes"])
for n in range(584, 772):
    receipt = core.read(HERE/f"research_round_{n}_checks.json")
    for key in ("new_file_hashes", "preserved_draft_hashes"):
        for name, digest in receipt[key].items():
            assert name not in history or history[name] == digest
            history[name] = digest
history.update(core.read(HERE/"cognitive_foundation_bridge_605_navigation.json")["supplementary_artifact_hashes"])
assert len(history) == 3773
for name, digest in history.items():
    assert core.digest(HERE/name) == digest, name
artifacts = (
    "round772_drafts/research_note_772_working.md",
    "round772_drafts/first_response_entry.py",
    "round772_drafts/first_response_entry_results.json",
    "round772_drafts/entry_scope_review.json",
    "round772_drafts/entry_condition_ledger.md",
)
spec = importlib.util.spec_from_file_location("entry772calibration", HERE/artifacts[1])
calibration = importlib.util.module_from_spec(spec)
spec.loader.exec_module(calibration)
assert calibration.run() == core.read(HERE/artifacts[2])
ast.parse((HERE/artifacts[1]).read_text("utf8"))
report_text = (HERE/artifacts[0]).read_text("utf8")
assert report_text.count("$$") == 14
for link in core.link_parser()(report_text):
    assert ((HERE/artifacts[0]).parent/link).resolve().exists(), link
for extra in artifacts[4:]:
    for link in core.link_parser()((HERE/extra).read_text("utf8")):
        assert ((HERE/extra).parent/link).resolve().exists(), link
text_checks = core.text_checks(HERE/artifacts[0])
assert text_checks["display_formulas"] == 7
paths = [ROOT/"README.md", RESEARCH/"README.md", RESEARCH/"research_direction.md",
         RESEARCH/"RESEARCH_STATE.md", HERE/"README.md",
         HERE/"spatial_premise_closure_audit.md", HERE/"three_dimensional_four_conditions_review.md"]
before = {p: p.read_bytes() for p in paths}
snapshot = HERE/"navigation_before_round772_entry_20261004"
assert not snapshot.exists()
summary = '**772首阶响应入口已核，正式仍771／3499：** [工作报告]({p}round772_drafts/research_note_772_working.md)把771完整来源接入旧732／754约束和Cauchy响应；这是直接应用，不新增轮次。保原自由态，下一项核有限阶共同插入与物理态/记录提升，不以全阶背景输送拖延首阶均值。[条件增量]({p}round772_drafts/entry_condition_ledger.md)、[入口核验]({p}round772_drafts/first_response_entry_checks.json)。目标及旧限定结论保持。'
planned = {}
for p, raw in before.items():
    enc = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf8"
    nl = "\r\n" if b"\r\n" in raw else "\n"
    txt = raw.decode(enc).replace("\r\n", "\n")
    assert "**772首阶响应入口已核" not in txt
    prefix = "research_cognition_physics/archive_231_/" if p.parent == ROOT else "archive_231_/" if p.parent == RESEARCH else ""
    addition = summary.format(p=prefix)
    if p in paths[:5]:
        head, tail = txt.split("\n\n", 1)
        txt = head+"\n\n"+addition+"\n\n"+tail
    else:
        txt += "\n\n"+addition+"\n"
    planned[p] = txt.replace("\n", nl).encode(enc)
links = 0
for p, raw in planned.items():
    for link in core.link_parser()(raw.decode("utf-8-sig")):
        target = (p.parent/link).resolve()
        assert target.exists() or target == TARGET.resolve(), (p, link)
        links += 1
assert all(p.read_bytes() == raw for p, raw in before.items())
snapshot.mkdir()
manifest = {}
for i, (p, raw) in enumerate(before.items()):
    name = f"{i}_{p.name}"
    with (snapshot/name).open("xb") as f:
        f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()] = dict(snapshot=name, sha256=hashlib.sha256(raw).hexdigest())
with (snapshot/"manifest.json").open("x", encoding="utf8") as f:
    json.dump(manifest, f, ensure_ascii=False, indent=2)
for p, raw in planned.items():
    assert p.read_bytes() == before[p]
    temp = p.with_name(p.name+".round772entry.tmp")
    with temp.open("xb") as f:
        f.write(raw)
    os.replace(temp, p)
report = dict(status="working entry only", latest_completed_round=771,
              cumulative_scientific_tests=3499, protected_scientific_evidence_files=3773,
              previous_evidence_unchanged=True, new_scientific_tests=0, bookkeeping_audit_reproduced=True, entry_calibrations_reproduced=0, entry_text_checks=text_checks,
              working_artifact_hashes={n: core.digest(HERE/n) for n in artifacts},
              navigation_links=links, broken_links=0,
              published_navigation_hashes={p.relative_to(ROOT).as_posix(): core.digest(p) for p in paths},
              inherited_767_768_free_state_preserved=True, one_loop_single_source_Ward_inherited=True, interacting_BRST_lift_proven=False, active_goal_unchanged=True,
              all_checks_passed=True)
with TARGET.open("x", encoding="utf8") as f:
    json.dump(report, f, ensure_ascii=False, indent=2)
print(json.dumps({k: report[k] for k in ("latest_completed_round", "navigation_links", "all_checks_passed")}))
