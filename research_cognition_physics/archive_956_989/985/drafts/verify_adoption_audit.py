"""Adoption audit after 984. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parents[1]
ROOT=HERE.parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"common_model_adoption_audit_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,985):
        old=read(STAGE/f"{n}/research_round_{n}_checks.json")
        for key in ("frozen_inputs","new_scientific_and_entry_files"):add(old[key])
    for path in ("875/drafts/clock_transport_working_checks.json",
        "878/drafts/working_checks.json","884/drafts/working_checks.json",
        "887/drafts/working_checks.json","894/drafts/working_checks.json",
        "896/drafts/working_checks.json","897/drafts/working_checks.json",
        "899/drafts/working_checks.json","900/drafts/working_checks.json",
        "904/drafts/working_checks.json","907/drafts/working_checks.json",
        "909/drafts/working_checks.json","912/drafts/working_checks.json",
        "913/drafts/working_checks.json","915/drafts/working_checks.json",
        "917/drafts/working_checks.json"):
        add(read(STAGE/path)["files"])
    extra={
      "909/drafts/scope_reaudit_checks.json":("preserved_files","evidence_and_audit_hashes"),
      "911/drafts/finite_scope_checkpoint_checks.json":("evidence_and_audit_hashes",),
      "912/drafts/effective_scope_reaudit_checks.json":("evidence_and_audit_hashes",),
      "914/drafts/finite_scope_decision_checks.json":("evidence_hashes","new_document_and_verifier_hashes"),
      "919/drafts/effective_scope_after_918_checks.json":("frozen_inputs_verified","new_document_and_verifier_hashes"),
      "953/drafts/priority_reaudit_checks.json":("frozen_evidence_hashes","audit_file_hashes")}
    extra["981/drafts/common_adoption_audit_checks.json"]=("frozen_evidence_hashes","audit_file_hashes")
    for path,keys in extra.items():
        value=read(STAGE/path)
        for key in keys:add(value[key])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    last=read(STAGE/"984/research_round_984_checks.json")
    assert last["all_delivery_checks_passed"] and last["formal_reports"]==984
    assert last["cumulative_numbered_test_groups_from_983"]==3769
    assert not (STAGE/"research_note_985.md").exists()
    r=read(STAGE/"984/smeared_source_geometry_results.json")
    assert r["scope"]["matter_induced_linear_Einstein_report_only"]
    assert not r["scope"]["finite_spacetime_detector_certified"]
    assert not r["scope"]["full_quantum_geometry_error_certified"]
    paths=[STAGE/f"research_note_{n}.md" for n in
           (898,936,939,940,950,951,953,957,961,964,969,971,973,978,979,980,981,982,983,984)]
    paths += [STAGE/"957/drafts/unified_operation_hypotheses_v0_2.md",
        STAGE/"951/drafts/joint_minimum_menu.md",STAGE/"981/drafts/common_parent_contract_v1.md",
        STAGE/"983/drafts/parent_recovery_scope_v2.md",STAGE/"984/drafts/geometry_adoption_scope.md",
        STAGE/"962/drafts/revised_goal_20261007.txt",STAGE/"962/drafts/app_goal_confirmation.json",
        STAGE/"984/research_round_984_checks.json",HERE/"STATUS.md"]
    for p in paths:assert p.exists(),p
    note=HERE/"research_note_985_working.md"
    s=note.read_text("utf-8")
    for term in ("不是新科学轮次","正式984／累计3769保持","U1：同一实际物理生成过程",
        "U2：实际来源与几何报告","U3：共同字典及采用域","不再接续984的噪声核细化",
        "不是新认知公理"):
        assert term in s,term
    assert "作用的每个量子／经典投影已匹配" in s
    audit=[note,HERE/"publish_adoption_audit.py",Path(__file__)]
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    links=0
    for p in audit:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    for p in [note]+nav:
        text=p.read_text("utf-8-sig")
        if p in nav:assert "984后整体采用审计：先接共同物理过程" in text
        content=re.sub(r"\$\$.*?\$\$","",text,flags=re.S)
        for link in re.findall(r"\]\(([^)]+)\)",content):
            if re.match(r"^[a-zA-Z]+://",link) or link.startswith("#"):continue
            target=(p.parent/link.split("#")[0].strip("<>")).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(p,link)
            links+=1
    assert "001—984轮共984份" in nav[0].read_text("utf-8-sig")
    assert "231—984的754份" in nav[5].read_text("utf-8-sig")
    ledger=nav[-1].read_text("utf-8-sig")
    section=ledger.split("## 六条共同协议：全局缺口对应与检验优先级（截至984）",1)[1].split("### 当前取舍",1)[0]
    assert re.findall(r"^\|(C\d\d) ",section,re.M)==[f"C{i:02d}" for i in range(1,28)]
    out=dict(date="2026-10-07",all_audit_checks_passed=True,formal_reports=984,
        cumulative_numbered_test_groups=3769,new_scientific_test_groups=0,
        previous_goal_turn="progress",current_audit="adoption_and_priority_decision",
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,full_goal_completed=False,app_goal_changed=False,
        next_focus="U1: P981 physical process and common task dictionary, before U2 classical geometry",
        frozen_evidence_hashes={str(p.relative_to(ROOT)):sha(p) for p in paths},
        audit_file_hashes={str(p.relative_to(ROOT)):sha(p) for p in audit})
    if not writing:
        saved=read(TARGET)
        for key in ("frozen_evidence_hashes","audit_file_hashes"):assert out[key]==saved[key],key
    return out

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    out=run(args.write)
    if args.write:
        with TARGET.open("x",encoding="utf-8") as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write("\n")
    print(json.dumps({k:v for k,v in out.items() if k not in
       ("frozen_evidence_hashes","audit_file_hashes")},ensure_ascii=False,indent=2))
