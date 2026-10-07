"""Delivery checks for 983. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_983_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,983):
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
    result=read(HERE/"sm_curvature_radiation_results.json")
    assert result["round"]==983 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    m=result["inherited_matter"];a=result["analytic"];b=result["benchmark"]
    assert (m["real_scalars"],m["Weyl_components"],m["gauge_vectors"])==(4,45,12)
    assert m["a"]["exact"]=="1991/720" and m["c"]["exact"]=="283/120"
    assert m["thermal_gstar"]["exact"]=="427/4"
    assert m["beta_curvature_times_16pi2"]==dict(bR="0",bW="283/120",bE="-1991/720")
    assert a["wrong_radiation_pressure_conservation_coefficient"]=="-4"
    assert a["missing_pressure_conservation_coefficient"]=="-5"
    assert a["proper_time_error_upper"]["exact"]=="63/25600000000"
    assert 0<b["reduced_vs_low_time_difference"]<a["proper_time_error_upper"]["decimal"]
    assert b["shared_equations_max_residual"]<2e-14
    assert float(b["branches"]["z_times_high"])>.98
    for key in ("full_interacting_SM_thermal_state","all_quantum_source_to_classical_geometry_error",
        "full_physical_truncation_error_certified","high_curvature_root_adopted",
        "cosmic_arrow_generated","full_goal_completed"):assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core983",HERE/"sm_curvature_radiation.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_983.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","不重新算成新发现","β_bR=0不表示有限bR=0",
        "不是全部物理遗漏的上界","同一份量子物质修正", "停止反常宇宙分支"):
        assert term in prose,term
    newfiles=[note,HERE/"sm_curvature_radiation.py",HERE/"sm_curvature_radiation_results.json",
        Path(__file__),HERE/"drafts/curvature_adoption_decision.md",
        HERE/"drafts/parent_recovery_scope_v2.md",HERE/"drafts/publish983.py",STAGE/"984/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至983）",1)[1].split("### 当前取舍",1)[0]
    assert re.findall(r"^\|(C\d\d) ",mapping,re.M)==[f"C{i:02d}" for i in range(1,28)]
    links=0
    for doc in [p for p in newfiles if p.suffix==".md"]+nav:
        content=re.sub(r"\$\$.*?\$\$","",doc.read_text("utf-8-sig"),flags=re.S)
        for link in re.findall(r"\]\(([^)]+)\)",content):
            if re.match(r"^[a-zA-Z]+://",link) or link.startswith("#"):continue
            target=(doc.parent/link.split("#")[0].strip("<>")).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for p in STAGE.parent.rglob("research_note_*.md"):
        match=re.fullmatch(r"research_note_(\d+).md",p.name)
        if match and p.parent.name.startswith("archive_"):nums.append(int(match[1]))
    assert sorted(n for n in nums if n<=983)==list(range(1,984))
    assert "001—983轮共983份" in nav[0].read_text("utf-8-sig")
    assert "231—983的753份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "983：标准模型曲率来源与宇宙有效阶" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"982/research_round_982_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=983,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=983,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_982=prev["cumulative_numbered_test_groups_from_981"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        inherited_SM_curvature_coefficients_mapped=True,
        same_source_density_pressure_trace_and_geometry_verified=True,
        finite_low_branch_time_bound_verified=True,
        full_physical_error_certified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
    if not writing:
        before=read(TARGET)
        for key in ("frozen_inputs","new_scientific_and_entry_files"):assert before[key]==out[key],key
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
         ("frozen_inputs","new_scientific_and_entry_files")},ensure_ascii=False,indent=2))

