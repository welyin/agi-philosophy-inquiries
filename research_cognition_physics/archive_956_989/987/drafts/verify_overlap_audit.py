"""Overlap audit after 986. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parents[1]
ROOT=HERE.parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"common_overlap_audit_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,987):
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
    extra["985/drafts/common_model_adoption_audit_checks.json"]=("frozen_evidence_hashes","audit_file_hashes")
    for path,keys in extra.items():
        value=read(STAGE/path)
        for key in keys:add(value[key])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    import importlib.util
    from fractions import Fraction
    spec=importlib.util.spec_from_file_location("audit987",HERE.parent/"overlap_adoption.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    result=read(HERE.parent/"overlap_adoption_results.json")
    mod.compare(mod.run(),result)
    for path,digest in result["source_hashes"].items():assert sha(ROOT/path)==digest,path
    assert Fraction(result["common_record"]["state_trace_distance_upper"]["exact"])==Fraction(860567,500000000)
    assert result["conditional_profile"]["k_sigma"]==20
    for row in result["conditional_profile"]["linear_source_second_moment_by_lapse"]:
        for constant,profiled,difference in zip(
            row["linear_vertex_force_second_moment_constant"],
            row["linear_vertex_force_second_moment_profiled"],row["absolute_difference"]):
            assert constant>1e-29 and 6e-30<profiled<7e-30
            assert math.isclose(profiled/constant,.5,rel_tol=1e-13,abs_tol=0)
            assert math.isclose(difference,constant-profiled,rel_tol=1e-13,abs_tol=0)
    assert not result["spectator_comparison"]["position_kinetic_energy_included"]
    assert not result["spectator_comparison"]["proves_join_with_985_motion"]
    assert not result["decision"]["actual_finite_record_failure_demonstrated"]
    assert not result["decision"]["full_goal_completed"]
    assert min(result["established_EFT"]["full_binary_k_r"])>9e7
    last=read(STAGE/"986/research_round_986_checks.json")
    assert last["all_delivery_checks_passed"] and last["formal_reports"]==986
    assert last["cumulative_numbered_test_groups_from_985"]==3771
    assert not (STAGE/"research_note_987.md").exists()
    note=HERE/"research_note_987_working.md"
    prose=note.read_text("utf-8")
    for text in ("不新增正式科学轮次","正式986／累计3771保持",
        "不是完整非线性力的全时间二阶矩","不增加微观连续或物理最小尺度假设",
        "不是985运动已被证明可忽略于完整合并模型","不是父理论的总物理误差"):
        assert text in prose,text
    assert prose.count("$$")==10
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,6)]
    audit=[note,HERE/"overlap_decision.md",HERE/"publish_overlap_audit.py",Path(__file__),
        HERE.parent/"overlap_adoption.py",HERE.parent/"overlap_adoption_results.json"]
    paths=[STAGE/"986/research_round_986_checks.json",HERE/"STATUS.md",
        STAGE/"985/drafts/research_note_985_working.md",
        STAGE/"962/drafts/revised_goal_20261007.txt",
        STAGE/"962/drafts/app_goal_confirmation.json"]
    for path in result["source_hashes"]:
        p=ROOT/path
        if p not in audit and p not in paths:paths.append(p)
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    links=0
    for p in audit:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    for p in [a for a in audit if a.suffix==".md"]+nav:
        text=p.read_text("utf-8-sig")
        if p in nav:assert "986后合并审计：记录交叠与来源边界" in text
        content=re.sub(r"\$\$.*?\$\$","",text,flags=re.S)
        for link in re.findall(r"\]\(([^)]+)\)",content):
            if re.match(r"^[a-zA-Z]+://",link) or link.startswith("#"):continue
            target=(p.parent/link.split("#")[0].strip("<>")).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(p,link)
            links+=1
    assert "001—986轮共986份" in nav[0].read_text("utf-8-sig")
    assert "231—986的756份" in nav[5].read_text("utf-8-sig")
    ledger=nav[-1].read_text("utf-8-sig")
    section=ledger.split("## 六条共同协议：全局缺口对应与检验优先级（截至986）",1)[1].split("### 当前取舍",1)[0]
    assert re.findall(r"^\|(C\d\d) ",section,re.M)==[f"C{i:02d}" for i in range(1,28)]
    out=dict(date="2026-10-07",all_audit_checks_passed=True,formal_reports=986,
        cumulative_numbered_test_groups=3771,new_scientific_test_groups=0,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,full_goal_completed=False,app_goal_changed=False,
        full_985_986_joint_model_claimed=False,
        next_focus="P981 common retained-order source and task dictionary; no mode optimization",
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
