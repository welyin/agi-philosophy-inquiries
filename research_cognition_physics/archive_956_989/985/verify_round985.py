"""Delivery checks for 985. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_985_checks.json"
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
    extra["985/drafts/common_model_adoption_audit_checks.json"]=("frozen_evidence_hashes","audit_file_hashes")
    for path,keys in extra.items():
        value=read(STAGE/path)
        for key in keys:add(value[key])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/"two_active_sources_results.json")
    assert result["round"]==985 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    import importlib.util
    from fractions import Fraction
    spec=importlib.util.spec_from_file_location("core985",HERE/"two_active_sources.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    budget=result["uniform_budget"]
    assert Fraction(budget["exact_sum_upper"]["exact"])<Fraction(1721,10**6)
    assert Fraction(budget["record_contrast_lower"]["exact"])>Fraction(482506,10**6)
    assert result["initial_moments"]["unknown_joint_two_sender_input_and_passive_reference"]
    assert result["initial_moments"]["all_Fock_occupations_retained_in_theorem"]
    assert result["stability"]["relative_negative_form_bound"]<.00125
    for key in ("entire_U1_completed","radiative_TT_matched","full_SM_matching",
                "GR_total_physical_error_certified","physical_minimum_scale_assumed",
                "microscopic_continuity_assumed","UV_completion_claimed",
                "unbounded_source_error_deduced_from_probability","full_goal_completed"):
        assert result["scope"][key] is False
    note=STAGE/"research_note_985.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","停止双体扩建和高能延拓优化",
        "能量无下界不等于量子态不正或不存在酉演化","本轮误差不是距完整GR的误差"):
        assert term in prose,term
    newfiles=[note,HERE/"two_active_sources.py",HERE/"two_active_sources_results.json",
        Path(__file__),HERE/"drafts/two_active_sources_decision.md",HERE/"drafts/two_source_adoption.md",
        HERE/"drafts/publish985.py",STAGE/"986/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至985）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=985)==list(range(1,986))
    assert "001—985轮共985份" in nav[0].read_text("utf-8-sig")
    assert "231—985的755份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "985：两个活动来源与有限域共同过程" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"984/research_round_984_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=985,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=985,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_984=prev["cumulative_numbered_test_groups_from_983"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,two_active_native_sources_joined=True,
        finite_domain_extension_and_record_bounds_verified=True,
        entire_U1_completed=False,full_physical_error_certified=False,
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

