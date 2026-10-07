"""Delivery checks for 972. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_972_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,972):
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
    for path,keys in extra.items():
        value=read(STAGE/path)
        for key in keys:add(value[key])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/"native_thermal_record_results.json")
    assert result["round"]==972 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    assert [x["commutant_dimension"] for x in result["algebra"]]==[2,1]
    assert result["algebra"][0]["sector_commutator"]==0
    assert result["algebra"][1]["sector_commutator"]>.02
    for row in result["rows"]:
        assert abs(sum(row["energy_changes"]))<1e-11
        assert abs(sum(row["energy_currents"]))<1e-12
        assert row["landauer_identity_error"]<1e-12
        assert row["entropy_production"]>row["entropy_production_error_bound"]
        assert row["certificate"]["material_trace_distance_bound"]<4e-8
        if row["eta"]==0:assert abs(row["record_information_loss"])<1e-12
    assert result["witness"]["information_loss_lower"]>2.15e-5
    assert result["rows"][3]["field_heat"]+result["rows"][3]["heat_error_bound"]<0
    assert abs(result["rows"][-1]["field_heat"])<result["rows"][-1]["heat_error_bound"]
    assert abs(result["isolated_material_control_information"]-math.log(2))<1e-12
    for key in ("complete_reset_proved","thermalization_proved",
                "monotone_irreversible_arrow_proved","cosmology_common_window_proved","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core972",HERE/"native_thermal_record.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_972.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","不领取正热流结论","不是材料熵减少的擦除过程",
                 "不是整个24维材料没有守恒量","没有证明完整复位"):
        assert term in prose,term
    newfiles=[note,HERE/"native_thermal_record.py",HERE/"native_thermal_record_results.json",
        Path(__file__),HERE/"drafts/thermal_record_decision.md",HERE/"drafts/probe_contact.py",
        HERE/"drafts/thermal_record_mechanism.md",HERE/"drafts/publish972.py",
        STAGE/"973/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至972）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=972)==list(range(1,973))
    assert "001—972轮共972份" in nav[0].read_text("utf-8-sig")
    assert "231—972的742份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "972：原生记录的热接触、保护扇区与内部交换" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"971/research_round_971_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=972,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=972,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_971=prev["cumulative_numbered_test_groups_from_970"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        arbitrary_electric_bath_sector_obstruction=True,
        old_native_exchange_lifts_local_protection=True,
        common_infinite_mode_finite_window_error_bound=True,
        finite_information_energy_entropy_checked=True,
        complete_reset_verified=False,irreversible_arrow_verified=False,
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

