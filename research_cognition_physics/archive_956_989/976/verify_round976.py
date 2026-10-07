"""Delivery checks for 976. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_976_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,976):
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
    result=read(HERE/"finite_internal_reference_results.json")
    assert result["round"]==976 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    contract=result["fixed_relational_readout"];cert=result["interval_certificate"]
    assert contract["effect"]=="(I+SWAP_BR)/2"
    assert contract["free_common_invariance_residual"]<1e-14
    assert contract["reference_preparation_is_resource"] and contract["isolation_is_additional_input"]
    assert contract["external_time_dependent_effect"] is False
    assert cert["whole_interval_not_only_grid_certified"]
    assert cert["uniform_infinite_contrast_lower"]>.48594
    assert cert["grid_points"]==257 and cert["mesh_gap_bound"]<.000956
    assert cert["binary_information_lower"]>.1232
    assert result["infinite_occupation_certificate"]["total_state_vector_error"]<2.733e-8
    account=result["reference_account"]
    assert account["source_includes_complete_reference_energy"]
    assert account["source_mass_positive_lower"]>199.68
    assert account["old_endpoint_check"]<1e-9
    assert account["old_path_predictions_unchanged"] is False
    assert account["internal_energy_variance"]>.00037
    for row in account["source_rows"]:
        assert row["old_path_visibility"]-row["new_path_visibility"]>.0012
        assert row["mass_mean_with_reference"]>200
    for key in ("phase_record_converted_to_protected_population","actual_joint_apparatus_derived",
        "arbitrary_noise_stability_proven","repeatable_nondestructive_readout_proven",
        "macroscopic_irreversible_memory_proven","full_microscopic_SM_matching_proven",
        "moving_support_or_quantum_metric_completed","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core976",HERE/"finite_internal_reference.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_976.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==14
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,8)]
    for term in ("整体目标未完成","隔离及共同支撑是输入","不是新的普遍认知公理",
        "尚未给具体仪器、读后能源账或复位","不得把其器件缺口升为全纲领门槛"):
        assert term in prose,term
    newfiles=[note,HERE/"finite_internal_reference.py",HERE/"finite_internal_reference_results.json",
        Path(__file__),HERE/"drafts/finite_reference_decision.md",
        HERE/"drafts/finite_record_adoption.md",HERE/"drafts/publish976.py",
        STAGE/"977/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至976）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=976)==list(range(1,977))
    assert "001—976轮共976份" in nav[0].read_text("utf-8-sig")
    assert "231—976的746份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "976：有限内部参考与关系记录" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"975/research_round_975_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=976,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=976,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_975=prev["cumulative_numbered_test_groups_from_974"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        finite_internal_reference_and_whole_window_verified=True,
        reference_mass_and_source_predictions_updated=True,
        macroscopic_irreversible_memory_proven=False,
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

