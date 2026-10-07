"""Delivery checks for 980. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_980_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,980):
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
    result=read(HERE/"finite_thermal_records_results.json")
    assert result["round"]==980 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["parameters"];c=result["channel"];r=result["finite_resources"]
    assert p["eta"]==.05 and p["theta"]==1 and p["Q_star_center"]==2
    assert c["contacts"]==600000 and c["full_charge_reset_bound"]<.001029
    assert c["comparison_thermal_reset_bound"]==13/14000
    assert c["arbitrary_unknown_input_and_passive_reference"]
    assert c["one_contact_dilation_error"]<1e-12 and c["stationary_error"]<1e-12
    assert c["final_choi_reset_bound"]<c["comparison_thermal_reset_bound"]
    cert=result["eigensystem_certificate"]["coarse_analytic_certificate"]
    assert cert["analytic_spectral_gap_lower"]>1/62500
    assert cert["nonresonant_gap_rational_lower"]>.03
    assert cert["minimum_population_diagonal_lower"]>.5
    assert cert["exp_9_6_rational_lower"]>14000
    assert r["prepared_thermal_copies"]==600000 and r["no_bath_reused_or_reset"]
    assert r["total_contact_time"]>5e19 and r["initial_bath_energy_above_ground"]>1e5
    assert r["thermal_state_to_best_pure_blank_distance"]>.58
    assert r["switching_work_absolute_bound"]<6.3e-9
    assert r["controller_preparation_not_constructed"]
    therm=result["thermodynamic_check"]
    assert therm["energy_balance_error"]<1e-12 and therm["landauer_identity_error"]<1e-12
    assert therm["endpoint_entropy_production"]>0
    for row in result["averaging_diagnostics"]:assert row["full_unitary_error"]<row["averaging_bound"]
    for key in ("infinite_bath_required","new_primitive_reset_operator_inserted",
        "constant_autonomous_controller_constructed","pure_blank_reset_proved",
        "universe_permanent_arrow_proved","actual_geometric_feedback_computed",
        "whole_lifecycle_with_976_verified","full_SM_matching_completed","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core980",HERE/"finite_thermal_records.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_980.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","不是宇宙永久时间箭头，也不是纯空白复位",
        "小开关平均功不是小控制成本证明","不保证实际每一步严格单调",
        "自主控制器未构造"):
        assert term in prose,term
    newfiles=[note,HERE/"finite_thermal_records.py",HERE/"finite_thermal_records_results.json",
        Path(__file__),HERE/"drafts/finite_thermal_decision.md",
        HERE/"drafts/thermal_arrow_adoption.md",HERE/"drafts/publish980.py",
        STAGE/"981/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至980）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=980)==list(range(1,981))
    assert "001—980轮共980份" in nav[0].read_text("utf-8-sig")
    assert "231—980的750份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "980：有限热材料与有效遗忘" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"979/research_round_979_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=980,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=980,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_979=prev["cumulative_numbered_test_groups_from_978"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        finite_native_thermal_reset_with_reference_verified=True,
        analytic_contraction_and_finite_averaging_verified=True,
        full_bath_and_mean_switching_work_accounts_verified=True,
        autonomous_controller_or_pure_blank_reset_completed=False,
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

