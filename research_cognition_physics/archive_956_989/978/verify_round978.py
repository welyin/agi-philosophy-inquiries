"""Delivery checks for 978. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_978_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,978):
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
    result=read(HERE/"moving_complete_source_results.json")
    assert result["round"]==978 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["parameters"];b=result["uniform_moving_error"];a=result["operator_bounds"]
    assert p["G"]==1e-30 and p["reference_radius"]==1e10 and p["position_sigma"]==1e3
    assert p["mu"]>2e35
    assert a["internal_occupation_not_physically_cutoff"]
    assert a["initial_rms_mass_upper"]==204 and a["analytic_triangle_mass_upper"]<204
    assert a["weak_potential_global_upper"]<.003
    assert a["full_H_lower_excluding_mu"]>198
    assert b["total"]<.001721
    assert b["actual_moving_record_contrast_lower"]>.4825
    assert b["moving_output_to_977_cq_bound"]<.010383
    assert b["arbitrary_sender_and_passive_reference"] and b["interval_starts_at_zero"]
    assert b["tail_probability_upper"]>0 and b["minimum_tail_exponent"]>1e11
    c=result["conservation"]
    for key in ("total_momentum_conserved_exactly","total_energy_conserved_exactly",
        "full_internal_mass_conserved_exactly","both_centres_are_dynamical",
        "source_is_complete_internal_energy"):assert c[key]
    assert c["force_identity_difference"]<1e-12 and c["speed_norm_upper"]<1e-5
    for row in result["spatial_checks"]:
        assert row["quadrature_second_moment"]<row["analytic_second_moment_upper"]
        assert row["force_Y_is_exact_opposite"]
    for key in ("actual_wavepacket_evolution_computed","force_signal_above_measurement_error_claimed",
        "full_post_Newtonian_budget","dynamical_Einstein_metric_derived",
        "entire_SM_matching_completed","actual_recombination_device_completed","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core978",HERE/"moving_complete_source.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_978.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==12
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,7)]
    for term in ("整体目标未完成","不直接继承945的旧误差数值",
        "尚未计入实体装置及读后能源","不是完整Einstein动态时空",
        "引力作用是物理输入，不是认知原则生成的定理"):
        assert term in prose,term
    newfiles=[note,HERE/"moving_complete_source.py",HERE/"moving_complete_source_results.json",
        Path(__file__),HERE/"drafts/moving_source_decision.md",
        HERE/"drafts/common_recovery_audit_v2.md",HERE/"drafts/publish978.py",
        STAGE/"979/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至978）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=978)==list(range(1,979))
    assert "001—978轮共978份" in nav[0].read_text("utf-8-sig")
    assert "231—978的748份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "978：完整原生能源与双体运动" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"977/research_round_977_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=978,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=978,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_977=prev["cumulative_numbered_test_groups_from_976"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        native_full_energy_and_two_moving_positions_verified=True,
        uniform_moment_based_transport_bound_verified=True,
        common_momentum_energy_and_force_identity_verified=True,
        full_GR_or_post_Newtonian_matching_completed=False,
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

