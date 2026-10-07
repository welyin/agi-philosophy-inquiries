"""Delivery checks for 977. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_977_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,977):
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
    result=read(HERE/"coarse_source_record_results.json")
    assert result["round"]==977 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["parameters"];b=result["uniform_budget"];a=result["algebra"]
    assert p["bin_width"]==.01 and p["joint_internal_dimension"]==352
    assert p["minimal_finite_spectrum_distance_to_bin_edge"]>2e-4
    assert b["uniform_on_entire_interval"]
    assert b["maximum_joint_output_distance"]<.008663
    assert b["record_contrast_lower"]>.46862
    assert b["infinite_occupation_certificate"]["total_state_vector_error"]<2.711e-8
    assert a["source_derivative_error"]<.0025
    assert a["full_energy_commutes_with_cq_generator_exactly"]
    assert a["finite_spectral_polar_basis_defines_positive_candidate"]
    assert a["minimum_cq_generator"]>199.67
    for row in result["rows"]:
        assert row["same_energy_mean_and_distribution_under_pinching"]
        assert row["pinching_menu_uniform_distance_bound"]<1.335e-5
        for endpoint in row["endpoints"]:
            assert endpoint["actual_total_joint_distance"]<row["full_interval_joint_distance_bound"]
    assert min(result["cq_endpoint_record_contrasts"])>.48
    assert max(result["rank_one_energy_diagonal_endpoint_record_contrasts"])<3.04e-6
    for row in result["characteristic_functions"]:
        assert .704<row["coarse_visibility"]<.705
        assert row["characteristic_error"]<math.pi*.01/2
    scope=result["scope"]
    assert scope["classical_bin_and_quantum_fibre_adopted"]
    assert scope["arbitrary_sender_reference_on_declared_joint_output"]
    for key in ("pinching_claimed_as_physical_irreversibility","all_internal_quantum_information_preserved",
        "exact_infinite_H_spectral_projectors_approximated","repeat_measurement_or_future_control_covered",
        "full_source_replaced_by_single_deterministic_mean","dynamical_classical_metric_derived",
        "Einstein_Langevin_noise_kernel_matched","full_goal_completed"):
        assert scope[key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core977",HERE/"coarse_source_record.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_977.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==14
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,8)]
    for term in ("整体目标未完成","不是声称原宇宙实际执行了谱测量或退相干",
        "并没有声称已逼近无限H的不连续谱投影" if False else "没有声称已逼近无限H的不连续谱投影",
        "本轮没有声称恢复Einstein—Langevin方程","宽度与时间没有在看到结果后调整"):
        assert term in prose,term
    newfiles=[note,HERE/"coarse_source_record.py",HERE/"coarse_source_record_results.json",
        Path(__file__),HERE/"drafts/coarse_source_decision.md",
        HERE/"drafts/coarse_source_adoption.md",HERE/"drafts/publish977.py",
        STAGE/"978/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至977）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=977)==list(range(1,978))
    assert "001—977轮共977份" in nav[0].read_text("utf-8-sig")
    assert "231—977的747份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "977：粗来源标签与块内量子记录" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"976/research_round_976_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=977,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=977,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_976=prev["cumulative_numbered_test_groups_from_975"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        common_record_and_path_output_bound_verified=True,
        source_derivative_error_and_energy_conservation_verified=True,
        classical_bin_quantum_fibre_distinguished_from_full_dephasing=True,
        dynamical_classical_metric_derived=False,
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

