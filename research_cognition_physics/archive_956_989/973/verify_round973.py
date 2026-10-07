"""Delivery checks for 973. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_973_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,973):
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
    result=read(HERE/"common_source_distribution_results.json")
    assert result["round"]==973 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    inputs=result["frozen_inputs"];stats=result["source_statistics"]
    assert inputs["modes"]==12 and inputs["mean_R"]==100
    assert .289<stats["coefficient_of_variation"]<.290
    assert .144<stats["initial_H_squared_relative_std"]<.145
    assert abs(stats["exact_variance"]-stats["numerical_variance"])<1e-8
    summ=result["summation"];diag=result["monopole_branch_diagnostic"]
    assert summ["chernoff_tail"]<6e-25
    assert abs(summ["normalization"]-1)<2e-12
    assert diag["report_certified_lower"]>.6725
    assert diag["report_under_mean_source"]==0
    assert diag["endpoint_time_margins"][0]>6 and diag["endpoint_time_margins"][1]>1
    assert diag["constraint_residual_at_mean_initial_momentum"]>1
    assert .0008<diag["mean_redshift_difference"]<.001
    assert result["mixture_witness"]["mixture_report"]==1
    assert result["mixture_witness"]["mean_source_report"]==0
    assert result["mixture_witness"]["initial_geometry_correlated_with_energy"]
    for key in ("full_stress_noise_or_Einstein_Langevin_solved",
                "original_969_mean_field_result_refuted","all_unknown_quantum_inputs_supported",
                "actual_geometric_detector_constructed","full_common_positive_quantum_model","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core973",HERE/"common_source_distribution.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_973.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","不是完整应力噪声或Einstein–Langevin解",
                 "保留969原结论","共同理论允许不同任务使用不同准备和边界",
                 "能源与初始几何相关的准备"):
        assert term in prose,term
    newfiles=[note,HERE/"common_source_distribution.py",HERE/"common_source_distribution_results.json",
        Path(__file__),HERE/"drafts/common_realization_decision.md",
        HERE/"drafts/common_realization_audit_v1.md",HERE/"drafts/publish973.py",
        STAGE/"974/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至973）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=973)==list(range(1,974))
    assert "001—973轮共973份" in nav[0].read_text("utf-8-sig")
    assert "231—973的743份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "973：共同实现的准备审计与平均来源边界" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"972/research_round_972_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=973,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=973,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_972=prev["cumulative_numbered_test_groups_from_971"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        actual_frozen_thermal_source_distribution_verified=True,
        bounded_report_in_monopole_diagnostic=True,
        source_geometry_constraint_correlation_retained=True,
        global_common_realization_audited=True,
        original_mean_field_refuted=False,full_stochastic_gravity_solved=False,
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

