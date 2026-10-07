"""Delivery checks for 975. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_975_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,975):
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
    result=read(HERE/"record_resource_audit_results.json")
    assert result["round"]==975 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    assert result["unchanged_physical_hamiltonian_preparation_and_endpoint"]
    ref=result["reference"];bal=result["resource_balance"];rec=result["record"]
    assert abs(ref["beta"]-2*math.log(2))<1e-14
    assert ref["reference_is_actual_joint_equilibrium"] is False
    assert ref["physical_heat_bath_added"] is False
    assert 6.78<bal["initial_readiness"]<6.79
    assert 5.25<bal["final_local_readiness"]<5.26
    assert 1.53<bal["final_total_correlation"]<1.531
    assert abs(bal["identity_residual"])<1e-10
    assert bal["interaction_energy_change"]<-1.6e-5
    assert rec["certified_binary_information_lower"]>.67507
    assert rec["sender_retains_label_exactly"]
    assert rec["arbitrary_unknown_quantum_state_copied"] is False
    assert max(abs(x-.5) for x in rec["receiver_conserved_sector_probabilities"])<1e-11
    assert rec["after_sector_dephasing_infinite_distance_upper"]<.000106
    assert rec["sector_dephasing_is_diagnostic_only"]
    errs=result["error_bounds"]
    assert errs["oscillator_mean_number_bound"]==1
    assert errs["derived_sqrt_mean_number_bound"]<.797
    assert errs["total_correlation"]<.000183 and errs["local_readiness"]<.000184
    assert errs["infinite_oscillator_entropy_controlled"]
    for key in ("global_entropy_production","actual_heat_or_work_cost_identified",
                "stable_macroscopic_record_dynamics_completed","full_cyclic_reset_completed",
                "cosmic_low_entropy_origin_derived","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core975",HERE/"record_resource_audit.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_975.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","乘积τ不是完整相互作用H的平衡态",
                 "没有把此退相位当成已经发生的环境噪声","不改变974的H、准备",
                 "新记录尚未写入现有作用所保护的变量"):
        assert term in prose,term
    newfiles=[note,HERE/"record_resource_audit.py",HERE/"record_resource_audit_results.json",
        Path(__file__),HERE/"drafts/record_resource_decision.md",
        HERE/"drafts/record_arrow_adoption.md",HERE/"drafts/publish975.py",
        STAGE/"976/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至975）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=975)==list(range(1,976))
    assert "001—975轮共975份" in nav[0].read_text("utf-8-sig")
    assert "231—975的745份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "975：原生记录的资源账与保持接口" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"974/research_round_974_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=975,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=975,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_974=prev["cumulative_numbered_test_groups_from_973"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        unchanged_joint_process_resource_balance_verified=True,
        same_record_distinguished_from_protected_sector=True,
        infinite_oscillator_entropy_error_controlled=True,
        global_entropy_increase_or_macroscopic_arrow_proven=False,
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

