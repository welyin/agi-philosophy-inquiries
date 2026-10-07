"""Delivery checks for 968. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_968_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,968):
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
    result=read(HERE/"internal_relay_results.json")
    assert result["round"]==968 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["parameters"]
    assert (p["physical_total_dimension"],p["exact_invariant_total_dimension"])==(13824,64)
    assert p["eta"]==.05 and p["kappa_AC"]==0.
    assert max(result["exact_restriction_checks"].values())<1e-14
    assert math.isclose(result["relay_initial_participation_acceleration"],3*.05**2/8)
    assert result["endpoint_record_projectors_conserved"]
    dark,active=result["fixed_menu_rows"]
    assert dark["contrast"]<1e-8 and not dark["accepted_nonzero_window"]
    assert active["window_contrast_lower"]>.338 and active["accepted_nonzero_window"]
    assert active["relay_participation_current_norm"]>.02
    assert abs(active["details"][0]["relay_spin_singlet_probability"]-
               active["details"][1]["relay_spin_singlet_probability"])>.068
    for row in (dark,active):
        assert row["certificate"]["total_operator_error_bound"]<1e-4
        assert row["energy_current_balance"]<1e-12
        for d in row["details"]:
            assert abs(sum(d["before_energies"])-sum(d["after_energies"]))<1e-11
    assert all(c["exact_contrast"]==0. for c in result["edge_cut_controls"])
    assert result["scope"]["signal_window_adopted"]
    for key in ("external_time_dependent_coupling_used","new_graph_register_added",
                "geometry_change_derived","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core968",HERE/"internal_relay.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_968.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","不是低能截断","固定接触布局仍然是输入",
                 "标准舍入模型","固定功能关系","时间箭头未恢复"):
        assert term in prose,term
    newfiles=[note,HERE/"internal_relay.py",HERE/"internal_relay_results.json",
        Path(__file__),HERE/"drafts/internal_relation_decision.md",
        HERE/"drafts/internal_organization_increment.md",HERE/"drafts/publish968.py",
        STAGE/"969/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至968）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=968)==list(range(1,969))
    assert "001—968轮共968份" in nav[0].read_text("utf-8-sig")
    assert "231—968的738份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "968：原生材料的内部组织更新与自主传递" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"967/research_round_967_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=968,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=968,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_967=prev["cumulative_numbered_test_groups_from_966"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        exact_invariant_material_sector=True,
        autonomous_internal_organization_updates=True,
        fixed_effect_nonzero_window_certified=True,
        numerical_certificate_scope="standard rounding model; not formal library verification",
        no_external_time_dependent_contact=True,
        same_H_energy_and_parameter_sources=True,
        geometry_change_derived=False,full_goal_completed=False,
        visual_checks_performed=False,app_goal_changed=False,
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

