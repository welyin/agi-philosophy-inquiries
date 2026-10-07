"""Delivery checks for 971. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_971_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,971):
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
    result=read(HERE/"native_covariant_source_results.json")
    assert result["round"]==971 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["fixed_native_inputs"];r=result["recoil"];c=result["covariance"]
    assert math.isclose(p["electric_dipole_matrix_element"]**2,p["d"],rel_tol=1e-13)
    assert math.isclose(p["electronic_mass_gap"],p["excited_mass"]-p["ground_mass"],abs_tol=1e-13)
    assert math.isclose(r["no_recoil_mass_shell_residual"],p["electronic_mass_gap"]**2,abs_tol=5e-12)
    assert r["recoil_shift_over_benchmark_width"]>120
    assert r["bare_mass_recoil_difference"]<.05*r["benchmark_970_half_width"]
    assert r["benchmark_is_not_a_waveguide_to_vacuum_identification"]
    assert c["maximum_lorentz_vertex_error"]<1e-13
    assert c["maximum_gauge_shift_error"]<1e-13
    assert c["maximum_ward_residual"]<1e-13
    assert c["maximum_mass_shell_residual"]<2e-11
    assert math.isclose(c["rows"][1]["electric_only_relative_error"],.5,abs_tol=1e-12)
    m=result["matching"]
    assert abs(m["actual_CAR_energy_curvature"]-m["resolved_polarizability"])<1e-7
    assert abs(m["wrong_explicit_plus_full_contact_curvature"]-2*m["resolved_polarizability"])<1e-7
    assert not m["self_polarization_of_965_removed"]
    for key in ("actual_transition_rates_certified","real_rotational_selection_rules_matched",
                "full_SM_to_material_matching","full_quantum_positive_parent_constructed",
                "round970_disproved","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core971",HERE/"native_covariant_source.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_971.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","970没有被推翻","不删除965的自极化项",
                 "运动学必要条件","不是已经认证的完整共同模型"):
        assert term in prose,term
    newfiles=[note,HERE/"native_covariant_source.py",HERE/"native_covariant_source_results.json",
        Path(__file__),HERE/"drafts/native_parent_decision.md",
        HERE/"drafts/native_parent_map_v0_1.md",HERE/"drafts/publish971.py",
        STAGE/"972/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至971）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=971)==list(range(1,972))
    assert "001—971轮共971份" in nav[0].read_text("utf-8-sig")
    assert "231—971的741份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "971：原生材料的共同父接口与重复计数检验" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"970/research_round_970_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=971,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=971,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_970=prev["cumulative_numbered_test_groups_from_969"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        same_native_mass_gap_and_electric_vertex=True,
        covariant_vertex_and_kinematic_recoil_verified=True,
        finite_matching_double_counting_counterexample=True,
        real_transition_rates_verified=False,
        full_quantum_parent_verified=False,full_SM_matching_verified=False,
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

