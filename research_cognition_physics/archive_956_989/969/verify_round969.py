"""Delivery checks for 969. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_969_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,969):
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
    result=read(HERE/"local_field_cosmology_results.json")
    assert result["round"]==969 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    prior=read(STAGE/"964/material_cosmology_results.json")
    assert result["shared_background_parameters"]["mu"]==prior["cosmology"]["mu"]
    assert result["shared_background_parameters"]["V0"]==prior["cosmology"]["coordinate_volume"]
    assert .25<result["local_sources"]["mass_increment"]<.25001
    assert result["local_sources"]["minimum_mass_lower"]>99
    assert result["local_sources"]["energy_change_error"]<1e-12
    bg=result["background"]
    assert 2.0006<bg["scale_at_T"]<2.0007
    assert bg["integration_error"]<2e-9 and bg["Raychaudhuri_constraint_max"]<2e-10
    joint=result["joint_task_transport"]
    assert joint["exact_factorization"] and joint["infinite_occupation_primary"]
    assert joint["inherited_receiver_contrast_lower"]>.9987
    assert joint["inherited_polarization_effect_contrast_lower"]>.0498
    assert joint["global_entropy_initial"]==joint["global_entropy_final"]
    bad=result["wrong_dictionary"]
    assert bad["nonzero_t_squared_coefficient"]<0
    assert bad["source_energy_offset_at_T"]<-.12
    for row in bad["initial_field_derivatives"]:
        assert abs(row["first"])<1e-15 and row["second"]>3.2e-7
    assert result["scope"]["one_nonselected_mean_background"]
    assert not any(result["scope"][key] for key in (
        "microscopic_cavity_matching","microscopic_Einstein_derivation",
        "external_radiation_absorption_recovered","irreversible_arrow_derived","full_goal_completed"))
    import importlib.util
    spec=importlib.util.spec_from_file_location("core969",HERE/"local_field_cosmology.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_969.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","不重调原引力系数","不是完整QED",
                 "精确分解","有限尺度","只反证上述混用","未完成"):
        assert term in prose,term
    newfiles=[note,HERE/"local_field_cosmology.py",HERE/"local_field_cosmology_results.json",
        Path(__file__),HERE/"drafts/common_background_decision.md",
        HERE/"drafts/common_recovery_v0_9.md",HERE/"drafts/publish969.py",
        STAGE/"970/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至969）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=969)==list(range(1,970))
    assert "001—969轮共969份" in nav[0].read_text("utf-8-sig")
    assert "231—969的739份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "969：局部量子交互与宇宙传播的共同来源" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"968/research_round_968_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=969,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=969,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_968=prev["cumulative_numbered_test_groups_from_967"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        same_965_preparation_interaction_and_effects=True,
        same_964_gravity_coefficients=True,
        complete_local_energy_in_mean_source=True,
        wrong_bound_free_source_dictionary_refuted=True,
        explicit_stable_support_and_monopole_inputs=True,
        real_cavity_matching=False,full_goal_completed=False,
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

