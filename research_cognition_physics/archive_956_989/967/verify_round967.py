"""Delivery checks for 967. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_967_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,967):
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
    result=read(HERE/"collective_material_results.json")
    assert result["round"]==967 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    assert result["parameters"]==dict(U=1.,v=.1,absolute_kappa=.2,
        active_dimension=216,full_with_old_spectator_spins=13824)
    from fractions import Fraction as F
    assert F(result["loop_sign_energy_difference_interval"]["lo"])>0
    bound=2*F(9,250)*F(347,1000)/F(219,500)
    assert result["fixed_menu_rows"][0]["connected_energy"]>0
    assert result["fixed_menu_rows"][1]["connected_energy"]<0
    for row in result["fixed_menu_rows"]:
        c=row["ground_interval"];z=row["connected_energy_interval"]
        assert c["lower_inertia"]==[0,27] and c["upper_inertia"]==[1,26]
        assert F(c["hi"])-F(c["lo"])==F(5,10**12)
        assert F(z["lo"])*F(z["hi"])>0
        assert F(row["all_time_uniform_input_error_bound"])==bound
        assert row["direct_uniform_error"]<float(bound)
        assert max(row["intertwining_error"],row["isometry_error"])<2e-13
        assert row["pair_only_prediction_probability"]==1.
        assert row["actual_receiver_probability"]<2e-5
        assert row["certified_prediction_discrepancy_lower"]>.885
        assert row["full_trace"]==216.
        assert row["energy_current_balance"]<1e-12
    assert not result["scope"]["independent_three_body_coupling_introduced"]
    assert not result["scope"]["autonomous_geometry_or_binding_proved"]
    # Recompute actual CAR-to-singlet mapping, rational LDL signs, all 64
    # dressed columns, the physical receiving effect and shared sources.
    import importlib.util
    spec=importlib.util.spec_from_file_location("core967",HERE/"collective_material.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_967.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","第三份材料","精确双体项","接触系数仍为外部固定输入",
                 "同一H的不同记录部门","不反证量子局部层析"):
        assert term in prose,term
    newfiles=[note,HERE/"collective_material.py",HERE/"collective_material_results.json",
        Path(__file__),HERE/"drafts/collective_material_decision.md",
        HERE/"drafts/collective_mechanism_increment.md",HERE/"drafts/publish967.py",
        STAGE/"968/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至967）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=967)==list(range(1,968))
    assert "001—967轮共967份" in nav[0].read_text("utf-8-sig")
    assert "231—967的737份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "967：原生材料的集体作用与共同来源" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"966/research_round_966_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=967,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=967,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_966=prev["cumulative_numbered_test_groups_from_965"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        exact_rational_ground_intervals=True,
        same_native_material_and_density_rule=True,
        no_independent_three_body_parameter=True,
        all_low_inputs_and_passive_references_bound=True,
        autonomous_geometry_certified=False,
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

