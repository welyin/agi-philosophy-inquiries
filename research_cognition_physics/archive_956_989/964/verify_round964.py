"""Delivery checks for 964. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_964_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,964):
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
    result=read(HERE/"material_cosmology_results.json")
    assert result["round"]==964 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    m=result["material"];c=result["cosmology"];th=result["thermal"];scope=result["scope"]
    old958=read(STAGE/"958/capacitive_material_write_results.json")
    assert abs(m["T"]-old958["exact_interface"]["controlled_phase_time"])<1e-10
    assert m["parameters"]==dict(U=1.,v=.1,kappa=.2,rest_offset=100.)
    assert abs(m["conditional_read_contrast"]-old958["exact_interface"]["finite_time_rows"][1]["receiver_contrast"])<1e-10
    assert m["rows"][2]["receiver_holevo"]>.69 and m["rows"][3]["receiver_holevo"]<1e-5
    assert abs(th["full_material_entropy"]-math.log(2))<1e-15
    assert th["photon_entropy_production"]==0
    for row in c["rows"]:
        a=row["a"];M=m["complete_mean_mass"];R=c["photon_comoving_energy"]
        assert abs(row["total_comoving_matter_energy"]+row["accumulated_pressure_work"]-(M+R))<1e-11
        assert abs(row["hubble"]**2-(M/a**3+R/a**4)/c["mu"])<1e-22
        assert abs(row["photon_to_material_gap_ratio"]-1/a)<1e-14
        assert abs(row["photon_entropy"]-th["photon_entropy"])<1e-9
    assert c["independent_Raychaudhuri"]["max_endpoint_error"]<2e-9
    assert c["scale_at_two_T"]>2
    assert scope["same_958_material_and_readout"] and scope["uniform_mean_source_and_pressure_work_closed"]
    assert scope["semiclassical_background_is_not_full_linear_quantum_process"]
    assert not scope["gravitational_source_fluctuations_certified"]
    assert not scope["full_SM_and_GR_common_recovery"] and not scope["full_goal_completed"]
    tidal=result["conditional_tidal_ledger"]
    assert not tidal["coefficient_matched_from_material_parent"] and not tidal["perturbed_background_error_certified"]
    note=STAGE/"research_note_964.md"
    prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("同一非选择态","半经典均匀平均场模型","不是完整共同量子模型",
                 "q没有由原材料父理论匹配","整体目标未完成"):
        assert term in prose,term
    newfiles=[note,HERE/"material_cosmology.py",HERE/"material_cosmology_results.json",
        Path(__file__),HERE/"drafts/cosmology_decision.md",HERE/"drafts/common_recovery_v0_7.md",
        HERE/"drafts/update_navigation964.py",STAGE/"965/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至964）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=964)==list(range(1,965))
    assert "001—964轮共964份" in nav[0].read_text("utf-8-sig")
    assert "231—964的734份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "964：同一材料的膨胀、记录与热账" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"963/research_round_963_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=964,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=964,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_963=prev["cumulative_numbered_test_groups_from_962"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,original_material_parameters_and_readout_verified=True,
        homogeneous_mean_source_not_quantum_gravity=True,
        expansion_does_not_automatically_supply_irreversible_arrow=True,
        tidal_error_is_conditional_not_parent_matching=True,
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

