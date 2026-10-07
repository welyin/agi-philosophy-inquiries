"""Delivery checks for 970. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_970_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,970):
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
    result=read(HERE/"native_photon_scattering_results.json")
    assert result["round"]==970 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["parameters"];b=result["analytic_bounds"]
    assert p["band_min"]>0 and p["lambda_charge"]==.05
    assert math.isclose(p["G"]**2,.05**2*p["d"],rel_tol=1e-14)
    assert max(result[k] for k in ("real_space_matching_error","two_port_unitarity_error",
                                   "record_probability_identity_error"))<1e-12
    assert b["reflection_and_asymptotic_record_lower"]==100/101
    assert b["endpoint_state_error_bound"]<b["reported_endpoint_error"]==.001
    assert b["finite_record_contrast_lower"]>.986
    for row in result["packet_rows"]:
        assert abs(row["normalization"]-1)<1e-12
        assert abs(row["mean_photon_energy"]-p["Delta"])<1e-12
        assert row["record_contrast"]>.999
    r=result["resources"]
    assert r["incoming_parities_same_free_photon_energy_distribution"]
    assert not r["finite_bare_preparations_full_spectral_equality_claimed"]
    assert not r["mechanical_recoil_dynamics_included"]
    assert r["total_energy_current_balance"]<1e-14
    assert result["overlap_with_969"]["cook_duration_over_969_T"]>9000
    for key in ("spatial_dimension_derived","same_full_965_Hamiltonian",
                "Pauli_Fierz_matching_error_proved","same_969_cosmology_window_proved",
                "irreversible_arrow_proved","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core970",HERE/"native_photon_scattering.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_970.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==18
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,10)]
    for term in ("整体目标未完成","机械反作用没有完成","不是965完整",
                 "最短时间","自由光子","旋波本身是新增有效模型输入"):
        assert term in prose,term
    newfiles=[note,HERE/"native_photon_scattering.py",HERE/"native_photon_scattering_results.json",
        Path(__file__),HERE/"drafts/scattering_adoption_decision.md",
        HERE/"drafts/scattering_mechanism_increment.md",HERE/"drafts/publish970.py",
        STAGE/"971/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至970）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=970)==list(range(1,971))
    assert "001—970轮共970份" in nav[0].read_text("utf-8-sig")
    assert "231—970的740份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "970：原生电荷与传播光子的共同记录接口" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"969/research_round_969_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=970,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=970,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_969=prev["cumulative_numbered_test_groups_from_968"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        native_charge_vertex_and_old_code=True,
        positive_band_exact_single_excitation_scattering=True,
        finite_wavepacket_Cook_bound=True,
        phase_inputs_same_free_energy_distribution=True,
        common_969_history_verified=False,
        rotating_wave_parent_error_verified=False,mechanical_recoil_verified=False,
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

