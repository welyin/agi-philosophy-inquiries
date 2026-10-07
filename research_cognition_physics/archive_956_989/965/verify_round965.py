"""Delivery checks for 965. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_965_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,965):
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
    result=read(HERE/"material_field_window_results.json")
    physical=read(HERE/"physical_mode_effect_results.json")
    assert result["round"]==965 and result["all_scientific_checks_passed"]
    assert physical["round"]==965 and physical["all_checks_passed"]
    for data in (result,physical):
        for rel,digest in data["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    previous=read(STAGE/"958/capacitive_material_write_results.json")
    assert result["T"]==previous["exact_interface"]["controlled_phase_time"]
    assert result["old_material_parameters"]==dict(U=1.,v=.1,kappa=.2)
    assert [r["cutoff"] for r in result["fixed_menu_rows"]]==[7,10]
    for row in result["fixed_menu_rows"]:
        eps=row["residual_bound"]["total"]
        assert eps<2e-6 and row["residual_bound"]["defining_matrix_operator_error"]<1e-12
        assert row["probabilities"]["receiver_contrast"]-4e-6>.998723
        assert row["algebra"]["low_charge_projection"]<1e-14
        assert row["algebra"]["projected_square_norm"]>.07
        assert row["algebra"]["charge_internal_commutator"]>.39
        for energy in row["sources"]["energy_rows"]:
            assert abs(sum(energy["after"])-sum(energy["before"]))<1e-11
        assert row["sources"]["current_identity_residual"]<1e-12
        assert not row["scope"]["local_propagation_or_full_Maxwell_recovered"]
        assert not row["scope"]["unbounded_energy_observable_error_transported"]
    for row in physical["rows"]:
        assert row["conservative_reported_lower"]>.049871
        assert row["maximum_compressed_effect_norm"]<=1+1e-12
    # Independently check that the near-unitary normalization defect is covered
    # by the explicit endpoint guard, on the WHOLE declared input isometry.
    import importlib.util
    import numpy as np
    spec=importlib.util.spec_from_file_location("core965",HERE/"material_field_window.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    m,Q,Hm,S,W0,local=mod.material()
    n=11;g=.003;omega=.5
    a=np.diag(np.sqrt(np.arange(1,n)),1)
    H=np.kron(Hm,np.eye(n))+np.kron(np.eye(36),omega*np.diag(np.arange(n)))
    H+=g*np.kron(S,a+a.T)+g*g/omega*np.kron(S@S,np.eye(n))
    ev,V=np.linalg.eigh(H);Vin=np.kron(W0,np.eye(n)[:,:2])
    output=V@(mod.phases(ev,result["T"])[:,None]*(V.T@Vin))
    gram_defect=float(np.linalg.norm(output.conj().T@output-np.eye(32)))
    assert gram_defect<1e-10
    note=STAGE/"research_note_965.md"
    prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("单模读口不等于空间传播","联合算符","整体目标未完成",
                 "没有用有界读数的全态误差自动认证无界"):
        assert term in prose,term
    newfiles=[note,HERE/"material_field_window.py",HERE/"material_field_window_results.json",
        HERE/"physical_mode_effect.py",HERE/"physical_mode_effect_results.json",
        Path(__file__),HERE/"drafts/electromagnetic_decision.md",
        HERE/"drafts/common_recovery_v0_8.md",HERE/"drafts/publish965.py",
        STAGE/"966/drafts/STATUS.md"]
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至965）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=965)==list(range(1,966))
    assert "001—965轮共965份" in nav[0].read_text("utf-8-sig")
    assert "231—965的735份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "965：原生材料的共同电磁读口" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"964/research_round_964_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=965,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=965,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_964=prev["cumulative_numbered_test_groups_from_963"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,uniform_input_gram_defect=gram_defect,
        old_material_and_receiver_reused=True,
        infinite_occupation_residual_bound_used=True,
        field_observable_includes_polarization=True,
        single_mode_not_local_propagation=True,
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

