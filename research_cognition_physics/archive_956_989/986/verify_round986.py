"""Delivery checks for 986. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_986_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,986):
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
    extra["981/drafts/common_adoption_audit_checks.json"]=("frozen_evidence_hashes","audit_file_hashes")
    extra["985/drafts/common_model_adoption_audit_checks.json"]=("frozen_evidence_hashes","audit_file_hashes")
    for path,keys in extra.items():
        value=read(STAGE/path)
        for key in keys:add(value[key])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/"metric_dipole_bridge_results.json")
    assert result["round"]==986 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    import importlib.util
    import numpy as np
    from fractions import Fraction
    spec=importlib.util.spec_from_file_location("core986",HERE/"metric_dipole_bridge.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    bound=result["infinite_domain_bound"]
    assert Fraction(bound["rational_Duhamel_upper"]["exact"])<Fraction(134,10**9)
    assert Fraction(bound["record_contrast_lower"]["exact"])>Fraction(4859478,10**7)
    assert bound["all_times_in_original_interval"] and bound["unknown_input_and_passive_reference"]
    # Independent second source moment check: retain cross terms explicitly.
    material=mod.load("native_check986",STAGE/"968/internal_relay.py")
    m,h,q,_,w,_=material.material();n=5
    a=np.diag(np.sqrt(np.arange(1,n)),1);I=np.eye(4)
    S=mod.kron(q,I)+mod.kron(I,q)
    Cbare=.5*mod.kron(np.eye(16),a@a+a.T@a.T)
    mixed=.003*mod.kron(S,a+a.T);C=Cbare+mixed
    photon=np.zeros(n);photon[:2]=1/math.sqrt(2)
    plus=w@np.ones(2)/math.sqrt(2)
    W=np.column_stack([np.kron(np.kron(w[:,i],plus),photon) for i in (0,1)])
    difference=(C@W).T@(C@W)-(Cbare@W).T@(Cbare@W)
    d=(1-1/math.sqrt(1.16))/2
    expected=.003**2*d*np.diag([3.,1.])
    assert np.linalg.norm(difference-expected)<1e-14
    assert np.linalg.norm(W.T@(Cbare@mixed+mixed@Cbare)@W)<1e-14
    for key in ("independent_quadrupole_matching_derived","full_U1_completed",
        "spatial_mode_matching_error_certified","full_GR_error_certified",
        "Newton_motion_985_and_TT_glued","all_SM_and_cosmology_completed",
        "physical_minimum_scale_assumed","full_goal_completed"):
        assert result["scope"][key] is False
    note=STAGE/"research_note_986.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","停止本模式、耦合和精度优化",
        "不声称该延拓恢复强场GR","独立固有四极在所选TT偏振上的投影为零",
        "不等于所有无界高阶物理余项都有同样误差界"):
        assert term in prose,term
    newfiles=[note,HERE/"metric_dipole_bridge.py",HERE/"metric_dipole_bridge_results.json",
        Path(__file__),HERE/"drafts/metric_dipole_decision.md",HERE/"drafts/metric_dipole_adoption.md",
        HERE/"drafts/publish986.py",STAGE/"987/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至986）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=986)==list(range(1,987))
    assert "001—986轮共986份" in nav[0].read_text("utf-8-sig")
    assert "231—986的756份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "986：原电磁交互与传播几何的共同生成元" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"985/research_round_985_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=986,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=986,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_985=prev["cumulative_numbered_test_groups_from_984"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,native_EM_metric_source_joined=True,
        independent_mixed_source_second_moment_verified=True,
        full_Fock_finite_time_record_bound_verified=True,
        entire_U1_completed=False,full_physical_error_certified=False,
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

