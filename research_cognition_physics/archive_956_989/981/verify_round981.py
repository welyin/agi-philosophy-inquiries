"""Delivery checks for 981. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_981_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,981):
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
    for path,keys in extra.items():
        value=read(STAGE/path)
        for key in keys:add(value[key])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/"native_response_source_results.json")
    assert result["round"]==981 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["parameters"];r=result["response"];pulse=result["finite_pulse"]
    assert p["native_dimension"]==6 and p["U"]==1 and p["v"]==.1 and p["eta"]==0
    assert r["spectral_error"]<1e-13 and r["lapse_five_point_error"]<1e-10
    assert abs(r["frozen_proper_frequency_missing_fraction_at_point6"]-9/17)<1e-13
    assert r["derivative_expansion"][1]["alpha_relative_remainder"]==.0001
    assert r["derivative_expansion"][1]["lapse_relative_remainder"]<.0005
    assert pulse["certified_contrast_lower"]>5.1e-5 and pulse["lapse_slope_difference"]<1e-12
    assert result["static_source"]["hellmann_feynman_error"]<1e-11
    for key in ("full_SM_matching_completed","Einstein_background_solved",
        "dynamic_quantum_geometry_constructed","all_unknown_input_pulse_contrast",
        "finite_pulse_is_autonomous_instrument","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core981",HERE/"native_response_source.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    spec=importlib.util.spec_from_file_location("cert981",HERE/"certify_response_bounds.py")
    cert=importlib.util.module_from_spec(spec);spec.loader.exec_module(cert)
    assert cert.run()==read(HERE/"response_bounds_certificate.json")
    note=STAGE/"research_note_981.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==14
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,8)]
    for term in ("整体目标未完成","不使用上述低频导数截断",
        "不是总Maxwell—材料能源或完整T00","不是新发现Kubo公式",
        "停止本混合响应的器件和更高频优化"):
        assert term in prose,term
    newfiles=[note,HERE/"native_response_source.py",HERE/"native_response_source_results.json",
        HERE/"certify_response_bounds.py",HERE/"response_bounds_certificate.json",Path(__file__),
        HERE/"drafts/parent_matching_decision.md",HERE/"drafts/common_parent_contract_v1.md",
        HERE/"drafts/publish981.py",STAGE/"982/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至981）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=981)==list(range(1,982))
    assert "001—981轮共981份" in nav[0].read_text("utf-8-sig")
    assert "231—981的751份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "981：共同父描述与电磁—时间来源匹配" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"980/research_round_980_checks.json",HERE/"drafts/STATUS.md",
        HERE/"drafts/common_adoption_audit_checks.json"]
    prev=read(oldfiles[0])
    out=dict(round=981,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=981,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_980=prev["cumulative_numbered_test_groups_from_979"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        common_parent_contract_specified=True,
        native_causal_response_and_lapse_insertion_verified=True,
        low_frequency_relative_remainders_verified=True,
        rational_finite_pulse_contrast_verified=True,
        dynamic_quantum_geometry_completed=False,
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

