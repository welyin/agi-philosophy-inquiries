"""Delivery checks for 966. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_966_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,966):
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
    result=read(HERE/"comparison_loop_results.json")
    assert result["round"]==966 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    assert result["physical_sector_dimension"]==64
    assert result["low_organization_dimension"]==16
    assert sum(r["survives"] for r in result["exact_third_order_surviving_words"])==6
    assert result["tree_third_order_zero"]
    assert result["fixed_joint_dictionary_factorizes_H_exactly"]
    assert result["same_dictionary_makes_original_loop_effect_record_dependent"]
    assert not result["content_dependent_spectrum_generated"]
    from fractions import Fraction
    for row in result["fixed_menu_rows"]:
        e=Fraction(str(row["epsilon"]))
        eta=Fraction(63,8)*e+Fraction(117,16)*e*e+Fraction(3,2)*e**3
        assert eta==Fraction(row["all_input_uniform_error_bound"])
        assert Fraction(row["certified_loop_contrast_lower"])==Fraction(4,5)-4*eta
        assert Fraction(row["certified_mean_replacement_discrepancy_lower"])==Fraction(4,5)-2*eta
        assert row["actual_loop_contrast"]>=float(Fraction(4,5)-4*eta)>0
        assert row["direct_uniform_operator_error"]<float(eta)
        assert row["mixed_record_negative_flux_probability"]>float(Fraction(4,5)-2*eta)
        for entry in row["energy_parts"]:
            assert abs(sum(entry["changes"]))<1e-12
        assert row["full_energy_current_residual"]<1e-12
    assert result["fixed_menu_rows"][1]["certified_loop_contrast_lower_float"]>=.641768
    assert not result["scope"]["SM_or_GR_derived"]
    assert not result["scope"]["actual_spacetime_propagation_proved"]
    # Re-execute the complete physical-sector reconstruction and algebra.
    import importlib.util
    spec=importlib.util.spec_from_file_location("core966",HERE/"comparison_loop.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.compute(),result)
    note=STAGE/"research_note_966.md"
    prose=note.read_text("utf-8")
    assert prose.count("$$")==18
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,10)]
    for term in ("整体目标未完成","降为说明性例","固定联合变换","额外建模输入",
                 "不是对FUCP","不能保持原准备和全部效果不变而删除记录因子"):
        assert term in prose,term
    newfiles=[note,HERE/"comparison_loop.py",HERE/"comparison_loop_results.json",
        Path(__file__),HERE/"drafts/mechanism_priority_entry.md",
        HERE/"drafts/mechanism_map_v0_2.md",HERE/"drafts/publish966.py",
        STAGE/"967/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至966）",1)[1].split("### 当前取舍",1)[0]
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
    assert sorted(n for n in nums if n<=966)==list(range(1,967))
    assert "001—966轮共966份" in nav[0].read_text("utf-8-sig")
    assert "231—966的736份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "966：局部比较的闭路生成与表示筛选" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"965/research_round_965_checks.json",HERE/"drafts/STATUS.md"]
    prev=read(oldfiles[0])
    out=dict(round=966,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=966,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_965=prev["cumulative_numbered_test_groups_from_964"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        exact_comparison_induced_record_loop=True,
        finite_time_all_input_bound=True,
        exact_joint_dictionary_factorization=True,
        candidate_deprioritized_after_representation_screen=True,
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

