"""Publish an adoption audit, without incrementing formal research rounds."""
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
STAGE=HERE.parents[1]
RESEARCH=STAGE.parent
old=(STAGE/"986/verify_round986.py").read_text("utf-8")
prefix=old.split('    result=read(HERE/"metric_dipole_bridge_results.json")',1)[0]
for a,b in (
    ('Delivery checks for 986. Passing is not proof of complete physical unification.',
     'Overlap audit after 986. Passing is not proof of complete physical unification.'),
    ('STAGE=HERE.parent\n','STAGE=HERE.parents[1]\n'),
    ('ROOT=HERE.parents[2]\n','ROOT=HERE.parents[3]\n'),
    ('research_round_986_checks.json','common_overlap_audit_checks.json'),
    ('range(776,986)','range(776,987)')):
    assert a in prefix,a
    prefix=prefix.replace(a,b)
checks=r'''    import importlib.util
    from fractions import Fraction
    spec=importlib.util.spec_from_file_location("audit987",HERE.parent/"overlap_adoption.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    result=read(HERE.parent/"overlap_adoption_results.json")
    mod.compare(mod.run(),result)
    for path,digest in result["source_hashes"].items():assert sha(ROOT/path)==digest,path
    assert Fraction(result["common_record"]["state_trace_distance_upper"]["exact"])==Fraction(860567,500000000)
    assert result["conditional_profile"]["k_sigma"]==20
    for row in result["conditional_profile"]["linear_source_second_moment_by_lapse"]:
        for constant,profiled,difference in zip(
            row["linear_vertex_force_second_moment_constant"],
            row["linear_vertex_force_second_moment_profiled"],row["absolute_difference"]):
            assert constant>1e-29 and 6e-30<profiled<7e-30
            assert math.isclose(profiled/constant,.5,rel_tol=1e-13,abs_tol=0)
            assert math.isclose(difference,constant-profiled,rel_tol=1e-13,abs_tol=0)
    assert not result["spectator_comparison"]["position_kinetic_energy_included"]
    assert not result["spectator_comparison"]["proves_join_with_985_motion"]
    assert not result["decision"]["actual_finite_record_failure_demonstrated"]
    assert not result["decision"]["full_goal_completed"]
    assert min(result["established_EFT"]["full_binary_k_r"])>9e7
    last=read(STAGE/"986/research_round_986_checks.json")
    assert last["all_delivery_checks_passed"] and last["formal_reports"]==986
    assert last["cumulative_numbered_test_groups_from_985"]==3771
    assert not (STAGE/"research_note_987.md").exists()
    note=HERE/"research_note_987_working.md"
    prose=note.read_text("utf-8")
    for text in ("不新增正式科学轮次","正式986／累计3771保持",
        "不是完整非线性力的全时间二阶矩","不增加微观连续或物理最小尺度假设",
        "不是985运动已被证明可忽略于完整合并模型","不是父理论的总物理误差"):
        assert text in prose,text
    assert prose.count("$$")==10
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,6)]
    audit=[note,HERE/"overlap_decision.md",HERE/"publish_overlap_audit.py",Path(__file__),
        HERE.parent/"overlap_adoption.py",HERE.parent/"overlap_adoption_results.json"]
    paths=[STAGE/"986/research_round_986_checks.json",HERE/"STATUS.md",
        STAGE/"985/drafts/research_note_985_working.md",
        STAGE/"962/drafts/revised_goal_20261007.txt",
        STAGE/"962/drafts/app_goal_confirmation.json"]
    for path in result["source_hashes"]:
        p=ROOT/path
        if p not in audit and p not in paths:paths.append(p)
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    links=0
    for p in audit:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    for p in [a for a in audit if a.suffix==".md"]+nav:
        text=p.read_text("utf-8-sig")
        if p in nav:assert "986后合并审计：记录交叠与来源边界" in text
        content=re.sub(r"\$\$.*?\$\$","",text,flags=re.S)
        for link in re.findall(r"\]\(([^)]+)\)",content):
            if re.match(r"^[a-zA-Z]+://",link) or link.startswith("#"):continue
            target=(p.parent/link.split("#")[0].strip("<>")).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(p,link)
            links+=1
    assert "001—986轮共986份" in nav[0].read_text("utf-8-sig")
    assert "231—986的756份" in nav[5].read_text("utf-8-sig")
    ledger=nav[-1].read_text("utf-8-sig")
    section=ledger.split("## 六条共同协议：全局缺口对应与检验优先级（截至986）",1)[1].split("### 当前取舍",1)[0]
    assert re.findall(r"^\|(C\d\d) ",section,re.M)==[f"C{i:02d}" for i in range(1,28)]
    out=dict(date="2026-10-07",all_audit_checks_passed=True,formal_reports=986,
        cumulative_numbered_test_groups=3771,new_scientific_test_groups=0,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,full_goal_completed=False,app_goal_changed=False,
        full_985_986_joint_model_claimed=False,
        next_focus="P981 common retained-order source and task dictionary; no mode optimization",
        frozen_evidence_hashes={str(p.relative_to(ROOT)):sha(p) for p in paths},
        audit_file_hashes={str(p.relative_to(ROOT)):sha(p) for p in audit})
    if not writing:
        saved=read(TARGET)
        for key in ("frozen_evidence_hashes","audit_file_hashes"):assert out[key]==saved[key],key
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
        ("frozen_evidence_hashes","audit_file_hashes")},ensure_ascii=False,indent=2))
'''
# The surrounding raw string should leave a single regex escape for a literal backslash.
checks=checks.replace('r"\\\\tag\\{(\\d+)\\}"','r"\\\\tag\{(\d+)\}"')
target=HERE/"verify_overlap_audit.py"
assert not target.exists()
target.write_text(prefix+checks,encoding="utf-8")

paths=[RESEARCH/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
    STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                     "_shared/notes/unified_physics_condition_ledger_current.md")]
original={p:p.read_bytes() for p in paths}
updates={}
for p,raw in original.items():
    s=raw.decode("utf-8-sig").replace("\r\n","\n")
    pre="archive_764_/" if p.parent==RESEARCH else "../../" if p==paths[-1] else ""
    heading="## 986后合并审计：记录交叠与来源边界"
    previous="## 986：原电磁交互与传播几何的共同生成元"
    assert heading not in s and s.count(previous)==1
    block=(heading+"\n\n"
        +f'[987工作报告]({pre}987/drafts/research_note_987_working.md)将985／986运输到同一记录读口，'
        +'输出迹距离上界0.001721134；空间来源不能随之直接移植。实际位置宽度与一个相容驻波给kσ=20，'
        +'来源二阶矩约差50%，其绝对影响及任务范围单独列账。'
        +f'[结果]({pre}987/overlap_adoption_results.json) · '
        +f'[审计核验]({pre}987/drafts/common_overlap_audit_checks.json)。\n\n'
        +'本次复用旧证明作合并审计，不增正式轮次；正式986／累计3771保持，整体未完成。'
        +'停止局部模式／参数优化，下一项先用成熟有效理论核P981共同保留阶的来源及任务字典。'
        +'应用目标已核为现行版本，保持active；历史和957假说不改。\n\n')
    s=s.replace(previous,block+previous,1)
    if p==paths[-1]:
        replacement=('4. 957假说v0.2保持；[987工作审计](../../987/drafts/research_note_987_working.md)'
            +'签收985／986共同记录的交叠界，不把该界当作空间来源匹配或完整联合模型。'
            +'下一项先用成熟物理结果核P981共同保留阶的来源及任务字典；停止当前模式、位置宽度和波长优化。')
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M)
        assert count==1
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Published eight live navigation entries. Formal round count remains 986.")
