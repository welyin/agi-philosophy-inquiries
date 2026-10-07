"""Publish an adoption audit, without adding a scientific round or experiment."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
STAGE=HERE.parents[1];RESEARCH=STAGE.parent
previous=(STAGE/'984/verify_round984.py').read_text('utf-8')
prefix=previous.split('    result=read(HERE/"smeared_source_geometry_results.json")',1)[0]
prefix=prefix.replace('Delivery checks for 984','Adoption audit after 984')
prefix=prefix.replace('STAGE=HERE.parent','STAGE=HERE.parents[1]')
prefix=prefix.replace('ROOT=HERE.parents[2]','ROOT=HERE.parents[3]')
prefix=prefix.replace('TARGET=HERE/"research_round_984_checks.json"',
                      'TARGET=HERE/"common_model_adoption_audit_checks.json"')
prefix=prefix.replace('range(776,984)','range(776,985)')
tail=r'''    last=read(STAGE/"984/research_round_984_checks.json")
    assert last["all_delivery_checks_passed"] and last["formal_reports"]==984
    assert last["cumulative_numbered_test_groups_from_983"]==3769
    assert not (STAGE/"research_note_985.md").exists()
    r=read(STAGE/"984/smeared_source_geometry_results.json")
    assert r["scope"]["matter_induced_linear_Einstein_report_only"]
    assert not r["scope"]["finite_spacetime_detector_certified"]
    assert not r["scope"]["full_quantum_geometry_error_certified"]
    paths=[STAGE/f"research_note_{n}.md" for n in
           (898,936,939,940,950,951,953,957,961,964,969,971,973,978,979,980,981,982,983,984)]
    paths += [STAGE/"957/drafts/unified_operation_hypotheses_v0_2.md",
        STAGE/"951/drafts/joint_minimum_menu.md",STAGE/"981/drafts/common_parent_contract_v1.md",
        STAGE/"983/drafts/parent_recovery_scope_v2.md",STAGE/"984/drafts/geometry_adoption_scope.md",
        STAGE/"962/drafts/revised_goal_20261007.txt",STAGE/"962/drafts/app_goal_confirmation.json",
        STAGE/"984/research_round_984_checks.json",HERE/"STATUS.md"]
    for p in paths:assert p.exists(),p
    note=HERE/"research_note_985_working.md"
    s=note.read_text("utf-8")
    for term in ("不是新科学轮次","正式984／累计3769保持","U1：同一实际物理生成过程",
        "U2：实际来源与几何报告","U3：共同字典及采用域","不再接续984的噪声核细化",
        "不是新认知公理"):
        assert term in s,term
    assert "作用的每个量子／经典投影已匹配" in s
    audit=[note,HERE/"publish_adoption_audit.py",Path(__file__)]
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    links=0
    for p in audit:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    for p in [note]+nav:
        text=p.read_text("utf-8-sig")
        if p in nav:assert "984后整体采用审计：先接共同物理过程" in text
        content=re.sub(r"\$\$.*?\$\$","",text,flags=re.S)
        for link in re.findall(r"\]\(([^)]+)\)",content):
            if re.match(r"^[a-zA-Z]+://",link) or link.startswith("#"):continue
            target=(p.parent/link.split("#")[0].strip("<>")).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(p,link)
            links+=1
    assert "001—984轮共984份" in nav[0].read_text("utf-8-sig")
    assert "231—984的754份" in nav[5].read_text("utf-8-sig")
    ledger=nav[-1].read_text("utf-8-sig")
    section=ledger.split("## 六条共同协议：全局缺口对应与检验优先级（截至984）",1)[1].split("### 当前取舍",1)[0]
    assert re.findall(r"^\|(C\d\d) ",section,re.M)==[f"C{i:02d}" for i in range(1,28)]
    out=dict(date="2026-10-07",all_audit_checks_passed=True,formal_reports=984,
        cumulative_numbered_test_groups=3769,new_scientific_test_groups=0,
        previous_goal_turn="progress",current_audit="adoption_and_priority_decision",
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,full_goal_completed=False,app_goal_changed=False,
        next_focus="U1: P981 physical process and common task dictionary, before U2 classical geometry",
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
target=HERE/'verify_adoption_audit.py';assert not target.exists()
target.write_text(prefix+tail,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 984后整体采用审计：先接共同物理过程'
    previous='## 984：共同量子来源与有限尺度的经典几何报告'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
       +f'[985工作报告]({pre}985/drafts/research_note_985_working.md)核对957—984与旧恢复工具：'
       +'物理作用、系数和准备允许明示输入，完整微观制造与UV完成退出默认队列；'
       +'实际生成元、来源及交叠预测仍须同一字典。先接P981共同物理过程，再选择经典几何恢复域。'
       +f'[审计回执]({pre}985/drafts/common_model_adoption_audit_checks.json)。'
       +'本次无新增科学实验，正式984／累计3769保持，985尚未科学结项。目标与957假说保持，整体未完成。\n\n')
    s=s.replace(previous,block+previous,1)
    if p==paths[-1]:
        rows={
        'C18':'|C18 共同作用、来源与反作用|1、4、5、H2；936、940、971、981—984；985采用审计|共同父合同及限定来源接口已有；当前先接实际物理生成过程U1，再取几何报告U2|不能把实际正过程与形式物理系数分属两份对象后直接合并；不要求全UV模型|',
        'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；854、898、899、936、950、957、971、981—984；985审计|有效系数可作为明示物理输入；同任务、同准备、同源字典及重复计数须核|全微观材料制造不作默认门槛；模型内界不自动等于父物理保留阶或全物理误差|',
        'C21':'|C21 完整物质与相互作用|1、4、5；531—553、735、939、951、957、971、981、983—984；985审计|P981三代SM及可指定有效系数保留；自由共形阶的4／45／12已有具体核|完整物理字段表不等于已完成全部共同过程；不以微观唯一计算每个材料参数作为前置关卡|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[985工作审计](../../985/drafts/research_note_985_working.md)把剩余工作聚焦U1共同物理过程、U2实际来源与几何报告、U3共同字典。停止近期局部例子的默认续修，先接940／971到P981的保留阶对象。正式984／累计3769保持，审计不算新科学轮。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published the adoption audit; scientific round remains 984.')
