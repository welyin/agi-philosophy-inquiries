"""Publish the radiation-adoption update while preserving frozen work."""
from pathlib import Path
import re

STAGE=Path(__file__).resolve().parents[2]
RESEARCH=STAGE.parent
HERE=STAGE/'1005'


def main():
    old=(STAGE/'1004/verify_round1004.py').read_text('utf-8')
    prefix=old.split('    import importlib.util',1)[0]
    prefix=prefix.replace('finite strong-field response 1004','mechanism update 1005')
    prefix=prefix.replace('research_round_1004_checks.json','research_round_1005_checks.json')
    prefix=prefix.replace('range(776,1004)','range(776,1005)')
    checks=r'''
    import importlib.util
    spec=importlib.util.spec_from_file_location('mechanism1005',HERE/'mechanism_update_audit.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'mechanism_update_results.json')
    assert core.run()==saved
    assert saved['all_audit_checks_passed'] and saved['fresh_physical_test_groups']==0
    assert saved['cumulative_test_groups']==3786
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('radiation_transport_solution_computed',
        'same_cosmic_material_formation_history_certified',
        'overall_mechanism_map_is_complete_physical_recovery',
        'mathematical_axioms_957_changed','app_goal_changed','full_goal_completed'):
        assert not saved[key],key
    assert saved['newly_closed_cognitive_axiom_gaps']==[]
    note=STAGE/'research_note_1005.md'
    prose=note.read_text('utf-8')
    assert prose.count('$$')==4
    for term in ('整体目标未完成','累计3786','新增物理试验组0','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'overall_operation_hypothesis_v2_1.md',
        HERE/'radiation_formation_adoption_v1.md',HERE/'mechanism_update_audit.py',
        HERE/'mechanism_update_results.json',Path(__file__),
        HERE/'drafts/adoption_decision.md',HERE/'drafts/review_notes.md',
        HERE/'drafts/publish1005.py',STAGE/'1006/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1005）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    links=0
    for doc in [p for p in newfiles if p.suffix=='.md']+nav:
        content=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for pth in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',pth.name)
        if m and pth.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=1005)==list(range(1,1006))
    assert '001—1005轮共1005份' in nav[0].read_text('utf-8-sig')
    assert '231—1005的775份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '1005：辐射反馈与整体假说v2.1' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'1004/research_round_1004_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=1005,date='2026-10-08',all_delivery_checks_passed=True,
        formal_reports=1005,fresh_test_groups=0,
        cumulative_numbered_test_groups_from_1004=prev['cumulative_numbered_test_groups_from_1003'],
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,mechanism_coverage_sections=len(saved['coverage']),
        mechanism_version='v2.1',radiation_solver_run=False,
        physical_unification_certified=False,full_goal_completed=False,
        visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
    if not writing:
        before=read(TARGET)
        for key in ('frozen_inputs','new_scientific_and_entry_files'):assert before[key]==out[key],key
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    out=run(args.write)
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    print(json.dumps({k:v for k,v in out.items() if k not in
        ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
'''
    paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    originals={p:p.read_bytes() for p in paths};updates={}
    for p,raw in originals.items():
        s=raw.decode('utf-8-sig').replace('\r\n','\n')
        title='## 1005：辐射反馈与整体假说v2.1';mark='## 1004：几何、场态与有限强场记录'
        assert title not in s and s.count(mark)==1,p
        pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
        block=(title+'\n\n'
            +f'[1005报告]({pre}research_note_1005.md)采用辐射输运、再吸收及热／动量反馈；'
            +'发射、净失能和温度下降分开，同一来源与接收者共同记账。'
            +f'[整体假说v2.1]({pre}1005/overall_operation_hypothesis_v2_1.md)整合任务供给、有限强场及形成机制，'
            +'保16类现象和输入边界。正式1005／累计3786，新增物理试验0，整体未完成。\n\n'
            +f'[机制补充]({pre}1005/radiation_formation_adoption_v1.md) · '
            +f'[审计结果]({pre}1005/mechanism_update_results.json) · '
            +f'[核验]({pre}1005/research_round_1005_checks.json)。停止输运求解器细化；'
            +f'接[1006]({pre}1006/drafts/STATUS.md)补宏观量子有序相的成熟机制。'
            +'应用目标、957和条件性空间接口保持。\n\n')
        s=s.replace(mark,block+mark,1)
        if p==paths[0]:
            assert s.count('001—1004轮共1004份')==1
            s=s.replace('001—1004轮共1004份','001—1005轮共1005份',1)
        if p==paths[4]:
            before='当前正式1004／累计3786，1004已结项';assert s.count(before)==1
            s=s.replace(before,'当前正式1005／累计3786，1005已结项',1)
        if p==paths[5]:
            assert s.count('231—1004的774份')==1
            s=s.replace('# 231—1004轮阶段成果总览','# 231—1005轮阶段成果总览',1)
            s=s.replace('231—1004的774份','231—1005的775份',1)
        if p==paths[-1]:
            before='## 六条共同协议：全局缺口对应与检验优先级（截至1004）';assert s.count(before)==1
            s=s.replace(before,before.replace('1004','1005'),1)
            updates_by_id={
                'C17':('整体v2.1分存在、形成、保存、运行；1005补辐射搬运及再吸收反馈',
                    '实际材料／组成及共同形成史未核；宏观量子有序机制下一项补'),
                'C20':('整体v2.1保16类机制；辐射与材料共用交换及来源，1004保场态／仪器',
                    '机制采用非全部共同预测；强场反作用和真实形成过程仍未认证'),
                'C23':('供给按任务分类；发射、净失能、温变和记录可用性分开',
                    '实际任务准备／访问仍输入；初始箭头及唯一结果未推出'),
                'C26':('1000配对热—力；1005采用物质—辐射输运及环境反馈',
                    '谱、闭合和物态匹配仍输入；未核现实冷却量或全部宏观量子相')}
            for cid,(state,boundary) in updates_by_id.items():
                pattern=r'^\|'+cid+r' [^\n]*$';matches=re.findall(pattern,s,re.M);assert matches
                cells=matches[0].split('|');assert len(cells)==6
                cells[2]+='；1005';cells[3]=state;cells[4]=boundary
                s,count=re.subn(pattern,lambda _:'|'.join(cells),s,count=1,flags=re.M);assert count==1
            anchor='### 当前取舍\n';assert s.count(anchor)==1
            s=s.replace(anchor,anchor+'\n正式1005／累计3786。成熟辐射交换采用已补形成反馈，'
                '整体v2.1不把各自合法的过程拼成同一现实历史。下一项优先补常规超导等宏观量子有序机制，'
                '随后审中微子产生—传播—探测；不扩大候选技术门槛。\n',1)
            replacement=('4. 957数学假说v0.2保持；[1005整体v2.1](../../1005/overall_operation_hypothesis_v2_1.md)'
                '整合供给、强场和辐射形成反馈。下一项先补常规超导的有序与电磁响应机制，'
                '不计算全相图或扩建装置。')
            s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M)
            assert count==1
        if b'\r\n' in raw:s=s.replace('\n','\r\n')
        updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
    assert all(p.read_bytes()==raw for p,raw in originals.items()),'Concurrent navigation edit'
    with (HERE/'verify_round1005.py').open('x',encoding='utf-8') as dest:dest.write(prefix+checks)
    for p,data in updates.items():p.write_bytes(data)
    print('Published 1005 mechanism adoption and overall v2.1.')


if __name__=='__main__':main()
