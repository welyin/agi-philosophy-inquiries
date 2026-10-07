"""Publish mechanism synthesis v2 without upgrading evidence to physical proof."""
from pathlib import Path
import re

STAGE=Path(__file__).resolve().parents[2]
RESEARCH=STAGE.parent
HERE=STAGE/'1003'


def main():
    old=(STAGE/'1002/verify_round1002.py').read_text('utf-8')
    prefix=old.split('    import importlib.util',1)[0]
    prefix=prefix.replace('thermal record reuse interface 1002','mechanism synthesis 1003')
    prefix=prefix.replace('research_round_1002_checks.json','research_round_1003_checks.json')
    prefix=prefix.replace('range(776,1002)','range(776,1003)')
    checks=r'''
    import importlib.util
    spec=importlib.util.spec_from_file_location('synthesis1003',HERE/'mechanism_synthesis_audit.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'mechanism_synthesis_results.json');assert core.run()==saved
    assert saved['all_audit_checks_passed'] and saved['round']==1003
    assert saved['fresh_physical_test_groups']==0 and saved['cumulative_test_groups']==3785
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('new_physics_theorem','joint_lifecycle_certified','all_physical_phenomena_explained',
        'full_goal_completed','coverage_is_physical_proof','mathematical_axioms_957_changed',
        'application_goal_changed','old_numerical_models_rerun'):assert not saved[key],key
    assert len(saved['coverage'])==16
    assert {c for row in saved['coverage'] for c in row['conditions']}=={f'C{i:02d}' for i in range(1,28)}
    assert sum(j['state_connection_verified'] for j in saved['connections'])==1
    join=saved['connections'][0]
    assert join['requires_new_resources'] and not join['arbitrary_phase_lifecycle_verified']
    note=STAGE/'research_note_1003.md';prose=note.read_text('utf-8')
    for term in ('整体目标未完成','累计3785','新增物理试验0','386或425','同一初态'):
        assert term in prose,term
    newfiles=[note,HERE/'overall_operation_hypothesis_v2.md',HERE/'mechanism_synthesis_audit.py',
        HERE/'mechanism_synthesis_results.json',Path(__file__),HERE/'drafts/adoption_decision.md',
        HERE/'drafts/review_notes.md',HERE/'drafts/publish1003.py',STAGE/'1004/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1003）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=1003)==list(range(1,1004))
    assert '001—1003轮共1003份' in nav[0].read_text('utf-8-sig')
    assert '231—1003的773份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '1003：整体假说v2与共同解释审计' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'1002/research_round_1002_checks.json',HERE/'drafts/STATUS.md',
              STAGE/'999/mechanism_adoption_results.json']
    prev=read(oldfiles[0])
    out=dict(round=1003,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=1003,fresh_test_groups=0,
        cumulative_numbered_test_groups_from_1002=prev['cumulative_numbered_test_groups_from_1001'],
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,phenomenon_groups=16,obligations_per_group=4,conditions_retained=27,
        mechanism_version='v2',coverage_is_physical_proof=False,joint_lifecycle_certified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
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
        title='## 1003：整体假说v2与共同解释审计';mark='## 1002：热材料再写与功能分工'
        assert title not in s and s.count(mark)==1,p
        pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
        block=(title+'\n\n'
            +f'[1003报告]({pre}research_note_1003.md)发布[整体假说v2]({pre}1003/overall_operation_hypothesis_v2.md)，'
            +'将核组成、材料相和有限记录接入16类现象说明；区分四项组织许可、任务资源与共同历史。'
            +f'[审计结果]({pre}1003/mechanism_synthesis_results.json) · '
            +f'[核验]({pre}1003/research_round_1003_checks.json)。正式1003／累计3785；新增物理试验0，整体未完成。\n\n'
            +'近期仅980→1002有新增资源下的实际热态接续，未签收相位或全功能生命周期。'
            +f'接[1004]({pre}1004/drafts/STATUS.md)补有限非平衡供给到任务准备；复用旧资源与参考结论，'
            +'停止写入器优化。应用目标、957及条件性空间接口保持。\n\n')
        s=s.replace(mark,block+mark,1)
        if p==paths[0]:
            assert s.count('001—1002轮共1002份')==1
            s=s.replace('001—1002轮共1002份','001—1003轮共1003份',1)
        if p==paths[4]:
            before='当前正式1002／累计3785，1002已结项';assert s.count(before)==1
            s=s.replace(before,'当前正式1003／累计3785，1003已结项',1)
        if p==paths[5]:
            assert s.count('231—1002的772份')==1
            s=s.replace('# 231—1002轮阶段成果总览','# 231—1003轮阶段成果总览',1)
            s=s.replace('231—1002的772份','231—1003的773份',1)
        if p==paths[-1]:
            before='## 六条共同协议：全局缺口对应与检验优先级（截至1002）';assert s.count(before)==1
            s=s.replace(before,before.replace('1002','1003'),1)
            rows={
                'C17':('|C17 质量、耦合、稳定物质结构|2、3、6；434、941—958、971、981、998—1003|'
                    '整体v2将存在、形成、保存、运行分开；核组成、选定相和热记录机制共同登记|'
                    '机制采用未变成同一材料形成史，现实相／参数与操作准备仍输入|'),
                'C20':('|C20 跨尺度、误差与共同来源|4、6、H3；854、898、899、936、950、957、971、981—1003|'
                    '16类现象逐项保来源与证据层；近期980→1002仅人口合同有实际后态接续|'
                    '共同作用和说明覆盖不等于共同实现；不同实验允许不同准备|'),
                'C23':('|C23 统计、退相干与时间箭头|1、4、6；305、354、522—523、929、947、961—1003|'
                    '功能按变量、时段和精度分工；供能者的熵、资料与参考资源不得遗漏|'
                    '非平衡资源到实际准备尚需机制接通；唯一结果及初始箭头未推得|')}
            for cid,row in rows.items():
                s,count=re.subn(r'^\|'+cid+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M)
                assert count==1,cid
            anchor='### 当前取舍\n';assert s.count(anchor)==1
            s=s.replace(anchor,anchor+'\n正式1003／累计3785。整体v2采用近期物质、热与记录机制，'
                '保留16类现象及27条件；文档审计不计物理试验、不签收共同历史。'
                '下一项补非平衡供给到任务准备，整体仍未完成。以下旧取舍按当时范围保留。\n',1)
            replacement=('4. 957数学假说v0.2保持；[1003整体v2](../../1003/overall_operation_hypothesis_v2.md)'
                '区分机制采用、条件实例、共同实现与现实恢复。下一项补任务准备的资源分类及供给，'
                '不继续修理1002写入器。')
            s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M)
            assert count==1
        if b'\r\n' in raw:s=s.replace('\n','\r\n')
        updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
    assert all(p.read_bytes()==raw for p,raw in originals.items()),'Concurrent navigation edit'
    with (HERE/'verify_round1003.py').open('x',encoding='utf-8') as dest:dest.write(prefix+checks)
    for p,data in updates.items():p.write_bytes(data)
    print('Published 1003 mechanism synthesis v2.')


if __name__=='__main__':main()
