"""Publish a working adoption audit without incrementing formal research counts."""
from pathlib import Path
import re

STAGE=Path(__file__).resolve().parents[2]
RESEARCH=STAGE.parent
HERE=STAGE/'1004'


def main():
    old=(STAGE/'1003/verify_round1003.py').read_text('utf-8')
    prefix=old.split('    import importlib.util',1)[0]
    prefix=prefix.replace('Verify mechanism synthesis 1003','Verify working resource adoption')
    prefix=prefix.replace("TARGET=HERE/'research_round_1003_checks.json'", "TARGET=HERE/'drafts/resource_adoption_checks.json'")
    prefix=prefix.replace('range(776,1003)','range(776,1004)')
    checks=r'''
    import importlib.util
    spec=importlib.util.spec_from_file_location('resourceaudit',HERE/'resource_preparation_audit.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    result=read(HERE/'resource_preparation_audit_results.json');assert core.run()==result
    assert result['all_audit_checks_passed'] and result['new_scientific_test_groups']==0
    assert result['formal_reports']==1003 and result['cumulative_numbered_test_groups']==3785
    for rel,digest in result['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert len(result['resource_classification'])==3
    assert not any(row['natural_supply_certified'] for row in result['resource_classification'])
    assert not result['full_energy_conservation_alone_implies_bare_covariance']
    assert result['covariance_requires_additive_energy_and_invariant_joint_input']
    assert not result['spatial_reference_522_523_generates_time_asymmetry']
    assert not result['natural_preparation_of_1002_auxiliary_certified']
    assert not result['full_goal_completed']
    note=STAGE/'research_note_1004_working.md'
    for term in ('正式仍1003／累计3785','新增物理试验0','整体目标未完成','H₀+V'):
        assert term in note.read_text('utf-8'),term
    files=[note,HERE/'resource_preparation_adoption_v1.md',HERE/'resource_preparation_audit.py',
        HERE/'resource_preparation_audit_results.json',HERE/'drafts/adoption_decision.md',
        HERE/'drafts/NEXT.md',HERE/'drafts/publish_resource_adoption.py',Path(__file__)]
    for pth in files:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    links=0
    for doc in [p for p in files if p.suffix=='.md']+nav:
        content=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    assert '001—1003轮共1003份' in nav[0].read_text('utf-8-sig')
    assert '231—1003的773份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '1004工作审计：任务资源与关系参考' in pth.read_text('utf-8-sig')
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1003）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    evidence=[STAGE/'1003/research_round_1003_checks.json',HERE/'drafts/STATUS.md']
    out=dict(date='2026-10-07',kind='working_resource_adoption_audit',
        all_audit_checks_passed=True,formal_reports=1003,cumulative_numbered_test_groups=3785,
        new_scientific_test_groups=0,historical_unique_files_verified=len(frozen),
        historical_manifest_evidence=layout,local_links_checked=links,
        app_goal_changed=False,joint_lifecycle_certified=False,full_goal_completed=False,
        visual_checks_performed=False,
        frozen_evidence_hashes={str(p.relative_to(ROOT)):sha(p) for p in evidence},
        audit_file_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})
    if not writing:
        before=read(TARGET)
        for key in ('frozen_evidence_hashes','audit_file_hashes'):assert before[key]==out[key],key
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
        ('frozen_evidence_hashes','audit_file_hashes')},ensure_ascii=False,indent=2))
'''
    paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    originals={p:p.read_bytes() for p in paths};updates={}
    for p,raw in originals.items():
        s=raw.decode('utf-8-sig').replace('\r\n','\n')
        title='## 1004工作审计：任务资源与关系参考';mark='## 1003：整体假说v2与共同解释审计'
        assert title not in s and s.count(mark)==1,p
        pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
        block=(title+'\n\n'
            +f'[工作报告]({pre}research_note_1004_working.md)按人口记录、关系参考与裸非协变操作区分供给义务；'
            +'直接复用421／447／457等结果。992给人口失配，未制造纯辅助、相位参考或控制；'
            +'关系任务不要求外部绝对相位，但实际准备和访问仍需交代。'
            +f'[机制补充]({pre}1004/resource_preparation_adoption_v1.md) · '
            +f'[审计结果]({pre}1004/resource_preparation_audit_results.json) · '
            +f'[核验]({pre}1004/drafts/resource_adoption_checks.json)。\n\n'
            +'正式仍1003／累计3785，新增物理试验0；整体未完成。供给分类已可采用，'
            +f'接[下一项]({pre}1004/drafts/NEXT.md)补有限强场的场态、探测记录与共同来源；'
            +'不扩建写入器，不改应用目标。\n\n')
        s=s.replace(mark,block+mark,1)
        if p==paths[-1]:
            anchor='### 当前取舍\n';assert s.count(anchor)==1
            s=s.replace(anchor,anchor+'\n1004工作审计后，正式1003／累计3785保持。'
                '任务准备分为人口、关系参考和非协变操作；不把能源失配直接等同全部操作资源，'
                '不额外要求外部绝对相位。下一项补有限强场中场态—探测—来源机制。\n',1)
            replacement=('4. 957数学假说v0.2保持；[1003整体v2](../../1003/overall_operation_hypothesis_v2.md)'
                '及[1004供给审计](../../research_note_1004_working.md)区分机制采用、任务资源和共同实现。'
                '下一项按声明有限域补场态、几何访问与探测，不把完整蒸发当当前门槛。')
            s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M)
            assert count==1
        if b'\r\n' in raw:s=s.replace('\n','\r\n')
        updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
    assert all(p.read_bytes()==raw for p,raw in originals.items()),'Concurrent navigation edit'
    with (HERE/'verify_resource_adoption.py').open('x',encoding='utf-8') as dest:dest.write(prefix+checks)
    for p,data in updates.items():p.write_bytes(data)
    print('Published resource working audit; formal count remains 1003.')


if __name__=='__main__':main()
