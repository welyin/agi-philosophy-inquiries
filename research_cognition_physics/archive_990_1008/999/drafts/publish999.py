"""Publish a mechanism synthesis; do not count document audit as physics tests."""
from pathlib import Path
import re

STAGE=Path(__file__).resolve().parents[2]
RESEARCH=STAGE.parent
HERE=STAGE/'999'


def main():
    old=(STAGE/'998/verify_round998.py').read_text('utf-8')
    prefix=old.split('    import importlib.util',1)[0]
    prefix=prefix.replace('Verify finite receiver 998','Verify mechanism synthesis 999')
    prefix=prefix.replace('research_round_998_checks.json','research_round_999_checks.json')
    prefix=prefix.replace('range(776,998)','range(776,999)')
    checks=r'''
    import importlib.util
    spec=importlib.util.spec_from_file_location('audit999',HERE/'mechanism_adoption_audit.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'mechanism_adoption_results.json')
    assert core.run()==saved
    assert saved['all_audit_checks_passed'] and saved['fresh_physical_test_groups']==0
    assert saved['cumulative_test_groups']==3782
    for key in ('new_physics_theorem','joint_lifecycle_certified','all_physical_phenomena_explained',
                'full_goal_completed','coverage_is_physical_proof','old_numerical_models_rerun',
                'mathematical_axioms_957_changed','application_goal_changed'):
        assert not saved[key],key
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert [r['id'] for r in saved['conditions']]==[f'C{i:02d}' for i in range(1,28)]
    assert all(not r['newly_proved_closed'] for r in saved['conditions'])
    assert len(saved['coverage'])==16 and len(saved['adoption_changes'])==3
    assert all(not r['certified_common_history'] for r in saved['non_joinable_evidence'])
    note=STAGE/'research_note_999.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==4
    assert re.findall(r'\\tag\{(\d+)\}',prose)==['1','2']
    for term in ('整体目标未完成','新增物理试验组0','累计3782','不是物理验证','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'overall_operation_hypothesis_v1.md',HERE/'mechanism_adoption_audit.py',
        HERE/'mechanism_adoption_results.json',Path(__file__),HERE/'drafts/adoption_decision.md',
        HERE/'drafts/review_notes.md',HERE/'drafts/publish999.py',STAGE/'1000/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至999）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    goal=(STAGE/'990/drafts/current_goal_20261007.txt').read_text('utf-8').strip()
    direction=nav[1].read_text('utf-8-sig').replace('\r\n','\n')
    assert goal.split('\n',1)[1].strip() in direction
    links=0
    for doc in [pth for pth in newfiles if pth.suffix=='.md']+nav:
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
    assert sorted(n for n in nums if n<=999)==list(range(1,1000))
    assert '001—999轮共999份' in nav[0].read_text('utf-8-sig')
    assert '231—999的769份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '999：选择性组织与整体假说' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'998/research_round_998_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=999,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=999,research_kind='mechanism_synthesis',fresh_test_groups=0,
        cumulative_numbered_test_groups_from_998=prev['cumulative_numbered_test_groups_from_997'],
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,condition_ids_accounted_for=27,phenomenon_groups_documented=16,
        condition_coverage_is_physical_proof=False,new_physics_theorem=False,
        joint_lifecycle_certified=False,full_goal_completed=False,
        visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(pth.relative_to(ROOT)):sha(pth) for pth in oldfiles},
        new_scientific_and_entry_files={str(pth.relative_to(ROOT)):sha(pth) for pth in newfiles})
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
    target=HERE/'verify_round999.py'
    assert not target.exists()
    paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    originals={p:p.read_bytes() for p in paths};updates={}
    for p,raw in originals.items():
        s=raw.decode('utf-8-sig').replace('\r\n','\n')
        title='## 999：选择性组织与整体假说';mark='## 998：组织形成与有限接收者'
        assert title not in s and s.count(mark)==1,p
        pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
        block=(title+'\n\n'
            +f'[999报告]({pre}research_note_999.md)采用选择性开放组织，直接复用旧写入／保持／热更新；'
            +'补普通Coulomb物质的统计稳定性，区分有效视界任务和高能完成。'
            +f'[整体假说v1]({pre}999/overall_operation_hypothesis_v1.md)整合16类现象及其输入和空缺。\n\n'
            +f'[审计结果]({pre}999/mechanism_adoption_results.json) · [核验]({pre}999/research_round_999_checks.json)。'
            +'正式999／累计3782；本轮新增物理试验0，整体未完成。不同旧准备未拼成共同生命周期；'
            +f'接[1000]({pre}1000/drafts/STATUS.md)补宏观相、集体模式和输运。应用目标及957保持。\n\n')
        s=s.replace(mark,block+mark,1)
        if p==paths[0]:
            assert s.count('001—998轮共998份')==1
            s=s.replace('001—998轮共998份','001—999轮共999份',1)
        if p==paths[4]:
            before='当前正式998／累计3782，998已结项';assert s.count(before)==1
            s=s.replace(before,'当前正式999／累计3782，999已结项',1)
        if p==paths[5]:
            assert s.count('231—998的768份')==1
            s=s.replace('# 231—998轮阶段成果总览','# 231—999轮阶段成果总览',1)
            s=s.replace('231—998的768份','231—999的769份',1)
        if p==paths[-1]:
            before='## 六条共同协议：全局缺口对应与检验优先级（截至998）'
            assert s.count(before)==1
            s=s.replace(before,before.replace('998','999'),1)
            rows={
                'C03':('|C03 事件身份、关系记录、访问|1、2、3、4；105—106、391、434、929、958—980、988、999|'
                       '999直接采用选择性交互与实际访问，保护、写入、接热和可读记录分开|'
                       '不同η、参考、环境和控制历史未接成共同生命周期；热输出不是纯空白|'),
                'C15':('|C15 表示、粒子谱与统计|1、3、4；942、947、950—951、981、999|'
                       '951角色不增基本物种；999将P981费米统计接到成熟Coulomb稳定性|'
                       '统计与三维仍物理输入；不是从认知或角色数推出SM谱|'),
                'C17':('|C17 质量、耦合、稳定物质结构|2、3、6；434、941—958、981、998—999|'
                       '原材料接口保留；999采用有明确非相对论工作域的普通物质规模线性能量下界|'
                       '下界不证明具体成键、相、形成、寿命或实际记录；共同材料匹配未核|'),
                'C25':('|C25 强引力、黑洞、紫外完成|1、4、6；990、993、999在明确任务内|'
                       '999区分GR外部／视界与半经典辐射任务，可按各自有效域采用成熟机制|'
                       '未核本项目完整黑洞恢复；奇点、完整蒸发及高能完成另列，不自动阻断有限任务|'),
                'C26':('|C26 宏观经典、流体等|1、3、4、6；350、354—355、758—763、973—977、998—999|'
                       '来源箱、热交换和选择性记录机制已保留；999定位相与集体输运为下一项|'
                       '物质稳定下界与粗粒化不能自动给固体刚性、声波或输运系数；相与时窗需明示|')}
            for cid,row in rows.items():
                s,count=re.subn(r'^\|'+cid+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M)
                assert count==1,cid
            anchor='### 当前取舍\n';assert s.count(anchor)==1
            s=s.replace(anchor,anchor+'\n正式999／累计3782。采用选择性组织、普通物质统计稳定性及视界范围拆分；'
                '整体v1覆盖主要现象的机制与空缺，覆盖不是物理证明。下一项宏观相与输运；'
                '不再扩建记录器或热设备，以下旧取舍保留当时范围。\n',1)
            replacement=('4. 957数学假说v0.2保持；[整体假说v1](../../999/overall_operation_hypothesis_v1.md)'
                '汇集实际量子过程、选择性组织、同源几何与有限资源。具体候选保B993及后续采用范围，'
                '下一项先补相、集体模式和输运的中间机制，不把全部材料制造或高能完成设为门槛。')
            s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M)
            assert count==1
        if b'\r\n' in raw:s=s.replace('\n','\r\n')
        updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
    assert all(p.read_bytes()==raw for p,raw in originals.items()),'Concurrent navigation edit'
    with target.open('x',encoding='utf-8') as dest:dest.write(prefix+checks)
    for p,data in updates.items():p.write_bytes(data)
    print('Published 999 synthesis to eight navigation files; physical test count unchanged.')


if __name__=='__main__':main()
