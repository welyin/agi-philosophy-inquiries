"""Publish neutrino observation adoption while preserving frozen work."""
from pathlib import Path
import re

STAGE = Path(__file__).resolve().parents[2]
RESEARCH = STAGE.parent
HERE = STAGE / '1007'


def main():
    old = (STAGE / '1006/verify_round1006.py').read_text('utf-8')
    assert old.count('    import importlib.util') == 1
    prefix = old.split('    import importlib.util', 1)[0]
    for before, after in (
        ('superconductivity adoption 1006', 'neutrino observation adoption 1007'),
        ('research_round_1006_checks.json', 'research_round_1007_checks.json'),
        ('range(776,1006)', 'range(776,1007)')):
        assert prefix.count(before) == 1, before
        prefix = prefix.replace(before, after, 1)
    checks = r'''
    import importlib.util
    spec=importlib.util.spec_from_file_location('mechanism1007',HERE/'neutrino_adoption_audit.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'neutrino_adoption_results.json')
    assert core.run()==saved
    assert saved['round']==1007 and saved['kind']=='mechanism_adoption_audit'
    assert saved['all_audit_checks_passed'] and saved['fresh_physical_test_groups']==0
    assert saved['cumulative_test_groups']==3786 and saved['section']=='P13'
    assert tuple(saved['evidence_obligations'])==(
        '参与者','机制','资源、记录与反作用','输入与边界')
    assert saved['evidence_obligation_count']==len(saved['evidence_obligations'])==4
    assert all(saved['evidence_obligations'].values())
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('new_cognitive_axioms','masses_predicted','PMNS_fit_performed',
        'actual_source_detector_simulated','dynamic_W995_oscillation_certified',
        'physical_unification_certified','full_goal_completed','visual_checks_performed'):
        assert saved[key] is False,key
    note=STAGE/'research_note_1007.md'
    prose=note.read_text('utf-8-sig')
    assert prose.count('$$')%2==0
    for term in ('整体目标未完成','累计3786','新增物理试验组0','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'neutrino_observation_adoption_v1.md',
        HERE/'neutrino_adoption_audit.py',HERE/'neutrino_adoption_results.json',Path(__file__),
        HERE/'drafts/adoption_decision.md',HERE/'drafts/review_notes.md',
        HERE/'drafts/publish1007.py',STAGE/'1008/drafts/STATUS.md']
    assert len(newfiles)==len(set(newfiles))==9
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8-sig'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1007）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=1007)==list(range(1,1008))
    assert '001—1007轮共1007份' in nav[0].read_text('utf-8-sig')
    assert '231—1007的777份' in nav[5].read_text('utf-8-sig')
    for pth in nav:
        navigation=pth.read_text('utf-8-sig')
        assert navigation.count('## 1007：中微子观测与介质味转换')==1,pth
        assert navigation.index('## 1007：中微子观测与介质味转换')<navigation.index(
            '## 1006：常规超导与集体相干'),pth
    oldfiles=[STAGE/'1006/research_round_1006_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    assert prev['all_delivery_checks_passed']
    assert prev['cumulative_numbered_test_groups_from_1005']==3786
    out=dict(round=1007,date='2026-10-08',all_delivery_checks_passed=True,
        formal_reports=1007,fresh_test_groups=0,
        cumulative_numbered_test_groups_from_1006=prev['cumulative_numbered_test_groups_from_1005'],
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,mechanism_evidence_obligations=len(saved['evidence_obligations']),
        base_mechanism_version='v2.1',mechanism_adoption='P13_neutrino_observation',
        neutrino_solver_run=False,new_cognitive_axioms=False,masses_predicted=False,
        PMNS_fit_performed=False,actual_source_detector_simulated=False,
        dynamic_W995_oscillation_certified=False,
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
    paths = [RESEARCH / name for name in ('README.md', 'research_direction.md', 'RESEARCH_STATE.md')] + [
        STAGE / name for name in ('README.md', '文件索引.md', '阶段成果总览.md', '跨阶段主题索引.md',
                                 '_shared/notes/unified_physics_condition_ledger_current.md')]
    originals = {path: path.read_bytes() for path in paths}
    updates = {}
    for path, raw in originals.items():
        text = raw.decode('utf-8-sig').replace('\r\n', '\n')
        title = '## 1007：中微子观测与介质味转换'
        mark = '## 1006：常规超导与集体相干'
        assert title not in text and text.count(mark) == 1, path
        pre = 'archive_764_/' if path.parent == RESEARCH else '../../' if path == paths[-1] else ''
        block = (title + '\n\n'
            + f'[1007报告]({pre}research_note_1007.md)采用中微子产生—传播—探测及普通介质味转换机制；'
            + '区分质量分支辨认、包重叠、宏观平均与绝热转换，来源、反冲和记录共同入账。'
            + '正式1007／累计3786，新增物理试验0，整体未完成。\n\n'
            + f'[机制补充]({pre}1007/neutrino_observation_adoption_v1.md) · '
            + f'[采用审计]({pre}1007/neutrino_adoption_audit.py) · '
            + f'[结果]({pre}1007/neutrino_adoption_results.json) · '
            + f'[核验]({pre}1007/research_round_1007_checks.json)。'
            + f'补[整体v2.1]({pre}1005/overall_operation_hypothesis_v2_1.md)的P13观测链，'
            + '保539质量接口、994初始净荷及995—996可选生成分支；普通味振荡不等于净荷生成。'
            + '不预测质量或PMNS参数，不认证现实源—探测器全链。'
            + f'接[1008]({pre}1008/drafts/STATUS.md)整合整体v2.2并审计四项任务义务，停止振荡细化。'
            + '应用目标、957和条件性空间接口保持。\n\n')
        text = text.replace(mark, block + mark, 1)
        if path == paths[0]:
            before = '001—1006轮共1006份'
            assert text.count(before) == 1
            text = text.replace(before, '001—1007轮共1007份', 1)
        if path == paths[4]:
            before = '当前正式1006／累计3786，1006已结项'
            assert text.count(before) == 1
            text = text.replace(before, '当前正式1007／累计3786，1007已结项', 1)
        if path == paths[5]:
            assert text.count('# 231—1006轮阶段成果总览') == 1
            assert text.count('231—1006的776份') == 1
            text = text.replace('# 231—1006轮阶段成果总览', '# 231—1007轮阶段成果总览', 1)
            text = text.replace('231—1006的776份', '231—1007的777份', 1)
        if path == paths[-1]:
            before = '## 六条共同协议：全局缺口对应与检验优先级（截至1006）'
            assert text.count(before) == 1
            text = text.replace(before, before.replace('1006', '1007'), 1)
            original_rows = dict(re.findall(r'^(\|C\d\d )([^\n]*)$', text, re.M))
            assert len(original_rows) == 27
            additions = {
                'C15': ('1007区分传播质量分支与弱相互作用味标签，采用产生—传播—探测链',
                        '现实质量差、PMNS、源谱与探测响应仍输入，未推出粒子谱或混合参数'),
                'C17': ('1007复用539质量接口并补普通介质味转换，保994初始荷及995—996生成分支',
                        '普通味振荡不生成净轻子／重子荷，未认证动态W995振荡或现实源产额'),
                'C20': ('1007将源、反冲、传播与探测记录共同登记，区分相干损失和经典平均',
                        '测试背景有适用域，未认证现实全链共同演化或全部背景反作用')}
            for cid, (state, boundary) in additions.items():
                pattern = r'^\|' + cid + r' [^\n]*$'
                matches = re.findall(pattern, text, re.M)
                assert len(matches) == 1, cid
                cells = matches[0].split('|')
                assert len(cells) == 6
                cells[2] += '；1007'
                cells[3] += '；' + state
                cells[4] += '；' + boundary
                text, count = re.subn(pattern, lambda _: '|'.join(cells), text, count=1, flags=re.M)
                assert count == 1
            changed_rows = dict(re.findall(r'^(\|C\d\d )([^\n]*)$', text, re.M))
            assert set(changed_rows) == set(original_rows)
            for key, row in original_rows.items():
                if key not in ('|C15 ', '|C17 ', '|C20 '):
                    assert changed_rows[key] == row, key
                else:
                    old_cells = row.split('|')
                    new_cells = changed_rows[key].split('|')
                    assert len(old_cells) == len(new_cells)
                    assert all(new.startswith(old) for old, new in zip(old_cells, new_cells)), key
            anchor = '### 当前取舍\n'
            assert text.count(anchor) == 1
            text = text.replace(anchor, anchor + '\n正式1007／累计3786。中微子观测和介质味转换已补成熟机制；'
                'P13的不对称、初始荷与生成／洗出分支保持，普通味转换不替代净荷来源。'
                '停止新增振荡曲线或现实PMNS拟合；下一项1008整合整体v2.2，并审计参与者、'
                '中间机制、资源记录反作用、输入边界四项任务义务及跨部门相容性。\n', 1)
            replacement = ('4. 957数学假说v0.2保持；[1005整体v2.1](../../1005/overall_operation_hypothesis_v2_1.md)'
                '及[1006超导补充](../../1006/conventional_superconductivity_adoption_v1.md)保持；'
                '[1007补充](../../1007/neutrino_observation_adoption_v1.md)纳入中微子观测链，'
                '不替换原P13不对称与生成分支。下一项1008整合整体v2.2及四项任务义务审计，'
                '不继续细化振荡求解。')
            text, count = re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',
                lambda _: replacement, text, count=1, flags=re.M)
            assert count == 1
        if b'\r\n' in raw:
            text = text.replace('\n', '\r\n')
        updates[path] = (b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'') + text.encode('utf-8')
    assert all(path.read_bytes() == raw for path, raw in originals.items()), 'Concurrent navigation edit'
    with (HERE / 'verify_round1007.py').open('x', encoding='utf-8') as dest:
        dest.write(prefix + checks)
    for path, data in updates.items():
        assert path.read_bytes() == originals[path], ('Concurrent navigation edit', path)
        path.write_bytes(data)
    print('Published 1007 neutrino observation mechanism adoption.')


if __name__ == '__main__':
    main()
