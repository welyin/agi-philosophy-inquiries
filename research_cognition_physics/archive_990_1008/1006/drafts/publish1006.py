"""Publish conventional superconductivity adoption, preserving frozen work."""
from pathlib import Path
import re

STAGE = Path(__file__).resolve().parents[2]
RESEARCH = STAGE.parent
HERE = STAGE / '1006'


def main():
    old = (STAGE / '1005/verify_round1005.py').read_text('utf-8')
    assert old.count('    import importlib.util') == 1
    prefix = old.split('    import importlib.util', 1)[0]
    for before, after in (
        ('mechanism update 1005', 'superconductivity adoption 1006'),
        ('research_round_1005_checks.json', 'research_round_1006_checks.json'),
        ('range(776,1005)', 'range(776,1006)')):
        assert prefix.count(before) == 1, before
        prefix = prefix.replace(before, after, 1)
    checks = r'''
    import importlib.util
    spec=importlib.util.spec_from_file_location('mechanism1006',HERE/'superconductivity_adoption_audit.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'superconductivity_adoption_results.json')
    assert core.run()==saved
    assert saved['round']==1006 and saved['kind']=='mechanism_adoption_audit'
    assert saved['all_audit_checks_passed'] and saved['fresh_physical_test_groups']==0
    assert saved['cumulative_test_groups']==3786
    assert tuple(saved['evidence_obligations'])==(
        '参与者与过程','中间机制','资源、记录与反作用','输入与范围')
    assert saved['evidence_obligation_count']==len(saved['evidence_obligations'])==4
    assert all(saved['evidence_obligations'].values())
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('new_cognitive_axioms','microscopic_pairing_derived_from_cognition',
        'real_material_Tc_predicted','full_quantum_record_lifecycle_certified',
        'physical_unification_certified','full_goal_completed','visual_checks_performed'):
        assert saved[key] is False,key
    note=STAGE/'research_note_1006.md'
    prose=note.read_text('utf-8-sig')
    assert prose.count('$$')>0 and prose.count('$$')%2==0
    for term in ('整体目标未完成','累计3786','新增物理试验组0','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'conventional_superconductivity_adoption_v1.md',
        HERE/'superconductivity_adoption_audit.py',
        HERE/'superconductivity_adoption_results.json',Path(__file__),
        HERE/'drafts/adoption_decision.md',HERE/'drafts/review_notes.md',
        HERE/'drafts/publish1006.py',STAGE/'1007/drafts/STATUS.md']
    assert len(newfiles)==len(set(newfiles))==9
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8-sig'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1006）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=1006)==list(range(1,1007))
    assert '001—1006轮共1006份' in nav[0].read_text('utf-8-sig')
    assert '231—1006的776份' in nav[5].read_text('utf-8-sig')
    for pth in nav:
        navigation=pth.read_text('utf-8-sig')
        assert navigation.count('## 1006：常规超导与集体相干')==1,pth
        assert navigation.index('## 1006：常规超导与集体相干')<navigation.index(
            '## 1005：辐射反馈与整体假说v2.1'),pth
    oldfiles=[STAGE/'1005/research_round_1005_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    assert prev['all_delivery_checks_passed']
    assert prev['cumulative_numbered_test_groups_from_1004']==3786
    out=dict(round=1006,date='2026-10-08',all_delivery_checks_passed=True,
        formal_reports=1006,fresh_test_groups=0,
        cumulative_numbered_test_groups_from_1005=prev['cumulative_numbered_test_groups_from_1004'],
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,mechanism_evidence_obligations=len(saved['evidence_obligations']),
        base_mechanism_version='v2.1',mechanism_adoption='P08_conventional_superconductivity',
        superconductivity_solver_run=False,new_cognitive_axioms=False,
        microscopic_pairing_derived_from_cognition=False,real_material_Tc_predicted=False,
        full_quantum_record_lifecycle_certified=False,
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
        title = '## 1006：常规超导与集体相干'
        mark = '## 1005：辐射反馈与整体假说v2.1'
        assert title not in text and text.count(mark) == 1, path
        pre = 'archive_764_/' if path.parent == RESEARCH else '../../' if path == paths[-1] else ''
        block = (title + '\n\n'
            + f'[1006报告]({pre}research_note_1006.md)补常规声子配对、集体相位刚性、规范不变电磁响应及实际关系相位；'
            + '宏观有序可以与其他变量中的耗散及经典记录共存，但零电阻不等于无成本生命周期。'
            + '正式1006／累计3786，新增物理试验0，整体未完成。\n\n'
            + f'[机制补充]({pre}1006/conventional_superconductivity_adoption_v1.md) · '
            + f'[采用审计]({pre}1006/superconductivity_adoption_audit.py) · '
            + f'[结果]({pre}1006/superconductivity_adoption_results.json) · '
            + f'[核验]({pre}1006/research_round_1006_checks.json)。'
            + f'本轮接[整体v2.1]({pre}1005/overall_operation_hypothesis_v2_1.md)的P08；'
            + '复用355中性凝聚态及1000材料接口，不计算现实Tc、完整相图或器件寿命。'
            + f'接[1007]({pre}1007/drafts/STATUS.md)补中微子产生—传播—探测。'
            + '应用目标、957和条件性空间接口保持。\n\n')
        text = text.replace(mark, block + mark, 1)
        if path == paths[0]:
            before = '001—1005轮共1005份'
            assert text.count(before) == 1
            text = text.replace(before, '001—1006轮共1006份', 1)
        if path == paths[4]:
            before = '当前正式1005／累计3786，1005已结项'
            assert text.count(before) == 1
            text = text.replace(before, '当前正式1006／累计3786，1006已结项', 1)
        if path == paths[5]:
            assert text.count('# 231—1005轮阶段成果总览') == 1
            assert text.count('231—1005的775份') == 1
            text = text.replace('# 231—1005轮阶段成果总览', '# 231—1006轮阶段成果总览', 1)
            text = text.replace('231—1005的775份', '231—1006的776份', 1)
        if path == paths[-1]:
            before = '## 六条共同协议：全局缺口对应与检验优先级（截至1005）'
            assert text.count(before) == 1
            text = text.replace(before, before.replace('1005', '1006'), 1)
            original_rows = dict(re.findall(r'^(\|C\d\d )([^\n]*)$', text, re.M))
            assert len(original_rows) == 27
            additions = {
                'C17': ('1006补常规超导的配对与电磁响应',
                        '采用常规支，不等于现实材料参数或全部量子相已核'),
                'C19': ('1006补有限制备及实际弱连接的关系配对相位',
                        '粒子数／电荷的关系资源不等同时间相干，超导体不自动提供任意任务参考'),
                'C26': ('1006采用常规超导配对、集体刚性和Meissner磁响应，区分集体有序与逻辑相干',
                        '动态噪声、现实Tc、器件寿命及完整量子记录生命周期未核')}
            for cid, (state, boundary) in additions.items():
                pattern = r'^\|' + cid + r' [^\n]*$'
                matches = re.findall(pattern, text, re.M)
                assert len(matches) == 1, cid
                cells = matches[0].split('|')
                assert len(cells) == 6
                if cid == 'C17':
                    cells[4] = cells[4].replace('宏观量子有序机制下一项补', '其他宏观量子相仍按各自适用域采用')
                cells[2] += '；1006'
                cells[3] += '；' + state
                cells[4] += '；' + boundary
                text, count = re.subn(pattern, lambda _: '|'.join(cells), text, count=1, flags=re.M)
                assert count == 1
            changed_rows = dict(re.findall(r'^(\|C\d\d )([^\n]*)$', text, re.M))
            assert set(changed_rows) == set(original_rows)
            for key, row in original_rows.items():
                if key not in ('|C17 ', '|C19 ', '|C26 '):
                    assert changed_rows[key] == row, key
            assert '354—355' in changed_rows['|C26 ']
            assert '1000配对热—力' in changed_rows['|C26 ']
            assert '1005采用物质—辐射输运及环境反馈' in changed_rows['|C26 ']
            anchor = '### 当前取舍\n'
            assert text.count(anchor) == 1
            text = text.replace(anchor, anchor + '\n正式1006／累计3786。常规超导的配对、刚性、'
                '规范不变响应与选择性记录已纳入解释，未把它们签收为同一器件的完整生命周期。'
                '停止超导求解器及材料Tc细化；下一项补中微子产生—传播—探测与介质味转换。'
                '355中性凝聚态、1000热—力和1005辐射反馈保持。\n', 1)
            replacement = ('4. 957数学假说v0.2保持；[1005整体v2.1](../../1005/overall_operation_hypothesis_v2_1.md)'
                '保供给、强场和辐射形成反馈；[1006补充](../../1006/conventional_superconductivity_adoption_v1.md)'
                '纳入常规超导与关系相位。下一项补中微子产生—传播—探测，'
                '不以现实Tc、完整噪声求解或全相图作为机制采用前提。')
            text, count = re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',
                lambda _: replacement, text, count=1, flags=re.M)
            assert count == 1
        if b'\r\n' in raw:
            text = text.replace('\n', '\r\n')
        updates[path] = (b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'') + text.encode('utf-8')
    assert all(path.read_bytes() == raw for path, raw in originals.items()), 'Concurrent navigation edit'
    with (HERE / 'verify_round1006.py').open('x', encoding='utf-8') as dest:
        dest.write(prefix + checks)
    for path, data in updates.items():
        assert path.read_bytes() == originals[path], ('Concurrent navigation edit', path)
        path.write_bytes(data)
    print('Published 1006 conventional superconductivity mechanism adoption.')


if __name__ == '__main__':
    main()
