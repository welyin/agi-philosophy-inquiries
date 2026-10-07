"""Publish the v2.2 explanatory deliverable, preserving frozen evidence."""
from pathlib import Path
import re

STAGE = Path(__file__).resolve().parents[2]
RESEARCH = STAGE.parent
HERE = STAGE / '1008'


def main():
    old = (STAGE / '1007/verify_round1007.py').read_text('utf-8')
    assert old.count('    import importlib.util') == 1
    prefix = old.split('    import importlib.util', 1)[0]
    for before, after in (
        ('neutrino observation adoption 1007', 'overall explanatory delivery 1008'),
        ('research_round_1007_checks.json', 'research_round_1008_checks.json'),
        ('range(776,1007)', 'range(776,1008)')):
        assert prefix.count(before) == 1, before
        prefix = prefix.replace(before, after, 1)
    checks = r'''
    import importlib.util
    spec=importlib.util.spec_from_file_location('completion1008',HERE/'overall_completion_audit.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'overall_completion_results.json')
    assert core.run()==saved
    assert saved['round']==1008 and saved['kind']=='overall_mechanism_completion_audit'
    assert saved['all_audit_checks_passed'] is True
    assert saved['fresh_physical_test_groups']==0 and saved['cumulative_test_groups']==3786
    assert saved['phenomenon_count']==16
    assert saved['explanation_obligation_count']==64
    assert saved['explanation_obligation_count']==4*saved['phenomenon_count']
    assert saved['interface_count']==9 and saved['requirement_count']==9
    assert saved['current_explanatory_deliverable_completed'] is True
    for key in ('new_cognitive_axioms','physical_unification_certified',
        'full_empirical_recovery_certified','common_realization_certified'):
        assert saved[key] is False,key
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    note=STAGE/'research_note_1008.md'
    prose=note.read_text('utf-8-sig')
    assert prose.count('$$')%2==0
    for term in ('累计3786','新增物理试验组0','386或425','物理统一尚未证明'):
        assert term in prose,term
    newfiles=[note,HERE/'overall_operation_hypothesis_v2_2.md',
        HERE/'completion_audit.md',HERE/'overall_completion_audit.py',
        HERE/'overall_completion_results.json',Path(__file__),
        HERE/'drafts/current_objective_20261008.md',HERE/'drafts/review_notes.md',
        HERE/'drafts/publish1008.py',HERE/'follow_on_questions.md']
    assert len(newfiles)==len(set(newfiles))==10
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8-sig'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1008）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=1008)==list(range(1,1009))
    assert '001—1008轮共1008份' in nav[0].read_text('utf-8-sig')
    assert '231—1008的778份' in nav[5].read_text('utf-8-sig')
    for pth in nav:
        navigation=pth.read_text('utf-8-sig')
        title='## 1008：整体认知操作假说v2.2与解释层验收'
        assert navigation.count(title)==1,pth
        assert navigation.index(title)<navigation.index('## 1007：中微子观测与介质味转换'),pth
    oldfiles=[STAGE/'1007/research_round_1007_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    assert prev['all_delivery_checks_passed']
    assert prev['cumulative_numbered_test_groups_from_1006']==3786
    out=dict(round=1008,date='2026-10-08',all_delivery_checks_passed=True,
        formal_reports=1008,fresh_test_groups=0,
        cumulative_numbered_test_groups_from_1007=prev['cumulative_numbered_test_groups_from_1006'],
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,phenomenon_count=saved['phenomenon_count'],
        explanation_obligation_count=saved['explanation_obligation_count'],
        interface_count=saved['interface_count'],requirement_count=saved['requirement_count'],
        mechanism_version='v2.2',current_explanatory_deliverable_completed=True,
        new_cognitive_axioms=False,physical_unification_certified=False,
        independent_physical_derivation_certified=False,
        full_empirical_recovery_certified=False,common_realization_certified=False,
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
        title = '## 1008：整体认知操作假说v2.2与解释层验收'
        mark = '## 1007：中微子观测与介质味转换'
        assert title not in text and text.count(mark) == 1, path
        pre = 'archive_764_/' if path.parent == RESEARCH else '../../' if path == paths[-1] else ''
        block = (title + '\n\n'
            + f'[1008报告]({pre}research_note_1008.md)与[整体假说v2.2]({pre}1008/overall_operation_hypothesis_v2_2.md)'
            + '整合16类主要现象、64项解释义务和9个跨部门接口，并按9项要求完成当前解释层验收。'
            + '当前解释层交付完成；完整物理统一、独立物理推导、全部现实共同预测及共同实现尚未完成。'
            + '正式1008／累计3786，新增物理试验0。\n\n'
            + f'[验收说明]({pre}1008/completion_audit.md) · '
            + f'[审计代码]({pre}1008/overall_completion_audit.py) · '
            + f'[结果]({pre}1008/overall_completion_results.json) · '
            + f'[核验]({pre}1008/research_round_1008_checks.json)。'
            + '各机制保参与者、作用、来源、记录和反作用及适用边界，不把条件性采用汇总为同一现实历史。'
            + f'后续见[待选问题]({pre}1008/follow_on_questions.md)，不自动开启重复细化或新轮次。'
            + '957及386或425等既有条件性空间接口保持。\n\n')
        text = text.replace(mark, block + mark, 1)
        if path == paths[0]:
            before = '001—1007轮共1007份'
            assert text.count(before) == 1
            text = text.replace(before, '001—1008轮共1008份', 1)
        if path == paths[4]:
            before = '当前正式1007／累计3786，1007已结项'
            assert text.count(before) == 1
            text = text.replace(before, '当前正式1008／累计3786，1008当前解释层交付已验收', 1)
        if path == paths[5]:
            assert text.count('# 231—1007轮阶段成果总览') == 1
            assert text.count('231—1007的777份') == 1
            text = text.replace('# 231—1007轮阶段成果总览', '# 231—1008轮阶段成果总览', 1)
            text = text.replace('231—1007的777份', '231—1008的778份', 1)
        if path == paths[-1]:
            current_mapping = text.split('## 六条共同协议：全局缺口对应与检验优先级（截至1007）', 1)[1].split('### 当前取舍', 1)[0]
            original_rows = re.findall(r'^\|C\d\d [^\n]*$', current_mapping, re.M)
            all_original_rows = re.findall(r'^\|C\d\d [^\n]*$', text, re.M)
            assert len(original_rows) == 27
            assert [row.split(' ', 1)[0] for row in original_rows] == [f'|C{i:02d}' for i in range(1, 28)]
            before = '## 六条共同协议：全局缺口对应与检验优先级（截至1007）'
            assert text.count(before) == 1
            text = text.replace(before, before.replace('1007', '1008'), 1)
            anchor = '### 当前取舍\n'
            assert text.count(anchor) == 1
            text = text.replace(anchor, anchor + '\n正式1008／累计3786。整体v2.2已整合16类现象、'
                '64项解释义务及9个跨部门接口，当前解释层交付按9项要求验收完成。'
                'C01—C27的原条件、成立范围和未核边界全部保持；物理统一尚未证明，'
                '独立推导、全部现实共同预测和共同实现仍未完成。后续按'
                '[待选问题](../../1008/follow_on_questions.md)明确下一项问题，'
                '不自动开启重复细化或新轮次。\n', 1)
            replacement = ('4. 957数学假说v0.2保持；[1008整体v2.2](../../1008/overall_operation_hypothesis_v2_2.md)'
                '整合1005辐射、1006超导和1007中微子观测机制，保原P13不对称及生成分支；'
                '[解释层验收](../../1008/completion_audit.md)完成当前交付，'
                '不代表独立推导、全部现实预测或共同实现已获认证。后续见'
                '[待选问题](../../1008/follow_on_questions.md)，不自动延长候选细化。')
            text, count = re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',
                lambda _: replacement, text, count=1, flags=re.M)
            assert count == 1
            assert re.findall(r'^\|C\d\d [^\n]*$', text, re.M) == all_original_rows
        if b'\r\n' in raw:
            text = text.replace('\n', '\r\n')
        updates[path] = (b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'') + text.encode('utf-8')
    assert all(path.read_bytes() == raw for path, raw in originals.items()), 'Concurrent navigation edit'
    with (HERE / 'verify_round1008.py').open('x', encoding='utf-8') as dest:
        dest.write(prefix + checks)
    for path, data in updates.items():
        assert path.read_bytes() == originals[path], ('Concurrent navigation edit', path)
        path.write_bytes(data)
    print('Published 1008 overall v2.2 explanatory delivery; physical unification is not certified.')


if __name__ == '__main__':
    main()
