"""Publish a mechanism synthesis, without recording a new physical experiment."""
from pathlib import Path
import hashlib
import json
import re

STAGE = Path(__file__).resolve().parents[2]
RESEARCH = STAGE.parent
ROOT = RESEARCH.parent
HERE = STAGE / '990'
paths = [RESEARCH / n for n in ('README.md', 'research_direction.md', 'RESEARCH_STATE.md')] + [
    STAGE / n for n in ('README.md', '文件索引.md', '阶段成果总览.md', '跨阶段主题索引.md',
                      '_shared/notes/unified_physics_condition_ledger_current.md')]
original = {p: p.read_bytes() for p in paths}
ledger = original[paths[-1]].decode('utf-8-sig').replace('\r\n', '\n')
table = ledger.split('## 六条共同协议：全局缺口对应与检验优先级（截至989）', 1)[1].split('### 当前取舍', 1)[0]
rows = re.findall(r'^\|(C\d\d) ([^|]+)\|([^|]+)\|([^|]+)\|([^|]+)\|$', table, re.M)
assert len(rows) == 27
groups = {
    'Q': [1, 2, 3, 4], 'X': [5, 6, 7, 8, 9, 10],
    'M': [14, 15, 16, 17, 21], 'G': [11, 12, 13, 18, 22],
    'T': [19, 23, 26], 'U': [24, 25], 'E': [20, 27]}
assert sorted(n for ns in groups.values() for n in ns) == list(range(1, 28))
status = {
    1: 'inherited_conditional_reconstruction', 2: 'inherited_partial_interfaces',
    3: 'record_and_poststate_mechanism', 4: 'specified_dynamics_input',
    5: 'conditional_propagation', 6: 'effective_lorentz_input',
    7: 'conditional_dimension_or_explicit_input', 8: 'conditional_coordinates_or_explicit_input',
    9: 'internal_reference_partial', 10: 'proposed_common_geometry_input',
    11: 'constraints_at_declared_order', 12: 'einstein_action_input',
    13: 'route_specific_entropy_conditions', 14: 'gauge_group_and_connection_input',
    15: 'representations_and_species_input', 16: 'conditional_anomaly_checks',
    17: 'stability_and_material_matching_input', 18: 'same_source_partial_interfaces',
    19: 'preparation_and_boundary_input', 20: 'conditional_matching_and_errors',
    21: 'specified_sm_parent_partial_recovery', 22: 'backreaction_partial_interfaces',
    23: 'finite_arrow_with_preparation_input', 24: 'unresolved_cosmological_mechanisms',
    25: 'range_extension_not_default_gate', 26: 'task_dependent_classical_reports',
    27: 'inherited_distinguishable_predictions'}
evidence = [
    STAGE / '957/drafts/unified_operation_hypotheses_v0_2.md',
    STAGE / 'research_note_961.md', STAGE / 'research_note_977.md',
    STAGE / '981/drafts/common_parent_contract_v1.md', STAGE / 'research_note_984.md',
    STAGE / '985/drafts/research_note_985_working.md',
    STAGE / '987/drafts/research_note_987_working.md', STAGE / 'research_note_989.md',
    STAGE / '989/drafts/mechanism_map_v0_5.md',
    STAGE / '989/shared_record_quantum_action_results.json',
    STAGE / '989/research_round_989_checks.json', HERE / 'drafts/STATUS.md']
inventory = {
    'round': 990, 'kind': 'mechanism_synthesis_not_physics_validation',
    'fresh_physical_test_groups': 0, 'prior_cumulative_test_groups': 3774,
    'goal_complete': False, 'all_physical_phenomena_explained': False,
    'app_goal_edited_by_assistant': False, 'mathematical_axioms_957_changed': False,
    'proposal': 'B: common geometric comparison, representation-dependent internal transport, same-source feedback',
    'proposal_status': 'effective_physical_closure_candidate_not_cognitive_derivation',
    'mechanism_version': 'v0.3', 'groups': groups,
    'conditions': [dict(id=key, title=title.strip(),
        group=next(g for g, ns in groups.items() if int(key[1:]) in ns),
        status=status[int(key[1:])], inherited_evidence=ev.strip(),
        previous_state=st.strip(), previous_test_or_boundary=test.strip(),
        newly_proved_closed=False) for key, title, ev, st, test in rows],
    'independent_input_families': ['structure', 'dynamics_and_matching',
                                 'state_and_boundary', 'operations_and_approximation'],
    'open_explanatory_mechanisms': ['early_nonequilibrium_boundary', 'dark_matter_identity',
        'dark_energy_and_lambda', 'inflation_if_adopted', 'matter_antimatter_asymmetry'],
    'not_automatic_gates': ['unique_hamiltonian', 'uv_completion',
        'exact_microscopic_continuum', 'all_instrument_manufacturing'],
    'source_hashes': {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in evidence},
    'external_sources_checked': [
        'https://arxiv.org/abs/gr-qc/9405057', 'https://arxiv.org/abs/0903.5082',
        'https://arxiv.org/abs/0810.2712']}
dest = HERE / 'mechanism_inventory.json'
assert not dest.exists()
dest.write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

# Inherit the existing hash/layout verification, without rerunning old experiments.
prefix = (STAGE / '989/verify_round989.py').read_text('utf-8').split(
    "    result=read(HERE/'shared_record_quantum_action_results.json')", 1)[0]
prefix = prefix.replace('Delivery checks for 989.', 'Delivery checks for mechanism synthesis 990.')
prefix = prefix.replace('research_round_989_checks.json', 'research_round_990_checks.json')
prefix = prefix.replace('range(776,989)', 'range(776,990)')
checks = r'''
    inventory=read(HERE/'mechanism_inventory.json')
    assert inventory['round']==990 and inventory['fresh_physical_test_groups']==0
    assert not inventory['goal_complete'] and not inventory['all_physical_phenomena_explained']
    assert not inventory['mathematical_axioms_957_changed']
    assert [r['id'] for r in inventory['conditions']]==[f'C{i:02d}' for i in range(1,28)]
    assert sorted(n for ns in inventory['groups'].values() for n in ns)==list(range(1,28))
    for row in inventory['conditions']:
        assert not row['newly_proved_closed']
        assert int(row['id'][1:]) in inventory['groups'][row['group']]
    for rel,digest in inventory['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    note=STAGE/'research_note_990.md'
    prose=note.read_text('utf-8')
    assert prose.count('$$')==6
    for term in ('整体目标未完成','新增数值实验0','累计试验组仍3774',
                 '386或425是替代路线','不是物理自洽性证明','暗物质是什么'):
        assert term in prose,term
    newfiles=[note,HERE/'mechanism_inventory.json',Path(__file__),
        HERE/'drafts/current_goal_20261007.txt',HERE/'drafts/publish990.py',STAGE/'991/drafts/STATUS.md']
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至990）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    goal=(HERE/'drafts/current_goal_20261007.txt').read_text('utf-8').strip()
    direction=nav[1].read_text('utf-8-sig').replace('\r\n','\n')
    assert goal.split('\n',1)[1].strip() in direction
    links=0
    for doc in [p for p in newfiles if p.suffix=='.md']+nav:
        content=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if m and p.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=990)==list(range(1,991))
    assert '001—990轮共990份' in nav[0].read_text('utf-8-sig')
    assert '231—990的760份' in nav[5].read_text('utf-8-sig')
    for p in nav:assert '990：整体假说与联合解释优先' in p.read_text('utf-8-sig')
    oldfiles=[STAGE/'989/research_round_989_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=990,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=990,research_kind='mechanism_synthesis',fresh_test_groups=0,
        cumulative_numbered_test_groups_from_989=prev['cumulative_numbered_test_groups_from_988'],
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,condition_ids_accounted_for=27,
        condition_coverage_is_physical_proof=False,new_physics_theorem=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        user_goal_revision_synced_to_documents=True,
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
target = HERE / 'verify_round990.py'
assert not target.exists()
target.write_text(prefix + checks, encoding='utf-8')

updates = {}
heading = '## 990：整体假说与联合解释优先'
for p, raw in original.items():
    s = raw.decode('utf-8-sig').replace('\r\n', '\n')
    pre = 'archive_764_/' if p.parent == RESEARCH else '../../' if p == paths[-1] else ''
    previous = '## 989：共享证据与不可替代的量子作用'
    assert heading not in s and s.count(previous) == 1
    block = (heading + '\n\n'
        + f'[990整体报告]({pre}research_note_990.md)按用户最新目标，先整合全部门机制与解释边界。'
        + '提出共同几何比较、内部表示输运和同源反作用的联合候选；物理闭合输入与认知推导分开，'
        + '保留宇宙边界、暗部门等尚无机制的空缺。'
        + f'[条件清单]({pre}990/mechanism_inventory.json) · [交付核验]({pre}990/research_round_990_checks.json)。'
        + '正式990；本轮新增数值／物理试验0，累计3774，整体目标未完成。\n\n'
        + f'用户追加的整体解释优先要求已同步至文档，见[目标快照]({pre}990/drafts/current_goal_20261007.txt)。'
        + '957数学假说与历史保持；本轮不续修单项经典化、器件或连续极限。'
        + f'接[991]({pre}991/drafts/STATUS.md)先完成传播、物质和几何反馈的共同依赖说明，再选择决定性检验。\n\n')
    s = s.replace(previous, block + previous, 1)
    s = s.replace('001—989轮共989份', '001—990轮共990份')
    if p == paths[1]:
        start = s.index('2026-10-07按用户“修改一下目标，然后开始”的授权修订。',
                        s.index('## 当前目标正文：'))
        end = s.index('## 952：共同交换核的局部时间来源', start)
        goal = (HERE / 'drafts/current_goal_20261007.txt').read_text('utf-8').strip()
        intro = ('2026-10-07同步用户在应用中更新的目标，包含最新追加的“先保证整个假说在直觉上解释是合理自洽的”。'
            '本次仅同步研究文档，没有再次修改应用目标。完整正文见'
            '[当前目标快照](archive_764_/990/drafts/current_goal_20261007.txt)；'
            '[962旧正文](archive_764_/962/drafts/revised_goal_20261007.txt)及历史报告保留。'
            '下列目标优先于冻结入口的旧下一步；数学假说957保持，990新增机制提案按其独立输入地位使用。\n\n')
        s = s[:start] + intro + '### ' + goal[2:] + '\n\n' + s[end:]
    if p == paths[2]:
        a = s.index('用户已授权修改目标并开始研究。')
        b = s.index('\n\n', a)
        s = s[:a] + ('用户已在应用目标中强调先补全整体认知操作假说、先检查直觉机制的合理自洽。'
            '[当前目标](research_direction.md#当前目标正文从认知操作压缩物理的独立假设)及'
            '[990目标快照](archive_764_/990/drafts/current_goal_20261007.txt)已同步；'
            '本次不再改应用目标，不创建任务或定时任务。保留共同模型及全部物理部门；'
            '细节验算须由整体机制的决定性问题触发。') + s[b:]
    if p == paths[4]:
        s = s.replace('当前正式989／累计3774，989已结项', '当前正式990／累计3774，990机制整合已交付')
    if p == paths[5]:
        s = s.replace('# 231—989轮阶段成果总览', '# 231—990轮阶段成果总览', 1)
        s = s.replace('231—989的759份', '231—990的760份', 1)
    if p == paths[-1]:
        s = s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至989）',
                      '## 六条共同协议：全局缺口对应与检验优先级（截至990）', 1)
        marker = '下表覆盖[原C01—C27定义]'
        intro = ('990将下表27项按Q/X/M/G/T/U/E七组机制共同整理，详见'
            '[整体报告](../../research_note_990.md)与[逐项清单](../../990/mechanism_inventory.json)。'
            '本次没有将任何条件新标为已证关闭；采用候选B的共同几何、内部比较与同源反馈作为联合解释框架，'
            '输入补齐与生成解释严格分开。先补解释链，停止自动深化单项技术候选。\n\n')
        assert s.count(marker) == 1
        s = s.replace(marker, intro + marker, 1)
        replacement = ('4. 957数学假说v0.2保持；[990整体机制提案](../../research_note_990.md)'
            '统一经典报告、量子后态、共同几何、内部比较和同源反馈的职责。'
            '宇宙边界与暗部门解释仍开放；先补联合机制和独立输入，不续修局部器件。')
        s, count = re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$', lambda _: replacement, s, count=1, flags=re.M)
        assert count == 1
    if b'\r\n' in raw:
        s = s.replace('\n', '\r\n')
    updates[p] = (b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'') + s.encode('utf-8')
assert all(p.read_bytes() == raw for p, raw in original.items()), 'Concurrent navigation edit'
for p, data in updates.items():
    p.write_bytes(data)
print('Published mechanism synthesis 990; physical experiment count remains 3774.')
