"""Publish the audited finite-effective completion after the real goal receipt."""
from pathlib import Path
import hashlib
import json
import re
import sys
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
from organize_research_1044_1062 import text_edits
from research_layout import relocation_original_bytes
from run_research_motion_active import MotionLayout

BASE = ROOT / 'research_cognition_physics'
STAGE = BASE / 'archive_1106_'
SHARED = STAGE / '_shared'
sha = lambda b: hashlib.sha256(b).hexdigest()
originals, updates = {}, {}

def read(p):
    p = Path(p)
    originals[p] = p.read_bytes()
    return originals[p].decode('utf-8-sig')

def put(p, text):
    p = Path(p)
    if p not in originals:
        originals[p] = p.read_bytes() if p.exists() else None
    if originals[p] and b'\r\n' in originals[p]:
        text = text.replace('\r\n', '\n').replace('\n', '\r\n')
    updates[p] = text.encode('utf-8')

def putjson(p, obj):
    put(p, json.dumps(obj, ensure_ascii=False, indent=2) + '\n')

receipt = json.loads(read(SHARED / 'goal_completion_receipt.json'))
assert receipt['goal']['status'] == 'complete'
assert not (SHARED / 'stage_completion_checks.json').exists(), 'Already published'
expected = {
 'stage_completion_audit.md': '91d4e4ff1db8ad64f7ffaeaf997efc50606ed3e1b7275cd4ac0a00469170379d',
 'goal_completion_adversarial_review.md': '598c5ac3f20dc0d9099fc4773bbb21784eac19961aadee48fcb4fc0150cb9c57',
 'model_dependency_summary.md': '6537a681b0947db83cbdec79c60ba6b6dd1f98b1a21f2eabd02b8572790ec158',
}
for name, h in expected.items():
    read(SHARED / name)
    assert sha(originals[SHARED / name]) == h, name
for name in ['stage_completion_audit.md', 'goal_completion_adversarial_review.md']:
    assert 'PASS_FOR_DECLARED_FINITE_EFFECTIVE_REALIZATION' in originals[SHARED / name].decode('utf-8-sig')

summary = ('**1106起的本阶段目标已完成，工具状态complete。** '
 '成果为[声明有限有效域中的共同实现](archive_1106_/_shared/stage_completion.md)：'
 '同一完整过程连接实际惯性参考、原CO位置、内部钟、未知Q及任意参考、全部已纳入影响和物质反冲，'
 '共同变换、四动量守恒及真实端口／后续相容。'
 '[1115](archive_1106_/research_note_1115.md)补齐对象和任务菜单证书，两份独立目标审计通过。'
 '1106—1112、1114—1115共九份报告、59组构造检查，科学校准增量0，累计3860。'
 '输入及剩余自由见[依赖总结](archive_1106_/_shared/model_dependency_summary.md)。'
 '未签所有原生材料的完整操作覆盖、经验真实性或完整物理SR；1113保留为可选加强候选。'
 '本目标停止自动续轮，未新建下一目标。历史复算使用 scripts/run_research_motion_active.py。')

lp = BASE / '_migration/motion_active_1106/layer.json'
layer = json.loads(read(lp))
marker = '## 历史导航（截至1105，保留原研究时点）'
for e in layer['entries']:
    p = ROOT / e['destination']
    if p.name not in {'README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md'}:
        continue
    before = read(p)
    baseline = relocation_original_bytes(originals[p], e)
    changed, count = re.subn(r'\*\*已完成1106—1112、1114—1115，共九份报告。\*\*[^\r\n]+', summary, before, count=1)
    assert count == 1
    changed = changed.replace('## 当前目标：实际狭义相对论的共同连接（1106起）',
                              '## 1106起阶段结项：有限有效域的实际狭义相对论共同连接', 1)
    changed = changed.replace('**新目标已按用户授权开启，工具状态active。**',
                              '**本阶段已完成原目标限定的共同实现，工具状态complete。**', 1)
    changed = changed.replace('推进次序：先建立实际钟与世界线的共同连接，再检查同一参考族的实际实现、全部可辨影响及物质交互的一致性。',
                              '已完成实际钟与世界线、同一参考族、全部纳入影响及物质交互的共同连接，适用域和材料输入按结项说明保留。', 1)
    changed = changed.replace('以下结项栏及历史“当前／下一步”保留其原研究时点；本栏和本阶段目标给出新的执行方向。',
                              '以下旧结项栏及历史“当前／下一步”保留原研究时点；当前完成状态与范围以上栏及1106起结项记录为准。', 1)
    assert before.split(marker, 1)[1] == changed.split(marker, 1)[1]
    put(p, changed)
    e.update(current_sha256=sha(updates[p]), current_bytes=len(updates[p]),
             edits=text_edits(baseline.decode('utf-8-sig'),updates[p].decode('utf-8')))
    assert relocation_original_bytes(updates[p], e) == baseline
layer.update(scope='Completed declared finite-effective realization stage; original navigation and paper baselines preserved',
             goal_status='complete', latest_completed_round=1115, new_formal_reports=9)
putjson(lp, layer)

p = STAGE / 'README.md'
t = read(p)
t = t.replace('2026-10-10启动。当前目标active，承接',
              '2026-10-10启动，2026-10-11完成声明有限有效域共同实现；目标工具状态complete。承接',1)
t, count = re.subn(r'已完成1106—1112、1114—1115九份报告[^\r\n]+',
 '已完成1106—1112、1114—1115九份报告、59组构造检查；科学校准增量0，累计3860。'
 '[阶段验收](_shared/stage_completion.md)及两份独立目标审计通过；[目标回执](_shared/goal_completion_receipt.json)为complete。'
 '结论限于声明有效域的共同实现，未签完整物理SR或全部原生材料的完整C。1113仍为可选加强候选，不计正式轮次。',t,count=1)
assert count == 1
t, count = re.subn(r'下一项按\[1114后目标范围审计\][^\r\n]+',
 '本目标的共同实现连接已验收；输入、职责合并及后续可选方向见[依赖总结](_shared/model_dependency_summary.md)。'
 '未开启下一目标。以下旧准入和每轮“下一步”均保留历史时点，不自动启动已完成研究。',t,count=1)
assert count == 1
t += '\n## 阶段验收资料\n\n- [主线验收与结论范围](_shared/stage_completion.md)\n- [独立原目标审计](_shared/stage_completion_audit.md)\n- [独立反方审阅](_shared/goal_completion_adversarial_review.md)\n- [输入、依赖消冗余与剩余自由](_shared/model_dependency_summary.md)\n- [冻结证据核验](_shared/stage_frozen_evidence.json) · [结项发布核验](_shared/stage_completion_checks.json)\n'
put(p,t)

p = STAGE / '文件索引.md'
t = read(p).replace('阶段验收继续，目标active。','阶段有限有效域验收通过，目标complete。',1)
t += '\n## 阶段完成\n\n- [验收范围与连接定理](_shared/stage_completion.md)\n- [逐项目标审计](_shared/stage_completion_audit.md) · [反方审阅](_shared/goal_completion_adversarial_review.md)\n- [依赖与剩余自由](_shared/model_dependency_summary.md)\n- [实际目标完成回执](_shared/goal_completion_receipt.json) · [最终检查](_shared/stage_completion_checks.json)\n'
put(p,t)

p = SHARED / 'stage_completion.md'
t = read(p).replace('这些阶段文件须全部定稿并通过主线检查后才请求工具标记完成。',
                    '这些阶段文件已全部定稿并经主线检查，目标工具已返回complete。',1)
put(p,t)

p = SHARED / 'dependency_ledger.json'
ledger = json.loads(read(p))
ledger['progress'].update(goal_complete=True, no_progress_streak=0,
    previous_goal_turn='progress: 1115 published; two independent original-goal audits accepted; declared-domain stage complete')
ledger['goal_status'] = 'complete'
ledger['stage_completion'] = dict(date='2026-10-11', scope='declared_finite_effective_complete_process_realization',
    mainline='research_cognition_physics/archive_1106_/_shared/stage_completion.md',
    two_independent_goal_audits=True, full_physical_SR=False, empirical_truth=False,
    all_native_material_C=False, generation_G_adopted=False, SR_S_adopted=False,
    receipt='research_cognition_physics/archive_1106_/_shared/goal_completion_receipt.json')
for c in ledger['open_connections']:
    c['stage_acceptance_status'] = 'satisfied_within_1115_declared_domain'
    c['stage_acceptance_source'] = 'research_cognition_physics/archive_1106_/_shared/stage_completion_audit.md'
    c['next'] = '原阶段责任已验收；更广覆盖仅按新声明域另行启动'
ledger['next_admission'] = dict(completed_round=False, automatic_next_round=False,
    path='research_cognition_physics/archive_1106_/_shared/model_dependency_summary.md',
    next='本目标完成；下一阶段方向待另行指定', candidate1113_optional=True)
ledger['dependency_consolidation'] = dict(new_axioms=0, logical_minimality_proved=False,
    material_packages=['K1106','R1111','K1112','R1114'],
    duties=['完整角色与准备','局部正时長完整响应','原端口与完整接续','整网参考与未来支持'],
    auxiliary_not_required_actual_vertices=['PT','B1109','E1110'],
    source='research_cognition_physics/archive_1106_/_shared/model_dependency_summary.md')
putjson(p,ledger)

for p, raw in originals.items():
    assert (p.read_bytes() if p.exists() else None)==raw, ('Concurrent change',p)
for p,raw in updates.items():
    p.write_bytes(raw)
compat=MotionLayout().verify_motion_layer()

checkpath = SHARED / 'stage_completion_checks.json'
files=[STAGE/'README.md',STAGE/'文件索引.md',*(SHARED/n for n in
    ['stage_completion.md','stage_completion_audit.md','goal_completion_adversarial_review.md','model_dependency_summary.md']),
    *(BASE/n for n in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md'])]
links=0
for p in files:
    txt=p.read_text(encoding='utf-8-sig')
    # Large historical navigation tails are preserved and covered by their old audits.
    if marker in txt: txt=txt.split(marker,1)[0]
    assert not re.search(r'[\x00-\x08\x0b\x0c\x0e-\x1f]',txt),p
    assert len(re.findall(r'^\s*\$\$\s*$',txt,re.M))%2==0,p
    for dest in re.findall(r'\]\(([^)]+)\)',txt):
        if '://' in dest or dest.startswith('#'): continue
        q=dest.strip('<>').split('#')[0]
        if not q: continue
        target=(p.parent/unquote(q)).resolve()
        assert target.exists() or target==checkpath.resolve(),(p,q)
        links+=1
evidence=json.loads((SHARED/'stage_frozen_evidence.json').read_text(encoding='utf-8'))
for a in evidence['rounds']:
    d=STAGE/str(a['round'])
    assert sha((d/'acceptance.json').read_bytes())==a['acceptance_sha256']
    acc=json.loads((d/'acceptance.json').read_text(encoding='utf-8-sig'))
    for name,h in acc['core_sha256'].items(): assert sha((d/name).read_bytes())==h
for name,h in expected.items(): assert sha((SHARED/name).read_bytes())==h
out=dict(schema='stage1106_completion_publication_v1',date='2026-10-11',goal_status='complete',
    acceptance_scope='declared_finite_effective_realization',reports=9,groups=59,
    full_physical_SR_claimed=False,new_scientific_calibration=0,scientific_count_total=3860,
    independent_audit_and_dependency_sha256=expected,all_round_core_hashes_preserved=True,
    compatibility=compat,markdown_files_checked=len(files),local_links_checked=links,
    historical_navigation_tails_preserved=True,new_goal_started=False,
    completion_receipt_sha256=sha((SHARED/'goal_completion_receipt.json').read_bytes()),
    completion_document_sha256=sha((SHARED/'stage_completion.md').read_bytes()))
checkpath.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,ensure_ascii=False))
