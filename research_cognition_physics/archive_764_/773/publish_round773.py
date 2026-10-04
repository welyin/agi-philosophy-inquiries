"""Publish 773 with byte snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
TARGET = HERE/'round773_navigation_checks.json'
checked = core.read(HERE/'research_round_773_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes', 'preserved_draft_hashes'):
    for name, digest in checked[key].items():
        assert core.digest(HERE/name) == digest, name
paths = [ROOT/'README.md', RESEARCH/'README.md', RESEARCH/'research_direction.md',
         RESEARCH/'RESEARCH_STATE.md', HERE/'README.md',
         HERE/'spatial_premise_closure_audit.md', HERE/'three_dimensional_four_conditions_review.md']
before = {p: p.read_bytes() for p in paths}
folder = HERE/'navigation_before_round773_20261004'
assert not folder.exists() and not TARGET.exists()
summary = '**第773轮完成（原物质参考与局部BRST）：** [局部规范左逆与形式修复]({p}research_note_773.md)原参考、Higgs及颜色电磁曲率给有限微分左逆；允许逆参考jet的正则片系数类中可逐阶求局部原始元。旧UV类别、共同QME和全局延伸未证。三组、十四式通过，最新773／3505，1631份编号科学文件、3806份保护证据。[核验]({p}research_round_773_checks.json)、[条件账]({p}unified_physics_condition_ledger_773.md)。'
order = '**当前执行顺序（773后，优先于下方历史安排）：** 接[774局部反项与来源首项]({p}round774_drafts/STATUS.md)，核原有限阶插入和771处方是否能共同实现；不重复参考秩或自由配对。[范围审计]({p}round773_drafts/scope_and_dedup_review.md)。补丁与系数类新增条件明示，原物理态、旧空间及目标保持。'
planned = {}
for p, raw in before.items():
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl = '\r\n' if b'\r\n' in raw else '\n'
    text = raw.decode(enc).replace('\r\n', '\n')
    assert '**第773轮完成（原物质参考与局部BRST）：**' not in text
    prefix = 'research_cognition_physics/archive_231_/' if p.parent == ROOT else 'archive_231_/' if p.parent == RESEARCH else ''
    if p in paths[:5]:
        head, rest = text.split('\n\n', 1)
        text = head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name == 'spatial_premise_closure_audit.md':
        assert '\n## 419.' not in text
        text += '\n\n## 419. 原物质参考与局部规范左逆\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 324.' not in text
        text += '\n\n## 324. 旧空间合同保持，局部参考与反项类别\n\n'+summary.format(p=prefix)+'\n'
    if p.name == 'research_direction.md':
        assert '完成231—772轮。' in text and '最新科学轮次与检查数为772／3502' in text
        text = text.replace('完成231—772轮。', '完成231—773轮。').replace('最新科学轮次与检查数为772／3502', '最新科学轮次与检查数为773／3505')
    if p == HERE/'README.md':
        text += '\n\n## 第773轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|773|[局部规范左逆与形式修复](research_note_773.md)|[代码](joint_material_local_brst.py)、[结果](joint_material_local_brst_results.json)、[核验](research_round_773_checks.json)|\n'
    planned[p] = text.replace('\n', nl).encode(enc)
links = 0
for p, raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(), (p, link)
        links += 1
assert all(p.read_bytes() == raw for p, raw in before.items())
folder.mkdir(exist_ok=False)
manifest = {}
for i, (p, raw) in enumerate(before.items()):
    name = f'{i}_{p.name}'
    with (folder/name).open('xb') as f:
        f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()] = dict(snapshot=name, sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x', encoding='utf8') as f:
    json.dump(manifest, f, ensure_ascii=False, indent=2)
for p, raw in planned.items():
    assert p.read_bytes() == before[p]
    temp = p.with_name(p.name+'.round773.tmp')
    with temp.open('xb') as f:
        f.write(raw)
    os.replace(temp, p)
report = dict(date='2026-10-04', latest_round=773, cumulative_tests=3505,
              numbered_scientific_files=1631, unique_protected_evidence_files=3806,
              navigation_files=7, navigation_links=links, broken_links=0,
              navigation_snapshots_preserved=True, complete_report_published=True,
              existing_results_preserved=True, next_round=774,
              active_goal_unchanged=True, stage_complete=False, all_checks_passed=True,
              published_navigation_hashes={p.relative_to(ROOT).as_posix(): core.digest(p) for p in paths})
with TARGET.open('x', encoding='utf8', newline='\n') as f:
    f.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
print(json.dumps({k: report[k] for k in ('latest_round', 'navigation_links', 'all_checks_passed')}))
