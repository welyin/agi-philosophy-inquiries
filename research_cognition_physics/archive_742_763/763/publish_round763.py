"""Publish 763 with byte snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
TARGET = HERE/'round763_navigation_checks.json'
checked = core.read(HERE/'research_round_763_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes', 'preserved_draft_hashes'):
    for name, digest in checked[key].items():
        assert core.digest(HERE/name) == digest, name
paths = [ROOT/'README.md', RESEARCH/'README.md', RESEARCH/'research_direction.md',
         RESEARCH/'RESEARCH_STATE.md', HERE/'README.md',
         HERE/'spatial_premise_closure_audit.md', HERE/'three_dimensional_four_conditions_review.md']
before = {p: p.read_bytes() for p in paths}
folder = HERE/'navigation_before_round763_20261004'
assert not folder.exists() and not TARGET.exists()
summary = '**第763轮完成（条件性初始连接）：** [同代数受约束响应与几何因子边界]({p}research_note_763.md)将旧完整来源右逆提升为同一源代数的线性初值，保正有序矩；原背景上的诊断来源给非对易共形响应。独立几何因子且精确保全部未知态的态依赖赋态被排除；原连续来源映射、内禀几何及实际完整过程仍开放。三组、十三式通过，最新763／3476，1601份编号科学文件、3666份保护证据。[核验]({p}research_round_763_checks.json)、[条件账]({p}unified_physics_condition_ledger_763.md)。'
order = '**当前执行顺序（763后，优先于下方历史安排）：** 接[764同一线性物理相空间]({p}round764_drafts/STATUS.md)，核原Gauss资料、完整引力约束、内禀初值与关系任务的共同映射；不重复右逆或交换子常数。[范围审计]({p}round763_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及总目标保持。'
planned = {}
for p, raw in before.items():
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl = '\r\n' if b'\r\n' in raw else '\n'
    text = raw.decode(enc).replace('\r\n', '\n')
    assert '**第763轮完成（条件性初始连接）：**' not in text
    prefix = 'research_cognition_physics/archive_231_/' if p.parent == ROOT else 'archive_231_/' if p.parent == RESEARCH else ''
    if p in paths[:5]:
        head, rest = text.split('\n\n', 1)
        text = head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name == 'spatial_premise_closure_audit.md':
        assert '\n## 409.' not in text
        text += '\n\n## 409. 同代数初始响应与几何因子边界\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 314.' not in text
        text += '\n\n## 314. 旧空间合同保持，同一来源与受约束响应\n\n'+summary.format(p=prefix)+'\n'
    if p.name == 'research_direction.md':
        assert '完成231—762轮。' in text and '最新科学轮次与检查数为762／3473' in text
        text = text.replace('完成231—762轮。', '完成231—763轮。').replace('最新科学轮次与检查数为762／3473', '最新科学轮次与检查数为763／3476')
    if p == HERE/'README.md':
        text += '\n\n## 第763轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|763|[同代数受约束响应与几何因子边界](research_note_763.md)|[代码](joint_operator_constraint_lift.py)、[结果](joint_operator_constraint_lift_results.json)、[核验](research_round_763_checks.json)|\n'
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
    temp = p.with_name(p.name+'.round763.tmp')
    with temp.open('xb') as f:
        f.write(raw)
    os.replace(temp, p)
report = dict(date='2026-10-04', latest_round=763, cumulative_tests=3476,
              numbered_scientific_files=1601, unique_protected_evidence_files=3666,
              navigation_files=7, navigation_links=links, broken_links=0,
              navigation_snapshots_preserved=True, complete_report_published=True,
              existing_results_preserved=True, next_round=764,
              active_goal_unchanged=True, stage_complete=False, all_checks_passed=True,
              published_navigation_hashes={p.relative_to(ROOT).as_posix(): core.digest(p) for p in paths})
with TARGET.open('x', encoding='utf8', newline='\n') as f:
    f.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
print(json.dumps({k: report[k] for k in ('latest_round', 'navigation_links', 'all_checks_passed')}))
