"""Publish 764 with byte snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
TARGET = HERE/'round764_navigation_checks.json'
checked = core.read(HERE/'research_round_764_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes', 'preserved_draft_hashes'):
    for name, digest in checked[key].items():
        assert core.digest(HERE/name) == digest, name
paths = [ROOT/'README.md', RESEARCH/'README.md', RESEARCH/'research_direction.md',
         RESEARCH/'RESEARCH_STATE.md', HERE/'README.md',
         HERE/'spatial_premise_closure_audit.md', HERE/'three_dimensional_four_conditions_review.md']
before = {p: p.read_bytes() for p in paths}
folder = HERE/'navigation_before_round764_20261004'
assert not folder.exists() and not TARGET.exists()
summary = '**第764轮完成（领先物理量子分支）：** [共同线性物理相空间与内禀量子态]({p}research_note_764.md)由原完整约束右逆构造辛截面、物理CCR及正则正态，接同背景领先费米部门。微扰量子化与准备是明示输入；Hadamard全来源、实际记录及原Q连续匹配仍开放。三组、十五式通过，最新764／3479，1604份编号科学文件、3679份保护证据。[核验]({p}research_round_764_checks.json)、[条件账]({p}unified_physics_condition_ledger_764.md)。'
order = '**当前执行顺序（764后，优先于下方历史安排）：** 接[765物理短距离态与同阶完整来源]({p}round765_drafts/STATUS.md)，核耦合物理态、玻色二次来源及费米来源的共同Ward和任务字典；不继续投影或自由度扫描。[范围审计]({p}round764_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及总目标保持。'
planned = {}
for p, raw in before.items():
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl = '\r\n' if b'\r\n' in raw else '\n'
    text = raw.decode(enc).replace('\r\n', '\n')
    assert '**第764轮完成（领先物理量子分支）：**' not in text
    prefix = 'research_cognition_physics/archive_231_/' if p.parent == ROOT else 'archive_231_/' if p.parent == RESEARCH else ''
    if p in paths[:5]:
        head, rest = text.split('\n\n', 1)
        text = head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name == 'spatial_premise_closure_audit.md':
        assert '\n## 410.' not in text
        text += '\n\n## 410. 共同线性物理相空间与正态\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 315.' not in text
        text += '\n\n## 315. 旧空间合同保持，共同量子相空间与同阶来源\n\n'+summary.format(p=prefix)+'\n'
    if p.name == 'research_direction.md':
        assert '完成231—763轮。' in text and '最新科学轮次与检查数为763／3476' in text
        text = text.replace('完成231—763轮。', '完成231—764轮。').replace('最新科学轮次与检查数为763／3476', '最新科学轮次与检查数为764／3479')
    if p == HERE/'README.md':
        text += '\n\n## 第764轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|764|[共同线性物理相空间与内禀量子态](research_note_764.md)|[代码](joint_linear_physical_phase.py)、[结果](joint_linear_physical_phase_results.json)、[核验](research_round_764_checks.json)|\n'
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
    temp = p.with_name(p.name+'.round764.tmp')
    with temp.open('xb') as f:
        f.write(raw)
    os.replace(temp, p)
report = dict(date='2026-10-04', latest_round=764, cumulative_tests=3479,
              numbered_scientific_files=1604, unique_protected_evidence_files=3679,
              navigation_files=7, navigation_links=links, broken_links=0,
              navigation_snapshots_preserved=True, complete_report_published=True,
              existing_results_preserved=True, next_round=765,
              active_goal_unchanged=True, stage_complete=False, all_checks_passed=True,
              published_navigation_hashes={p.relative_to(ROOT).as_posix(): core.digest(p) for p in paths})
with TARGET.open('x', encoding='utf8', newline='\n') as f:
    f.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
print(json.dumps({k: report[k] for k in ('latest_round', 'navigation_links', 'all_checks_passed')}))
