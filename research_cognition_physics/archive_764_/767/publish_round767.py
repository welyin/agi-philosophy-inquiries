"""Publish 767 with byte snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
TARGET = HERE/'round767_navigation_checks.json'
checked = core.read(HERE/'research_round_767_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes', 'preserved_draft_hashes'):
    for name, digest in checked[key].items():
        assert core.digest(HERE/name) == digest, name
paths = [ROOT/'README.md', RESEARCH/'README.md', RESEARCH/'research_direction.md',
         RESEARCH/'RESEARCH_STATE.md', HERE/'README.md',
         HERE/'spatial_premise_closure_audit.md', HERE/'three_dimensional_four_conditions_review.md']
before = {p: p.read_bytes() for p in paths}
folder = HERE/'navigation_before_round767_20261004'
assert not folder.exists() and not TARGET.exists()
summary = '**第767轮完成（同一自由物理态的正性）：** [共同平滑修正与Hadamard正态]({p}research_note_767.md)在原紧初片无稳定子背景上，对完整耦合线性玻色场共同满足正性、规范、CCR和短距离条件。态选择不唯一，未证原实际制备、全来源或相互作用。三组、十八式通过，最新767／3488，1613份编号科学文件、3714份保护证据。[核验]({p}research_round_767_checks.json)、[条件账]({p}unified_physics_condition_ledger_767.md)。'
order = '**当前执行顺序（767后，优先于下方历史安排）：** 接[768同阶全部门来源与共同Ward]({p}round768_drafts/STATUS.md)，复用734—735，核同一原作用、物理商及玻色/费米/ghost来源；不继续态的辅助常数。[范围审计]({p}round767_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及目标保持。'
planned = {}
for p, raw in before.items():
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl = '\r\n' if b'\r\n' in raw else '\n'
    text = raw.decode(enc).replace('\r\n', '\n')
    assert '**第767轮完成（同一自由物理态的正性）：**' not in text
    prefix = 'research_cognition_physics/archive_231_/' if p.parent == ROOT else 'archive_231_/' if p.parent == RESEARCH else ''
    if p in paths[:5]:
        head, rest = text.split('\n\n', 1)
        text = head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name == 'spatial_premise_closure_audit.md':
        assert '\n## 413.' not in text
        text += '\n\n## 413. 同一完整自由物理商的Hadamard正性\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 318.' not in text
        text += '\n\n## 318. 旧空间合同保持，自由物理态与完整来源\n\n'+summary.format(p=prefix)+'\n'
    if p.name == 'research_direction.md':
        assert '完成231—766轮。' in text and '最新科学轮次与检查数为766／3485' in text
        text = text.replace('完成231—766轮。', '完成231—767轮。').replace('最新科学轮次与检查数为766／3485', '最新科学轮次与检查数为767／3488')
    if p == HERE/'README.md':
        text += '\n\n## 第767轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|767|[共同平滑修正与Hadamard正态](research_note_767.md)|[代码](joint_physical_hadamard_positivity.py)、[结果](joint_physical_hadamard_positivity_results.json)、[核验](research_round_767_checks.json)|\n'
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
    temp = p.with_name(p.name+'.round767.tmp')
    with temp.open('xb') as f:
        f.write(raw)
    os.replace(temp, p)
report = dict(date='2026-10-04', latest_round=767, cumulative_tests=3488,
              numbered_scientific_files=1613, unique_protected_evidence_files=3714,
              navigation_files=7, navigation_links=links, broken_links=0,
              navigation_snapshots_preserved=True, complete_report_published=True,
              existing_results_preserved=True, next_round=768,
              active_goal_unchanged=True, stage_complete=False, all_checks_passed=True,
              published_navigation_hashes={p.relative_to(ROOT).as_posix(): core.digest(p) for p in paths})
with TARGET.open('x', encoding='utf8', newline='\n') as f:
    f.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
print(json.dumps({k: report[k] for k in ('latest_round', 'navigation_links', 'all_checks_passed')}))
