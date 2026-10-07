"""Publish 768 with byte snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
TARGET = HERE/'round768_navigation_checks.json'
checked = core.read(HERE/'research_round_768_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes', 'preserved_draft_hashes'):
    for name, digest in checked[key].items():
        assert core.digest(HERE/name) == digest, name
paths = [ROOT/'README.md', RESEARCH/'README.md', RESEARCH/'research_direction.md',
         RESEARCH/'RESEARCH_STATE.md', HERE/'README.md',
         HERE/'spatial_premise_closure_audit.md', HERE/'three_dimensional_four_conditions_review.md']
before = {p: p.read_bytes() for p in paths}
folder = HERE/'navigation_before_round768_20261004'
assert not folder.exists() and not TARGET.exists()
summary = '**第768轮完成（共同态与相对量子来源）：** [BRST扩展与完整玻色来源差]({p}research_note_768.md)保767同一物理态，平滑物理变化可保持辅助部门不变；原作用给有限联合守恒来源差，并接原形式首阶约束响应。绝对全部门Ward和有限强度闭合未证。三组、十八式通过，最新768／3491，1616份编号科学文件、3726份保护证据。[核验]({p}research_round_768_checks.json)、[条件账]({p}unified_physics_condition_ledger_768.md)。'
order = '**当前执行顺序（768后，优先于下方历史安排）：** 接[769绝对全部门来源与量子Ward]({p}round769_drafts/STATUS.md)，核同一原作用的局部EFT/BV规范化、离壳背景及有限项；不重复自由辅助分块。[范围审计]({p}round768_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及目标保持。'
planned = {}
for p, raw in before.items():
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl = '\r\n' if b'\r\n' in raw else '\n'
    text = raw.decode(enc).replace('\r\n', '\n')
    assert '**第768轮完成（共同态与相对量子来源）：**' not in text
    prefix = 'research_cognition_physics/archive_231_/' if p.parent == ROOT else 'archive_231_/' if p.parent == RESEARCH else ''
    if p in paths[:5]:
        head, rest = text.split('\n\n', 1)
        text = head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name == 'spatial_premise_closure_audit.md':
        assert '\n## 414.' not in text
        text += '\n\n## 414. 同一物理态的BRST扩展与完整相对来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 319.' not in text
        text += '\n\n## 319. 旧空间合同保持，同态来源与绝对规范化\n\n'+summary.format(p=prefix)+'\n'
    if p.name == 'research_direction.md':
        assert '完成231—767轮。' in text and '最新科学轮次与检查数为767／3488' in text
        text = text.replace('完成231—767轮。', '完成231—768轮。').replace('最新科学轮次与检查数为767／3488', '最新科学轮次与检查数为768／3491')
    if p == HERE/'README.md':
        text += '\n\n## 第768轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|768|[BRST扩展与完整玻色来源差](research_note_768.md)|[代码](joint_brst_relative_source.py)、[结果](joint_brst_relative_source_results.json)、[核验](research_round_768_checks.json)|\n'
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
    temp = p.with_name(p.name+'.round768.tmp')
    with temp.open('xb') as f:
        f.write(raw)
    os.replace(temp, p)
report = dict(date='2026-10-04', latest_round=768, cumulative_tests=3491,
              numbered_scientific_files=1616, unique_protected_evidence_files=3726,
              navigation_files=7, navigation_links=links, broken_links=0,
              navigation_snapshots_preserved=True, complete_report_published=True,
              existing_results_preserved=True, next_round=769,
              active_goal_unchanged=True, stage_complete=False, all_checks_passed=True,
              published_navigation_hashes={p.relative_to(ROOT).as_posix(): core.digest(p) for p in paths})
with TARGET.open('x', encoding='utf8', newline='\n') as f:
    f.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
print(json.dumps({k: report[k] for k in ('latest_round', 'navigation_links', 'all_checks_passed')}))
