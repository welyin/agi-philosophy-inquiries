"""Publish 774 with byte snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
TARGET = HERE/'round774_navigation_checks.json'
checked = core.read(HERE/'research_round_774_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes', 'preserved_draft_hashes'):
    for name, digest in checked[key].items():
        assert core.digest(HERE/name) == digest, name
paths = [ROOT/'README.md', RESEARCH/'README.md', RESEARCH/'research_direction.md',
         RESEARCH/'RESEARCH_STATE.md', HERE/'README.md',
         HERE/'spatial_premise_closure_audit.md', HERE/'three_dimensional_four_conditions_review.md']
before = {p: p.read_bytes() for p in paths}
folder = HERE/'navigation_before_round774_20261004'
assert not folder.exists() and not TARGET.exists()
summary = '**第774轮完成（来源与量子作用的条件性连接）：** [保首项及共同插入]({p}research_note_774.md)在明确N1—N3下，同一局部一圈作用修复可保771来源并共同提升关系观测；原N1—N3共同实现、相互作用正态及仪器仍未签收。三组、十五式通过，最新774／3508，1634份编号科学文件、3815份保护证据。[核验]({p}research_round_774_checks.json)、[条件账]({p}unified_physics_condition_ledger_774.md)。'
order = '**当前执行顺序（774后，优先于下方历史安排）：** 接[775共同规范化的原对象映射]({p}round775_drafts/STATUS.md)，优先核同一实际来源、因果处方及初态/接触边界；不重复本轮条件递推或773同调。[范围审计]({p}round774_drafts/scope_and_dedup_review.md)。原物理准备、旧空间与目标保持。'
planned = {}
for p, raw in before.items():
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl = '\r\n' if b'\r\n' in raw else '\n'
    text = raw.decode(enc).replace('\r\n', '\n')
    assert '**第774轮完成（来源与量子作用的条件性连接）：**' not in text
    prefix = 'research_cognition_physics/archive_231_/' if p.parent == ROOT else 'archive_231_/' if p.parent == RESEARCH else ''
    if p in paths[:5]:
        head, rest = text.split('\n\n', 1)
        text = head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name == 'spatial_premise_closure_audit.md':
        assert '\n## 420.' not in text
        text += '\n\n## 420. 来源与共同量子作用的条件性连接\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 325.' not in text
        text += '\n\n## 325. 旧空间合同保持，条件性一圈共同插入\n\n'+summary.format(p=prefix)+'\n'
    if p.name == 'research_direction.md':
        assert '完成231—773轮。' in text and '最新科学轮次与检查数为773／3505' in text
        text = text.replace('完成231—773轮。', '完成231—774轮。').replace('最新科学轮次与检查数为773／3505', '最新科学轮次与检查数为774／3508')
    if p == HERE/'README.md':
        text += '\n\n## 第774轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|774|[保来源首项与共同插入](research_note_774.md)|[代码](joint_source_action_lift.py)、[结果](joint_source_action_lift_results.json)、[核验](research_round_774_checks.json)|\n'
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
    temp = p.with_name(p.name+'.round774.tmp')
    with temp.open('xb') as f:
        f.write(raw)
    os.replace(temp, p)
report = dict(date='2026-10-04', latest_round=774, cumulative_tests=3508,
              numbered_scientific_files=1634, unique_protected_evidence_files=3815,
              navigation_files=7, navigation_links=links, broken_links=0,
              navigation_snapshots_preserved=True, complete_report_published=True,
              existing_results_preserved=True, next_round=775,
              active_goal_unchanged=True, stage_complete=False, all_checks_passed=True,
              published_navigation_hashes={p.relative_to(ROOT).as_posix(): core.digest(p) for p in paths})
with TARGET.open('x', encoding='utf8', newline='\n') as f:
    f.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
print(json.dumps({k: report[k] for k in ('latest_round', 'navigation_links', 'all_checks_passed')}))
