"""Publish round 586, preserving all navigation bytes and the active goal."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
TARGET = HERE/'round586_navigation_checks.json'
assert not TARGET.exists(), 'already published'
checked = core.read(HERE/'research_round_586_checks.json')
assert checked['round'] == 586 and checked['all_reported_checks_passed']
for group in ('new_file_hashes', 'preserved_draft_hashes'):
    for name, digest in checked[group].items():
        assert core.digest(HERE/name) == digest, name
paths = [ROOT/'README.md', RESEARCH/'README.md', RESEARCH/'research_direction.md',
         RESEARCH/'RESEARCH_STATE.md', HERE/'README.md',
         HERE/'spatial_premise_closure_audit.md', HERE/'three_dimensional_four_conditions_review.md']
raws = {p: p.read_bytes() for p in paths}
folder = HERE/'navigation_before_round586_20261001'
assert not folder.exists()
summary = ('**第586轮完成：** [共同记录压缩、后继过程与同源反作用]({p}research_note_586.md)'
           '在原完整曲目标物质中证明：两次原读取只报告奇偶，与直接奇偶读同即时概率，却有不同源注能、原H下后继读数及几何响应。'
           '给本仪器类的无损替换条件和有限历史源预算，复用506／557，不新增物质H。'
           '五组、十四式及主代理核验通过；最新586／2989，1070份编号科学文件、1642份保护证据。'
           '[核验]({p}research_round_586_checks.json)、[条件账]({p}unified_physics_condition_ledger_586.md)。'
           '理想仪器与准备仍输入，未完成自主装置、动态量子几何或连续极限；未取得新的独立代理审查。')
order = ('**当前执行顺序（586后，优先于下方历史安排）：** 依[正向成果回顾]({p}cognition_forward_bridge_review_585.md)，'
         '用完整过程及反作用检验共同映射。接续[587能源流与切向几何候选]({p}round587_drafts/STATUS.md)，'
         '从同一H和局部能源分配推流，纳入曲目标、规范场及原仪器来源；不继续优化读口精度。统一目标不变。')
planned = {}
for p, raw in raws.items():
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline = '\r\n' if b'\r\n' in raw else '\n'
    body = raw.decode(enc).replace('\r\n', '\n')
    assert '**第586轮完成：**' not in body
    prefix = 'research_cognition_physics/archive_231_/' if p.parent == ROOT else 'archive_231_/' if p.parent == RESEARCH else ''
    if p in paths[:5]:
        first, rest = body.split('\n\n', 1)
        body = first+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name == 'spatial_premise_closure_audit.md':
        assert '\n## 232.' not in body
        body += '\n\n## 232. 共同记录须保持后继过程和来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 137.' not in body
        body += '\n\n## 137. 摘要概率并非共同物理过程的充分接口\n\n'+summary.format(p=prefix)+'\n'
    if p.name == 'research_direction.md':
        assert '最新科学轮次与检查数为585／2984' in body and '完成231—585轮。' in body
        body = body.replace('最新科学轮次与检查数为585／2984', '最新科学轮次与检查数为586／2989')
        body = body.replace('完成231—585轮。', '完成231—586轮。')
    if p == HERE/'README.md':
        body += ('\n\n## 第586轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
                 '|586|[共同记录压缩、后继过程与同源反作用](research_note_586.md)|'
                 '[代码](joint_record_source_compression.py)、[结果](joint_record_source_compression_results.json)、'
                 '[核验](research_round_586_checks.json)|\n\n'
                 '[587能源流与切向几何](round587_drafts/STATUS.md)已登记，尚未完成。\n')
    planned[p] = body.replace('\n', newline).encode(enc)
links = 0
for p, raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(), (p, link)
        links += 1
support = ('unified_physics_condition_ledger_586.md', 'round587_drafts/STATUS.md', 'cognition_forward_bridge_review_585.md')
support_links = 0
for name in support:
    for link in core.link_parser()((HERE/name).read_text('utf8')):
        assert (HERE/name).parent.joinpath(link).resolve().exists(), (name, link)
        support_links += 1
assert all(p.read_bytes() == raw for p, raw in raws.items()), 'concurrent edit'
folder.mkdir(exist_ok=False)
manifest = {}
for i, (p, raw) in enumerate(raws.items()):
    name = f'{i}_{p.name}'
    with (folder/name).open('xb') as stream: stream.write(raw)
    manifest[p.relative_to(ROOT).as_posix()] = dict(snapshot=name, sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x', encoding='utf8') as stream:
    json.dump(manifest, stream, ensure_ascii=False, indent=2)
for p, raw in planned.items():
    assert p.read_bytes() == raws[p]
    temp = p.with_name(p.name+'.round586.tmp')
    with temp.open('xb') as stream: stream.write(raw)
    os.replace(temp, p)
report = dict(date='2026-10-01', latest_round=586, cumulative_tests=2989,
              numbered_scientific_files=1070, unique_protected_evidence_files=1642,
              navigation_files=7, navigation_links=links, broken_links=0,
              supporting_document_links=support_links,
              supporting_document_hashes={name: core.digest(HERE/name) for name in support},
              navigation_snapshots_preserved=True, complete_report_published=True,
              existing_results_preserved=True, next_round=587, active_goal_unchanged=True,
              stage_complete=False, all_checks_passed=True)
with TARGET.open('x', encoding='utf8', newline='\n') as stream:
    stream.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
print(json.dumps({k: report[k] for k in ('latest_round', 'navigation_files', 'navigation_links', 'broken_links', 'all_checks_passed')}))
