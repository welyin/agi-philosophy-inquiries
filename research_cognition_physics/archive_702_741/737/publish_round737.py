"""Publish737 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round737_navigation_checks.json'
checked=core.read(HERE/'research_round_737_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round737_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第737轮完成：** [响应零方向的局部补足与实际径向闭合]({p}research_note_737.md)原领先局部作用补足谱零方向；R²有限项须按四阶组织。632匹配真空的共同缩放射线不闭合，须保实际中性混合响应。两组、十四式通过，最新737／3413，1523份编号科学文件、3327份保护证据。[核验]({p}research_round_737_checks.json)、[条件账]({p}unified_physics_condition_ledger_737.md)。完整约束与非线性自洽仍开放。'
order='**当前执行顺序（737后，优先于下方历史安排）：** 接[738实际中性混合质量响应]({p}round738_drafts/STATUS.md)，保632同一匹配背景及两个真实径向，核混合阈值、矩阵对数主部和共同局部项；不能把单一缩放射线当闭合动力学。旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第737轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 383.' not in text
        text+='\n\n## 383. 响应零方向的局部补足与实际径向闭合\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 288.' not in text
        text+='\n\n## 288. 旧空间合同保持，局部补足不替代完整约束\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—736轮。' in text and '最新科学轮次与检查数为736／3411' in text
        text=text.replace('完成231—736轮。','完成231—737轮。').replace('最新科学轮次与检查数为736／3411','最新科学轮次与检查数为737／3413')
    if p==HERE/'README.md':
        text+='\n\n## 第737轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|737|[响应零方向的局部补足与实际径向闭合](research_note_737.md)|[代码](joint_response_null_completion.py)、[结果](joint_response_null_completion_results.json)、[核验](research_round_737_checks.json)|\n'
    planned[p]=text.replace('\n',nl).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder.mkdir(exist_ok=False);manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round737.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=737,cumulative_tests=3413,numbered_scientific_files=1523,
    unique_protected_evidence_files=3327,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=738,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
