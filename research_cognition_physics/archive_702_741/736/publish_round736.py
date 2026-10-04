"""Publish736 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round736_navigation_checks.json'
checked=core.read(HERE/'research_round_736_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round736_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第736轮完成：** [原完整谱的短时因果逆]({p}research_note_736.md)原定常完整质量核精确分解为对数主部和可积记忆，活跃来源有短时唯一因果逆；有限局部项须满足明确有界合同。两组、十八式通过，最新736／3411，1520份编号科学文件、3318份保护证据。[核验]({p}research_round_736_checks.json)、[条件账]({p}unified_physics_condition_ledger_736.md)。完整约束、变化背景及非线性自洽仍开放。'
order='**当前执行顺序（736后，优先于下方历史安排）：** 接[737零方向、局部项与实际约束]({p}round737_drafts/STATUS.md)，把632共同接触和经典互补块放回七来源方程；不能仅逆六个通道就签收完整反馈。旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第736轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 382.' not in text
        text+='\n\n## 382. 原完整谱的短时因果逆与反馈适用条件\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 287.' not in text
        text+='\n\n## 287. 旧空间合同保持，短时逆不替代完整约束\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—735轮。' in text and '最新科学轮次与检查数为735／3409' in text
        text=text.replace('完成231—735轮。','完成231—736轮。').replace('最新科学轮次与检查数为735／3409','最新科学轮次与检查数为736／3411')
    if p==HERE/'README.md':
        text+='\n\n## 第736轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|736|[原完整谱的短时因果逆与反馈适用条件](research_note_736.md)|[代码](joint_causal_log_inverse.py)、[结果](joint_causal_log_inverse_results.json)、[核验](research_round_736_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round736.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=736,cumulative_tests=3411,numbered_scientific_files=1520,
    unique_protected_evidence_files=3318,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=737,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
