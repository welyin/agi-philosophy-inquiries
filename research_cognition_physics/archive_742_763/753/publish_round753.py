"""Publish753 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round753_navigation_checks.json'
checked=core.read(HERE/'research_round_753_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round753_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第753轮完成：** [原参考、可积扰动与共同背景]({p}research_note_753.md)原连续零色背景有不可积的一阶颜色切向；另以原三方向颜色场构造全经典约束、原局部参考且无连续联合稳定子的解族。原作用不改，未签收量子物理态。三组、十八式通过，最新753／3446，1571份编号科学文件、3541份保护证据。[核验]({p}research_round_753_checks.json)、[条件账]({p}unified_physics_condition_ledger_753.md)。'
order='**当前执行顺序（753后，优先于下方历史安排）：** 接[754共同背景与量子来源]({p}round754_drafts/STATUS.md)，复用731右逆和730—741同态来源，核颜色来源、引力约束及参考的共同接入；不继续任意色场扫描。[范围审计]({p}round753_drafts/scope_and_dedup_review.md)。旧空间、604、649／699与目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第753轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 399.' not in text
        text+='\n\n## 399. 原参考、可积扰动与共同背景\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 304.' not in text
        text+='\n\n## 304. 旧空间合同保持，共同约束背景与原参考\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—752轮。' in text and '最新科学轮次与检查数为752／3443' in text
        text=text.replace('完成231—752轮。','完成231—753轮。').replace('最新科学轮次与检查数为752／3443','最新科学轮次与检查数为753／3446')
    if p==HERE/'README.md':
        text+='\n\n## 第753轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|753|[原参考、可积扰动与共同背景](research_note_753.md)|[代码](joint_reference_constraint_strata.py)、[结果](joint_reference_constraint_strata_results.json)、[核验](research_round_753_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round753.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=753,cumulative_tests=3446,numbered_scientific_files=1571,
    unique_protected_evidence_files=3541,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=754,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
