"""Publish754 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round754_navigation_checks.json'
checked=core.read(HERE/'research_round_754_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round754_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第754轮完成：** [全颜色来源与共同约束]({p}research_note_754.md)原753背景的全模颜色Gauss逆接通任意光滑来源及完整动量/能量补偿；真实Hadamard态差可带非零颜色均值，并接条件首阶反作用。数值为外给诊断，未签收全量子物理态。三组、二十式通过，最新754／3449，1574份编号科学文件、3550份保护证据。[核验]({p}research_round_754_checks.json)、[条件账]({p}unified_physics_condition_ledger_754.md)。'
order='**当前执行顺序（754后，优先于下方历史安排）：** 接[755共同物理状态与过程]({p}round755_drafts/STATUS.md)，按[认知观察与共同候选工作报告]({p}round755_drafts/cognitive_joint_candidate_working_report.md)收拢三条假说、替代分支和失败判据；下一项只核同一量子过程、来源与受约束背景的连接。工作报告不增加科学轮次，不继续诊断源或读口精度支线。[范围审计]({p}round754_drafts/scope_and_dedup_review.md)。旧空间、604、649／699与应用目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第754轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 400.' not in text
        text+='\n\n## 400. 全颜色来源与共同约束\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 305.' not in text
        text+='\n\n## 305. 旧空间合同保持，量子来源与完整初始约束\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—753轮。' in text and '最新科学轮次与检查数为753／3446' in text
        text=text.replace('完成231—753轮。','完成231—754轮。').replace('最新科学轮次与检查数为753／3446','最新科学轮次与检查数为754／3449')
    if p==HERE/'README.md':
        text+='\n\n## 第754轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|754|[全颜色来源与共同约束](research_note_754.md)|[代码](joint_irreducible_source_completion.py)、[结果](joint_irreducible_source_completion_results.json)、[核验](research_round_754_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round754.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=754,cumulative_tests=3449,numbered_scientific_files=1574,
    unique_protected_evidence_files=3550,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=755,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
