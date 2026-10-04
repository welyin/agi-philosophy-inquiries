"""Publish760 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round760_navigation_checks.json'
checked=core.read(HERE/'research_round_760_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round760_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第760轮完成：** [混合Gauss过程与条件来源流]({p}research_note_760.md)精确条件振幅保原动态/记录/来源，无需能带隙；原有限能源两态同密度、条件矩阵及总流却异来源流。三组、十六式通过，最新760／3467，1592份编号科学文件、3627份保护证据。[核验]({p}research_round_760_checks.json)、[条件账]({p}unified_physics_condition_ledger_760.md)。这是精确表示，尚非受控宏观或量子引力。'
order='**当前执行顺序（760后，优先于下方历史安排）：** 接[761同一任务与实际压缩]({p}round761_drafts/STATUS.md)，复用旧加权过程估计，检验原记录与来源任务下真正可删除的资料；不把完整相干重写算作宏观化。[范围审计]({p}round760_drafts/scope_and_dedup_review.md)。认知共同候选、旧空间、604、649／699及目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第760轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 406.' not in text
        text+='\n\n## 406. 混合Gauss过程与条件来源流\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 311.' not in text
        text+='\n\n## 311. 旧空间合同保持，条件相干与来源\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—759轮。' in text and '最新科学轮次与检查数为759／3464' in text
        text=text.replace('完成231—759轮。','完成231—760轮。').replace('最新科学轮次与检查数为759／3464','最新科学轮次与检查数为760／3467')
    if p==HERE/'README.md':
        text+='\n\n## 第760轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|760|[混合Gauss过程与条件来源流](research_note_760.md)|[代码](joint_conditional_background.py)、[结果](joint_conditional_background_results.json)、[核验](research_round_760_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round760.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=760,cumulative_tests=3467,numbered_scientific_files=1592,
    unique_protected_evidence_files=3627,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=761,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
