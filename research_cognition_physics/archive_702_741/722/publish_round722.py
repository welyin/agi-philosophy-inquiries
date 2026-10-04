"""Publish722 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round722_navigation_checks.json'
checked=core.read(HERE/'research_round_722_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round722_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第722轮完成：** [原平方参考的能源障碍与可容许读取]({p}research_note_722.md)指定平方相位读法使正常Gauss有限能源态出原能源域；完整Q含空间差分的有理instrument却保原能源、实际记录及几何一阶总响应。三组、二十二式通过，最新722／3371，1478份编号科学文件、3116份保护证据。[核验]({p}research_round_722_checks.json)、[全条件账]({p}unified_physics_condition_ledger_722.md)。固定图结论，联合坐标、跨尺度及统一目标开放。'
order='**当前执行顺序（722后，优先于下方历史安排）：** 接[723同一完整参考与关系任务]({p}round723_drafts/STATUS.md)，复用554、647—652、704及707—708，核原四参考、共同态、区域和实际记录映射。停止谱函数、尾部及常数优化，工程设计后置；旧空间及699范围保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第722轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 368.' not in text
        text+='\n\n## 368. 原平方参考、读取与能源域\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 273.' not in text
        text+='\n\n## 273. 旧空间合同复用，单参考读取不自动给完整坐标\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—721轮。' in text and '最新科学轮次与检查数为721／3368' in text
        text=text.replace('完成231—721轮。','完成231—722轮。').replace('最新科学轮次与检查数为721／3368','最新科学轮次与检查数为722／3371')
    if p==HERE/'README.md':
        text+='\n\n## 第722轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|722|[原平方参考的能源障碍与可容许读取](research_note_722.md)|[代码](joint_squared_reference_readout.py)、[结果](joint_squared_reference_readout_results.json)、[核验](research_round_722_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round722.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=722,cumulative_tests=3371,numbered_scientific_files=1478,
    unique_protected_evidence_files=3116,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=723,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
