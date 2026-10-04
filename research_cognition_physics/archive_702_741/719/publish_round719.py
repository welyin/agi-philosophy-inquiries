"""Publish719 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round719_navigation_checks.json'
checked=core.read(HERE/'research_round_719_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round719_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第719轮完成：** [局部关系记录、内部参考与原质量反馈]({p}research_note_719.md)原sterile双自旋可承担关系读取及共享参考；能重建均值却不保原Luders历史。读取消去指定参考相干，Dirac／跳跃给反馈，Majorana在等待中生成关联。三组、二十式通过，最新719／3362，1469份编号科学文件、3074份保护证据。[核验]({p}research_round_719_checks.json)、[全条件账]({p}unified_physics_condition_ledger_719.md)。连续因果装置及共同尺度仍开放。'
order='**当前执行顺序（719后，优先于下方历史安排）：** 接[720联合自旋参考与空间传播]({p}round720_drafts/STATUS.md)，先回查604／614／633／666—667等，核同一参考结构与实际Weyl方向传播的兼容性。工程设计后置，停止一般局部相位和参考精度扫描；旧空间及699范围保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第719轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 365.' not in text
        text+='\n\n## 365. 局部关系记录、原物质参考及质量反馈\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 270.' not in text
        text+='\n\n## 270. 旧空间合同复用，关系统计不替代原仪器\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—718轮。' in text and '最新科学轮次与检查数为718／3359' in text
        text=text.replace('完成231—718轮。','完成231—719轮。').replace('最新科学轮次与检查数为718／3359','最新科学轮次与检查数为719／3362')
    if p==HERE/'README.md':
        text+='\n\n## 第719轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|719|[局部关系记录、内部参考与原质量反馈](research_note_719.md)|[代码](joint_local_relational_record.py)、[结果](joint_local_relational_record_results.json)、[核验](research_round_719_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round719.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=719,cumulative_tests=3362,numbered_scientific_files=1469,
    unique_protected_evidence_files=3074,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=720,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
