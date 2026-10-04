"""Publish742 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round742_navigation_checks.json'
checked=core.read(HERE/'research_round_742_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round742_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第742轮完成：** [原记录与二次实现的严格边界]({p}research_note_742.md)同一纯参考的真实记录有非高斯四点关联，确定二次费米演化加独立高斯探针的误差至少2p(1−p)/7；633连续例为1/14。仅排除该实现类，保741平均来源及完整量子相互作用候选。一组、十二式通过，最新742／3422，1538份编号科学文件、3387份保护证据。[核验]({p}research_round_742_checks.json)、[条件账]({p}unified_physics_condition_ledger_742.md)。'
order='**当前执行顺序（742后，优先于下方历史安排）：** 接[743原量子相互作用与同一记录]({p}round743_drafts/STATUS.md)，先查原Yukawa／Majorana耦合、实际读口和完整来源；不因高斯限制就添加新物种。旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第742轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 388.' not in text
        text+='\n\n## 388. 原实际记录与二次实现的严格边界\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 293.' not in text
        text+='\n\n## 293. 旧空间合同保持，平均来源不替代完整记录\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—741轮。' in text and '最新科学轮次与检查数为741／3421' in text
        text=text.replace('完成231—741轮。','完成231—742轮。').replace('最新科学轮次与检查数为741／3421','最新科学轮次与检查数为742／3422')
    if p==HERE/'README.md':
        text+='\n\n## 第742轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|742|[原实际记录与二次实现的严格边界](research_note_742.md)|[代码](joint_record_gaussian_boundary.py)、[结果](joint_record_gaussian_boundary_results.json)、[核验](research_round_742_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round742.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=742,cumulative_tests=3422,numbered_scientific_files=1538,
    unique_protected_evidence_files=3387,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=743,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
