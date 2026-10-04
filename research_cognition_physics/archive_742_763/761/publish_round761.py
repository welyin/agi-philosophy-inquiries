"""Publish761 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round761_navigation_checks.json'
checked=core.read(HERE/'research_round_761_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round761_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第761轮完成（新增缩放的条件结果）：** [有限时间输运与同源反作用]({p}research_note_761.md)在明示Hb,epsilon＋epsilon B分支中，完整矩阵、首阶背景响应及原字段报告有共同固定时间连接。旧主阶B族保留，未证epsilon=1或连续GR。三组、十五式通过，最新761／3470，1595份编号科学文件、3637份保护证据。[核验]({p}research_round_761_checks.json)、[条件账]({p}unified_physics_condition_ledger_761.md)。'
order='**当前执行顺序（761后，优先于下方历史安排）：** 接[762阶次依据与共同尺度]({p}round762_drafts/STATUS.md)，核相对能源阶次的资源解释及同一态/连续来源连接；不继续局部振子精度。[范围审计]({p}round761_drafts/scope_and_dedup_review.md)。旧空间、604、649／699、原缩放及总目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第761轮完成（新增缩放的条件结果）：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 407.' not in text
        text+='\n\n## 407. 有限时间输运与同源反作用\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 312.' not in text
        text+='\n\n## 312. 旧空间合同保持，新增阶次及同源响应\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—760轮。' in text and '最新科学轮次与检查数为760／3467' in text
        text=text.replace('完成231—760轮。','完成231—761轮。').replace('最新科学轮次与检查数为760／3467','最新科学轮次与检查数为761／3470')
    if p==HERE/'README.md':
        text+='\n\n## 第761轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|761|[有限时间输运与同源反作用](research_note_761.md)|[代码](joint_loop_scale_transport.py)、[结果](joint_loop_scale_transport_results.json)、[核验](research_round_761_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round761.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=761,cumulative_tests=3470,numbered_scientific_files=1595,
    unique_protected_evidence_files=3637,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=762,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
