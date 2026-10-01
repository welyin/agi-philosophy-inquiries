"""Publish verified 607, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round607_navigation_checks.json'
checked=core.read(HERE/'research_round_607_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round607_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第607轮完成：** [保完整历史与共同来源的投影动力学补全]({p}research_note_607.md)'
         '在明示改变H的分支中，同一分块形式保Gauss、热态、指定仪器历史及任意参考，'
         '联络、横向能源和几何来源由同一映射固定；原质量也须同步压缩，不能保留旧谱不核。'
         '三组、十四式通过，最新607／3066，1133份编号科学文件、1842份保护证据。'
         '[核验]({p}research_round_607_checks.json)、[条件账]({p}unified_physics_condition_ledger_607.md)。'
         '局域性、区域组合及完整物质／连续／GR仍开放；无新增独立代理审查。')
order=('**当前执行顺序（607后，优先于下方历史安排）：** 继续由共同过程映射合并条件，认知实现后置。'
       '接[608投影过程与区域组合、局域性]({p}round608_drafts/STATUS.md)，'
       '核原投影及补全的规模一致传播／拼接条件，保留质量和测度缺口，目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第607轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 253.' not in text
        text+='\n\n## 253. 投影补全的历史、Gauss与来源共同闭合\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 158.' not in text
        text+='\n\n## 158. 修改动力学后的交织不替代局域性\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—606轮。' in text and '最新科学轮次与检查数为606／3063' in text
        text=text.replace('完成231—606轮。','完成231—607轮。').replace('最新科学轮次与检查数为606／3063','最新科学轮次与检查数为607／3066')
    if p==HERE/'README.md':
        text+='\n\n## 第607轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|607|[保完整历史与共同来源的投影动力学补全](research_note_607.md)|[代码](joint_projected_process_completion.py)、[结果](joint_projected_process_completion_results.json)、[核验](research_round_607_checks.json)|\n'
    planned[p]=text.replace('\n',newline).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder.mkdir(exist_ok=False)
manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as stream:stream.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as stream:json.dump(manifest,stream,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round607.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=607,cumulative_tests=3066,numbered_scientific_files=1133,
            unique_protected_evidence_files=1842,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=608,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
