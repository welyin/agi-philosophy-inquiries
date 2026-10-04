"""Publish verified 595, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round595_navigation_checks.json'
checked=core.read(HERE/'research_round_595_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round595_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第595轮完成：** [正能源载体、完整记录与同一自主演化]({p}research_note_595.md)'
         '原联合H可接正电池守恒扩张和401内部时钟，有限预算下保完整记录、未知参考及旧来源；'
         '新增门的几何动量力有界且晚时可控。四组、十四式及主代理核验通过，'
         '最新595／3032，1097份编号科学文件、1746份保护证据。'
         '[核验]({p}research_round_595_checks.json)、[条件账]({p}unified_physics_condition_ledger_595.md)。'
         '新增硬件、全局门与准备仍输入，未实现原物质局部装置或完整GR来源；未取得新的独立代理审查。')
order=('**当前执行顺序（595后，优先于下方历史安排）：** 接续'
       '[596原物质与记录硬件的谱兼容性]({p}round596_drafts/STATUS.md)，'
       '核精确全时嵌入、有限窗口及大系统极限的不同条件，接回原能源和几何来源账；统一目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第595轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 241.' not in text
        text+='\n\n## 241. 正能源记录与自主保存的共同来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 146.' not in text
        text+='\n\n## 146. 能源守恒装置仍须核总几何来源\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—594轮。' in text and '最新科学轮次与检查数为594／3028' in text
        text=text.replace('完成231—594轮。','完成231—595轮。').replace('最新科学轮次与检查数为594／3028','最新科学轮次与检查数为595／3032')
    if p==HERE/'README.md':
        text+='\n\n## 第595轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|595|[正能源载体与同一自主演化](research_note_595.md)|[代码](joint_battery_clock_record.py)、[结果](joint_battery_clock_record_results.json)、[核验](research_round_595_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round595.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=595,cumulative_tests=3032,numbered_scientific_files=1097,
            unique_protected_evidence_files=1746,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=596,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
