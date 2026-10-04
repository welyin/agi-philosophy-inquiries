"""Publish verified 612, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round612_navigation_checks.json'
checked=core.read(HERE/'research_round_612_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round612_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第612轮完成：** [完整时间关联、正演化与热态条件的共同匹配]({p}research_note_612.md)'
         '继承自由overlap核的非零动量关联含正连续谱；可正演化重建，'
         '却不能同时精确全部时间匹配原固定图的紧谱与有限热迹。'
         '三组、十四式通过，最新612／3081，1148份编号科学文件、1877份保护证据。'
         '[核验]({p}research_round_612_checks.json)、[条件账]({p}unified_physics_condition_ledger_612.md)。'
         '仅排除明示自由强匹配分支；有限窗口、完整规范及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（612后，优先于下方历史安排）：** 先整合条件，认知系统设计后置。'
       '接[613共同匹配范围与剩余条件整合]({p}round613_drafts/STATUS.md)，'
       '回填598—612结果，比较量子过程、物质、热态和几何来源的同一尺度合同；'
       '不继续向记忆或辅助硬件设计下钻，目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第612轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 258.' not in text
        text+='\n\n## 258. 同一时间谱、正演化与热态合同\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 163.' not in text
        text+='\n\n## 163. 正演化与原热态条件须共同匹配\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—611轮。' in text and '最新科学轮次与检查数为611／3078' in text
        text=text.replace('完成231—611轮。','完成231—612轮。').replace('最新科学轮次与检查数为611／3078','最新科学轮次与检查数为612／3081')
    if p==HERE/'README.md':
        text+='\n\n## 第612轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|612|[完整时间关联、正演化与热态条件的共同匹配](research_note_612.md)|[代码](joint_overlap_time_spectrum.py)、[结果](joint_overlap_time_spectrum_results.json)、[核验](research_round_612_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round612.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=612,cumulative_tests=3081,numbered_scientific_files=1148,
            unique_protected_evidence_files=1877,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=613,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
