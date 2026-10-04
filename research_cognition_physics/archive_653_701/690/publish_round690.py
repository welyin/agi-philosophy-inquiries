"""Publish verified 690, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round690_navigation_checks.json'
checked=core.read(HERE/'research_round_690_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round690_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第690轮完成：** [空间传播、原物理电荷与完整平均]({p}research_note_690.md)'
         '原中性通道给完整平均后的精确电荷配对161−72√5，排除将689同片符号无改动当作同一守恒CAR奇偶；不判定Q0反射负性。'
         '两组、十六式通过，最新690／3287，1382份编号科学文件、2628份保护证据。'
         '[核验]({p}research_round_690_checks.json)、[全条件账]({p}unified_physics_condition_ledger_690.md)。'
         '扩大过程、守恒量重识别及受控尺度路线仍开放；旧空间与全目标不变。')
order=('**当前执行顺序（690后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[691守恒流、时间切口与物理观测]({p}round691_drafts/STATUS.md)，'
       '核真实守恒作用、同一来源与正时间支撑；继承612／646，停止纯时间链重证。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第690轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 336.' not in text
        text+='\n\n## 336. 空间传播后的同一物理电荷合同\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 241.' not in text
        text+='\n\n## 241. 旧空间合同复用，守恒量对象须共同匹配\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—689轮。' in text and '最新科学轮次与检查数为689／3285' in text
        text=text.replace('完成231—689轮。','完成231—690轮。').replace('最新科学轮次与检查数为689／3285','最新科学轮次与检查数为690／3287')
    if p==HERE/'README.md':
        text+='\n\n## 第690轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|690|[空间传播、原物理电荷与完整平均](research_note_690.md)|[代码](joint_spatial_charge_dictionary.py)、[结果](joint_spatial_charge_dictionary_results.json)、[核验](research_round_690_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round690.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=690,cumulative_tests=3287,numbered_scientific_files=1382,
            unique_protected_evidence_files=2628,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=691,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
