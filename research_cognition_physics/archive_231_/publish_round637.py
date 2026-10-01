"""Publish verified 637, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round637_navigation_checks.json'
checked=core.read(HERE/'research_round_637_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round637_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第637轮完成：** [原完整Gauss热态、区域熵与同尺度近似的共同能源条件]({p}research_note_637.md)'
         '原完整H能源下界经同一切分控制区域热迹与有限熵；'
         '原有界注能记录及保能源的同图全态近似共同覆盖。'
         '两组、十四式通过，最新637／3151，1223份编号科学文件、2059份保护证据。'
         '[核验]({p}research_round_637_checks.json)、[条件账]({p}unified_physics_condition_ledger_637.md)。'
         '限固定图与正背景；空间连续、面积律、Newton匹配及GR仍开放。')
order=('**当前执行顺序（637后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[638跨尺度区域态、规范边界与几何熵]({p}round638_drafts/STATUS.md)，'
       '联查旧细化障碍、原区域态及635局部几何项，核共同有限量；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第637轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 283.' not in text
        text+='\n\n## 283. 原完整能源与区域熵的共同成立\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 188.' not in text
        text+='\n\n## 188. 固定图有限区域熵不等于连续面积律\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—636轮。' in text and '最新科学轮次与检查数为636／3149' in text
        text=text.replace('完成231—636轮。','完成231—637轮。').replace('最新科学轮次与检查数为636／3149','最新科学轮次与检查数为637／3151')
    if p==HERE/'README.md':
        text+='\n\n## 第637轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|637|[原完整Gauss热态、区域熵与同尺度近似的共同能源条件](research_note_637.md)|[代码](joint_full_thermal_region_entropy.py)、[结果](joint_full_thermal_region_entropy_results.json)、[核验](research_round_637_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round637.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=637,cumulative_tests=3151,numbered_scientific_files=1223,
            unique_protected_evidence_files=2059,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=638,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
