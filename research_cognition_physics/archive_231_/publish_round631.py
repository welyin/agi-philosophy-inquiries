"""Publish verified 631, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round631_navigation_checks.json'
checked=core.read(HERE/'research_round_631_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round631_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第631轮完成：** [原整代物质的张量应力谱与共同曲率系数]({p}research_note_631.md)'
         '原同一连续真空的应力切割具有迹与五个无迹通道；'
         '张量正谱复现旧费米Weyl平方运行系数，质量和几何来源共用同一物质计数。'
         '三组、十四式通过，最新631／3133，1205份编号科学文件、2015份保护证据。'
         '[核验]({p}research_round_631_checks.json)、[条件账]({p}unified_physics_condition_ledger_631.md)。'
         '限平直非接触二点；背景自洽、接触项、图映射及GR仍开放。')
order=('**当前执行顺序（631后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[632共同态、背景方程与接触项]({p}round632_drafts/STATUS.md)，'
       '先核旧真空与同阶变分，再接入实时张量来源；不把诊断背景当自洽几何，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第631轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 277.' not in text
        text+='\n\n## 277. 原物质张量谱与曲率系数的共同连接\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 182.' not in text
        text+='\n\n## 182. 应力自旋2不是引力子或空间生成\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—630轮。' in text and '最新科学轮次与检查数为630／3130' in text
        text=text.replace('完成231—630轮。','完成231—631轮。').replace('最新科学轮次与检查数为630／3130','最新科学轮次与检查数为631／3133')
    if p==HERE/'README.md':
        text+='\n\n## 第631轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|631|[原整代物质的张量应力谱与共同曲率系数](research_note_631.md)|[代码](joint_tensor_stress_spectrum.py)、[结果](joint_tensor_stress_spectrum_results.json)、[核验](research_round_631_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round631.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=631,cumulative_tests=3133,numbered_scientific_files=1205,
            unique_protected_evidence_files=2015,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=632,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
