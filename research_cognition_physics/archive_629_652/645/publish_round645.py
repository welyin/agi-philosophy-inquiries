"""Publish verified 645, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round645_navigation_checks.json'
checked=core.read(HERE/'research_round_645_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round645_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第645轮完成：** [共同热参考、原质量与几何背景的联合条件]({p}research_note_645.md)'
         '原全代质量的同一热函数共同固定标量力、应力和几何接触项；'
         '固定β的lapse必须输送局部温度，仅调真空势不能保持平直常场热背景。'
         '三组、十六式通过，最新645／3174，1247份编号科学文件、2119份保护证据。'
         '[核验]({p}research_round_645_checks.json)、[条件账]({p}unified_physics_condition_ledger_645.md)。'
         '限指定连续自由费米分支；图态映射、完整熵与GR仍未完成。')
order=('**当前执行顺序（645后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[646跨分支共同区域、尺度与几何资料]({p}round646_drafts/STATUS.md)，'
       '核实际共同映射；不继续热势或装置优化，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第645轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 291.' not in text
        text+='\n\n## 291. 共同热参考的质量与几何背景条件\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 196.' not in text
        text+='\n\n## 196. 有限热参考仍须匹配真实几何来源\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—644轮。' in text and '最新科学轮次与检查数为644／3171' in text
        text=text.replace('完成231—644轮。','完成231—645轮。').replace('最新科学轮次与检查数为644／3171','最新科学轮次与检查数为645／3174')
    if p==HERE/'README.md':
        text+='\n\n## 第645轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|645|[共同热参考、原质量与几何背景的联合条件](research_note_645.md)|[代码](joint_thermal_background_source.py)、[结果](joint_thermal_background_source_results.json)、[核验](research_round_645_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round645.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=645,cumulative_tests=3174,numbered_scientific_files=1247,
            unique_protected_evidence_files=2119,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=646,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
