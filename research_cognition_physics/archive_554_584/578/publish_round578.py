"""Publish round 578 and the consolidated work order without changing the goal."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round578_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_578_checks.json')
assert checked['round']==578 and checked['all_reported_checks_passed']
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round578_20261001'
resuming=folder.exists()
if resuming:
    manifest=core.read(folder/'manifest.json')
    raws={p:(folder/manifest[str(p.relative_to(ROOT))]['snapshot']).read_bytes() for p in paths}
    for p,raw in raws.items():
        assert hashlib.sha256(raw).hexdigest()==manifest[str(p.relative_to(ROOT))]['sha256']
else:
    raws={p:p.read_bytes() for p in paths};folder.mkdir(exist_ok=False);manifest={}
    for p,raw in raws.items():
        key=('root_' if p.parent==ROOT else 'research_' if p.parent==RESEARCH else 'archive_')+p.name
        with (folder/key).open('xb') as stream:stream.write(raw)
        manifest[str(p.relative_to(ROOT))]=dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    with (folder/'manifest.json').open('x',encoding='utf8') as stream:
        json.dump(manifest,stream,ensure_ascii=False,indent=2)

summary='**第578轮完成：** [曲目标量子能源与共同几何尺度]({p}research_note_578.md)证明574原H对任意Gauss／纠缠态有E≥ℏ²N²/(3V)。固定ℏ、有限体积的无限细化，不能将该未减除能源直接匹配到有界曲率和外曲率的Einstein初值。固定背景谱底平移保记录，但几何变分不同；ξ=1/5只是待补的排序分支。四组、十四式及独立终审通过，最新578／2947，1046份编号科学文件，1543份保护证据；[核验]({p}research_round_578_checks.json)、[条件账]({p}unified_physics_condition_ledger_578.md)。不是完整基态、一般量子场论或统一计划的否定。'
order='**执行顺序：** 接续同一物质、能源与几何的共同尺度；[578条件账]({p}unified_physics_condition_ledger_578.md)已排除未经减除的直接细化。579检验减谱底后原正势的残余界，随后比较有限截止及完整有效作用／反项，不继续调整读口。'
next_steps='\n\n### 第578轮后：从自由谱底转向完整残余能源\n\n原完整量子H的全态下界将固定ℏ、有限体积、无限细化和常规几何直接匹配的冲突写成明确不等式。该反证只限声明的未经重整化处方；没有排除有限截止或其它量子处方。\n\nξ=1/5内部目标曲率排序恰可减掉该自由谱底，保持固定图记录，但不能保证完整正势和涨落的能源在细化时有界，也不能省去几何变分。579先检验原正势是否仍迫使残余能源增长；已有Hardy形式、指数势界及不均匀体积分配的入口证据，尚待完整独立核验。\n\n若单一减除仍不够，转向共同有效作用及应力匹配或明确有限截止。保持维数、手征物质、引力来源与观测相容的统一目标，不将一个细化分支失败扩成全计划失败。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第578轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+order.format(p=prefix)+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 224.' not in body
        body+='\n\n## 224. 曲目标能源与常规几何匹配的限定障碍\n\n'+block+next_steps
    else:
        assert '\n## 129.' not in body
        body+='\n\n## 129. 能源连续尺度不等于维数选择\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[579减谱底后完整势候选](archive_231_/round579_drafts/STATUS.md)已登记，尚未计完成轮次；[最新条件账](archive_231_/unified_physics_condition_ledger_578.md)。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为577／2943','最新科学轮次与检查数为578／2947')
        body=body.replace('完成231—577轮。','完成231—578轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第578轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|578|[曲目标量子能源与共同几何尺度](research_note_578.md)|[代码](joint_quantum_energy_geometry.py)、'
            '[结果](joint_quantum_energy_geometry_results.json)、[核验](research_round_578_checks.json)|\n')
        body+='\n[条件账](unified_physics_condition_ledger_578.md)已更新。**579候选：** [减谱底后的原完整势](round579_drafts/STATUS.md)，尚未计完成。\n'
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
support_links=0
support_files=('unified_physics_condition_ledger_578.md','round579_drafts/STATUS.md')
for name in support_files:
    path=HERE/name
    for link in core.link_parser()(path.read_text('utf8')):
        assert (path.parent/link).resolve().exists(),(path,link)
        support_links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body:continue
    temp=path.with_name(path.name+'.round578.tmp')
    if temp.exists():assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=578,cumulative_tests=2947,numbered_scientific_files=1046,
    unique_protected_evidence_files=1543,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=579,
    supporting_document_links=support_links,
    supporting_document_hashes={name:core.digest(HERE/name) for name in support_files},
    consolidation_not_counted_as_scientific_round=True,active_goal_unchanged=True,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
