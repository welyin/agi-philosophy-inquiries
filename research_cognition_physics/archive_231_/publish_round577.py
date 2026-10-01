"""Publish round 577 and the consolidated work order without changing the goal."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round577_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_577_checks.json')
assert checked['round']==577 and checked['all_reported_checks_passed']
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round577_20261001'
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

summary='**第577轮完成：** [原物质区域的相位传递与记录]({p}research_note_577.md)在574原完整有限图H中，严格Gauss来源的±相位准备同初始读者统计、同正能源成本，经原相互作用产生第三阶有界记录概率差。任意固定ℏ成立，本任务无需新增χ字段与传播连接。三组、十四式及独立终审通过，最新577／2943，1043份编号科学文件，1534份保护证据；[核验]({p}research_round_577_checks.json)、[条件账]({p}unified_physics_condition_ledger_577.md)。准备、终读及固定几何仍输入，非普适装置或量子几何，统一仍开放。'
order='**执行顺序：** 原物质承担限定接收功能，已减少两类额外装置输入；[577条件账]({p}unified_physics_condition_ledger_577.md)保留全部前提。接续量子能源、共同几何与连续尺度，先核原曲目标谱下界；不继续优化单读口灵敏度。'
next_steps='\n\n### 第577轮后：原物质记录接回量子能源与共同几何\n\n原有限图H已能将相邻端的相位信息传给有界读口，无需新χ或gs dT连接。该结果限明确初态、准备和末读；不是普适量子装置或稳定记忆。初始完整读者相同与同平均总能量也不保证相同应力或几何。\n\n下一步先核574完整双曲目标的量子动能下界，与w=epsilon³ psi⁶的图细化权是否允许固定ℏ、同一常规几何和未减除能源同时成立。成熟双曲谱工具直接复用；只有原模型的联合条件得到新限制才计完整轮次。\n\n明确区分固定图半经典极限、固定ℏ连续尺度及物理应力定义。若需要排序或反项，逐项登记，不把固定几何上的全局相位等价当作相同引力应力。不重做350一般均场反例，保持维数、物质、引力来源与统一目标开放。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第577轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+order.format(p=prefix)+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 223.' not in body
        body+='\n\n## 223. 原物质记录与同一量子能源来源\n\n'+block+next_steps
    else:
        assert '\n## 128.' not in body
        body+='\n\n## 128. 原物质转导不等于空间维数或连续几何\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[578量子能源与几何尺度候选](archive_231_/round578_drafts/STATUS.md)已登记，尚未计完成轮次；[最新条件账](archive_231_/unified_physics_condition_ledger_577.md)。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为576／2940','最新科学轮次与检查数为577／2943')
        body=body.replace('完成231—576轮。','完成231—577轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第577轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|577|[原物质区域的相位传递与记录](research_note_577.md)|[代码](joint_original_matter_record.py)、'
            '[结果](joint_original_matter_record_results.json)、[核验](research_round_577_checks.json)|\n')
        body+='\n[条件账](unified_physics_condition_ledger_577.md)已更新。**578候选：** [曲目标量子能源与共同几何](round578_drafts/STATUS.md)，尚未计完成。\n'
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
support_links=0
support_files=('unified_physics_condition_ledger_577.md','round578_drafts/STATUS.md')
for name in support_files:
    path=HERE/name
    for link in core.link_parser()(path.read_text('utf8')):
        assert (path.parent/link).resolve().exists(),(path,link)
        support_links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body:continue
    temp=path.with_name(path.name+'.round577.tmp')
    if temp.exists():assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=577,cumulative_tests=2943,numbered_scientific_files=1043,
    unique_protected_evidence_files=1534,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=578,
    supporting_document_links=support_links,
    supporting_document_hashes={name:core.digest(HERE/name) for name in support_files},
    consolidation_not_counted_as_scientific_round=True,active_goal_unchanged=True,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
