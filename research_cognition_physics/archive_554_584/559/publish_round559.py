"""Publish round 559; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_finite_time_gauge_probe as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round559_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_559_checks.json')
assert checked['round']==559 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round559_20260930'
resuming=folder.exists()
if resuming:
    manifest=core.read(folder/'manifest.json')
    raws={p:(folder/manifest[str(p.relative_to(ROOT))]['snapshot']).read_bytes() for p in paths}
    for p,raw in raws.items():
        assert hashlib.sha256(raw).hexdigest()==manifest[str(p.relative_to(ROOT))]['sha256']
else:
    raws={p:p.read_bytes() for p in paths}; folder.mkdir(exist_ok=False); manifest={}
    for p,raw in raws.items():
        key=('root_' if p.parent==ROOT else 'research_' if p.parent==RESEARCH else 'archive_')+p.name
        with (folder/key).open('xb') as stream: stream.write(raw)
        manifest[str(p.relative_to(ROOT))]=dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    with (folder/'manifest.json').open('x',encoding='utf8') as stream:
        json.dump(manifest,stream,ensure_ascii=False,indent=2)

summary='**第559轮完成：** [有限时间探针、规范约束与共同能源]({p}research_note_559.md)在完整正势物质源中给不含b的有限时注能预算，并保逐记录Gauss；有限质量指针有充分半有界时长条件，开关功及末端单读口明确入账。固定背景快基态另有实际记录误差界，两种状态范围不混同。6组、14式及独立终审通过，最新559／2829，989份编号科学文件，1360份保护证据；[核验]({p}research_round_559_checks.json)、[条件更新]({p}unified_physics_condition_ledger_559.md)。完整动态物质参考、自治装置、时空与统一仍开放。'
next_steps='\n\n### 第559轮后：把实际记录接到动态物质参考\n\n同一H和F已共同承载Gauss合法性、完整源能量预算和固定背景快读数匹配；其中完整源任意态与固定物质快基态的量词不同，不能合并成未经证明的动态坐标结论。\n\n后继560回查548—549在壳参考和历史有限时记录，检验同一演化／测量窗口保留哪些独立参考量，并与共同几何和应力账核对。优先解决这个跨部门接口，不继续细扫559探针参数；若仅复述时间平均可能丢信息，不计新轮次。\n\n图、规范群、维数、源准备、指针读口和控制来源仍为输入或开放条件。统一目标不改，阶段未完成，历史冻结材料保留。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第559轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 205.' not in body
        body+='\n\n## 205. 有限时间读数、Gauss约束与完整源预算\n\n'+block+next_steps
    else:
        assert '\n## 110.' not in body
        body+='\n\n## 110. 合法有限时探针尚未生成动态关系坐标\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[560动态物质参考的实际记录候选](archive_231_/round560_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为558／2823','最新科学轮次与检查数为559／2829')
        body=body.replace('完成231—558轮。','完成231—559轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第559轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|559|[有限时间探针与共同能源](research_note_559.md)|[代码](joint_finite_time_gauge_probe.py)、'
            '[结果](joint_finite_time_gauge_probe_results.json)、[核验](research_round_559_checks.json)|\n')
        body+='\n**560候选：** [实际记录与动态物质参考](round560_drafts/STATUS.md)，尚未完成；先回查548—549及历史时间平均／可观测性，不重做通用保护测量。\n'
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body: continue
    temp=path.with_name(path.name+'.round559.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=559,cumulative_tests=2829,numbered_scientific_files=989,
    unique_protected_evidence_files=1360,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=560,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))

