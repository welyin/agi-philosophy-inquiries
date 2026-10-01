"""Publish round 560; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_probe_time_resolution as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round560_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_560_checks.json')
assert checked['round']==560 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round560_20260930'
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

summary='**第560轮完成：** [完整物质源的时间辨识与慢测量预算]({p}research_note_560.md)证明同一完整无界H的慢探针在固定源预算下，对固定时间位移的全部输出响应至多O(∣s∣/τ)，不需基态、谱隙或冻结；固定单次判别优势要求初始源能量至少Ω(τ²)。给严格Gauss非平稳见证，并明确变分数值不是全H动力学。6组、15式及独立终审通过，最新560／2835，992份编号科学文件，1368份保护证据；[核验]({p}research_round_560_checks.json)、[条件更新]({p}unified_physics_condition_ledger_560.md)。只限制固定慢仪器族，有限时动态参考、时空与统一仍开放。'
next_steps='\n\n### 第560轮后：有限时动态记录的共同窗口\n\n完整源能量控制不自动给动态时间信息：固定来源和指针资源时，无限慢极限无法维持固定单次时间分辨。这不是平稳输入的假象，正规Gauss非平稳源同样受界。\n\n后继561回到同一物质模型的有限时直接或时序联合读口，先核实际独立信息，再与548—549共同几何相接。复用525—530、554—560；不重复一般可观测性，也不继续优化保护极限。必要预算不充当充分实现或全路线无解。\n\n图、群、维数、准备、读口和自治控制来源继续列账。统一目标不变，阶段未完成，历史冻结证据保留。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第560轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 206.' not in body
        body+='\n\n## 206. 完整源时间响应与实际参考条件\n\n'+block+next_steps
    else:
        assert '\n## 111.' not in body
        body+='\n\n## 111. 低注能不等于可读动态坐标\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[561有限时动态记录的共同窗口候选](archive_231_/round561_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为559／2829','最新科学轮次与检查数为560／2835')
        body=body.replace('完成231—559轮。','完成231—560轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第560轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|560|[完整物质源的时间辨识与预算](research_note_560.md)|[代码](joint_probe_time_resolution.py)、'
            '[结果](joint_probe_time_resolution_results.json)、[核验](research_round_560_checks.json)|\n')
        body+='\n**561候选：** [有限时动态记录的共同窗口](round561_drafts/STATUS.md)，尚未完成；先回查525—530和554—560，明确同一物质的实际读口与来源，不重做通用可观测性。\n'
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
    temp=path.with_name(path.name+'.round560.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=560,cumulative_tests=2835,numbered_scientific_files=992,
    unique_protected_evidence_files=1368,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=561,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))

