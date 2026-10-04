"""Publish round 561; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_finite_duration_clock_readout as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round561_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_561_checks.json')
assert checked['round']==561 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round561_20260930'
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

summary='**第561轮完成：** [完整物质源的有限时间记录与资源代价]({p}research_note_561.md)在同一全H、Gauss来源和有限时装置中，证明真实二值记录对源演化时刻有严格正局部响应，并保留整体半有界和源／指针能源账。充分构造成本高，补偿末读、来源和控制仍输入；单响应不等于四坐标。6组、14式及独立终审通过，最新561／2841，995份编号科学文件，1380份保护证据；[核验]({p}research_round_561_checks.json)、[条件账]({p}unified_physics_condition_ledger_561.md)、[合并现状]({p}joint_condition_compression_update_561.md)。后继优先对接物质读口与连续参考，C04／C21、时空与统一仍开放。'
next_steps='\n\n### 第561轮后：按共同对象合并条件，停止单读口常数优化\n\n完整源的有限时动态二值记录已有正面存在证书；560的固定资源无限慢限制仍保留，二者不矛盾。后继562先核548—553与556—561的物质读口、状态及尺度映射，区分径向差与内部取向，不把F直接当空间导数，也不把单时间响应当四坐标。\n\n继承用户“先合并已有证据充分的条件，再攻关键缺口”的顺序。势常数和非最小耦合接引力时须共同匹配；C04自治控制与C21同物质探针仍开放。只有新共同实现、减少独立输入、明确反例或可区分预测再新增轮次。统一目标未改，历史证据保留。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第561轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 207.' not in body
        body+='\n\n## 207. 有限时动态记录与条件合并\n\n'+block+next_steps
    else:
        assert '\n## 112.' not in body
        body+='\n\n## 112. 单一动态记录与四坐标的范围\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[562物质读口与关系参考映射候选](archive_231_/round562_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为560／2835','最新科学轮次与检查数为561／2841')
        body=body.replace('完成231—560轮。','完成231—561轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第561轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|561|[完整物质源的有限时间记录](research_note_561.md)|[代码](joint_finite_duration_clock_readout.py)、'
            '[结果](joint_finite_duration_clock_readout_results.json)、[核验](research_round_561_checks.json)|\n')
        body+='\n**562候选：** [实际物质读口与关系参考映射](round562_drafts/STATUS.md)，尚未完成；先按[合并现状](joint_condition_compression_update_561.md)核对径向／内部取向、连续菜单及共同尺度，不重做单探针常数优化。\n'
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
    temp=path.with_name(path.name+'.round561.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=561,cumulative_tests=2841,numbered_scientific_files=995,
    unique_protected_evidence_files=1380,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=562,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))

