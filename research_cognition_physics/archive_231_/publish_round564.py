"""Publish round 564; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_radial_prediction_closure as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round564_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_564_checks.json')
assert checked['round']==564 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round564_20260930'
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

summary='**第564轮完成：** [同能量径向资料的预测障碍与精确联合约化]({p}research_note_564.md)在同一全H中构造同径向／中性完整约化态、同平均总能量但未来实际径向记录不同的Gauss来源；精确约化明确需保留联合取向资料，有限时区分另有逐来源存在论证。五组、十四式及独立终审通过，最新564／2858，1004份编号科学文件，1408份保护证据；[核验]({p}research_round_564_checks.json)、[条件账]({p}unified_physics_condition_ledger_564.md)。只限制所列径向预测合同，不否定受控有效理论或认定空间、GR及统一完成。'
next_steps='\n\n### 第564轮后：检验保留物质动力学的受控取向约化\n\n径向实际读取和平均总能源已可共同交代，但还不能删除影响未来的联合取向。完整Gauss变量及准确H已经明确；不重复四参数读口、单探针常数或内部维数类比。\n\n下一项在固定a、k、W而b增大的明确分支中，检验全物质动力学及同一R仪器的共同极限。558的固定端点谱仅作前置；要另核闭形式、实时间传播及仪器。当舍去取向时，还须检查真正的邻接物质相互作用是否保留，不能只看低能子空间更小。\n\n新尺度或匹配条件须单列。连续背景、C04／C21、维数与完整几何反馈仍开放，统一目标不改。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第564轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 210.' not in body
        body+='\n\n## 210. 径向预测与联合内部资料\n\n'+block+next_steps
    else:
        assert '\n## 115.' not in body
        body+='\n\n## 115. 径向记录与独立动力学的区别\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[565全物质快链路极限候选](archive_231_/round565_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为563／2853','最新科学轮次与检查数为564／2858')
        body=body.replace('完成231—563轮。','完成231—564轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第564轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|564|[同能径向预测与联合约化](research_note_564.md)|[代码](joint_radial_prediction_closure.py)、'
            '[结果](joint_radial_prediction_closure_results.json)、[核验](research_round_564_checks.json)|\n')
        body+='\n**565候选：** [完整物质的快链路极限](round565_drafts/STATUS.md)，尚未完成；先核同一无界H及实际R仪器的受控极限，再问邻接相互作用是否保留，不用固定端点谱替代全动力学。\n'
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
    temp=path.with_name(path.name+'.round564.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=564,cumulative_tests=2858,numbered_scientific_files=1004,
    unique_protected_evidence_files=1408,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=565,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))

