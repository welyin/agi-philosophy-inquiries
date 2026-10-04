"""Publish round 567; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_scalar_propagation_matching as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round567_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_567_checks.json')
assert checked['round']==567 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round567_20261001'
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

summary='**第567轮完成：** [同一物质的图传播与连续参考]({p}research_note_567.md)核原图同势双标量的经典主部缺项；明示补singlet邻接规则及共同尺度后，保留同一真空和质量谱。门户允许单模有效传播，不能替代两场完整匹配。六组、十四式及独立终审通过，最新567／2875，1013份编号科学文件，1436份保护证据；[核验]({p}research_round_567_checks.json)、[条件账]({p}unified_physics_condition_ledger_567.md)。新增梯度及尺度仍输入，未建立量子连续极限、三维或统一。'
next_steps='\n\n### 第567轮后：合并现有条件，核规范与物质共同尺度\n\n当前h＋s的势、真空与质量矩阵可沿用；两场共同非零主部需要singlet邻接规则及相对速度匹配。原固定kappa的细化速度趋零，非零速度尺度族另列输入。低能单模、量子及引力诱导梯度仍是替代路线；不把具体缺项扩大成全理论失败。\n\n后继568优先核同一给定图的电项、闭环磁项和物质梯度能否共用参数尺度。先回查历史与成熟格点规范理论；原Gauss和群链路不自动产生磁作用或其系数。图、群、维数及新增闭环作用清楚列账。\n\n执行顺序保持：复用已相容条件，补能连接多部门的接口，再攻剩余核心缺口；不继续扫描567低p斜率或优化单探针常数。冻结历史、统一目标和C04/C21及完整几何反馈的开放状态保持。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第567轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 213.' not in body
        body+='\n\n## 213. 同一物质的传播补全与共同尺度\n\n'+block+next_steps
    else:
        assert '\n## 118.' not in body
        body+='\n\n## 118. 给定格图的双标量主部与维数输入\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[568闭环规范作用与共同尺度候选](archive_231_/round568_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为566／2869','最新科学轮次与检查数为567／2875')
        body=body.replace('完成231—566轮。','完成231—567轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第567轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|567|[同势物质的传播补全](research_note_567.md)|[代码](joint_scalar_propagation_matching.py)、'
            '[结果](joint_scalar_propagation_matching_results.json)、[核验](research_round_567_checks.json)|\n')
        body+='\n**568候选：** [闭环规范作用与共同尺度](round568_drafts/STATUS.md)，尚未计算；先核磁作用与原图电项的准确映射，再联查物质尺度，不把成熟公式复述计成新轮次。\n'
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
    temp=path.with_name(path.name+'.round567.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=567,cumulative_tests=2875,numbered_scientific_files=1013,
    unique_protected_evidence_files=1436,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=568,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))

