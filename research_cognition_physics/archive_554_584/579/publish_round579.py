"""Publish round 579 and the consolidated work order without changing the goal."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round579_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_579_checks.json')
assert checked['round']==579 and checked['all_reported_checks_passed']
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round579_20261001'
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

summary='**第579轮完成：** [单一谱底减除后的残余量子能源]({p}research_note_579.md)证明保留原完整正势时，ξ=1/5后仍有任意Gauss／纠缠态和不均匀节点体积适用的显式能源下界B_N，至少按N²/log²N增长。固定ℏ、有界体积的无限细化仍不能直接接入所列有界常规几何。五组、十二式及独立终审通过，最新579／2952，1049份编号科学文件，1554份保护证据；[核验]({p}research_round_579_checks.json)、[条件账]({p}unified_physics_condition_ledger_579.md)。此为限定接法的障碍，不是实际基态渐近或一般重整化理论的反证。'
order='**执行顺序：** [先合并已兼容条件，再处理剩余连接]({p}joint_condition_compression_update_579.md)。原物质已有限定记录功能，不再反复细化读口；579排除单一谱底减除后，接续同一量子处方、有效作用与几何来源。580已有入口计算和预审，尚未完成。'
next_steps='\n\n### 第579轮后：共同量子处方、有效作用与几何来源\n\n原正势使单一自由谱底减除仍不足；全态有限N下界及其体积条件已在579证明。保留有限截止、完整反项／有效作用匹配和其它受控宏观近似，不继续调整单一能源常数。\n\n先复用已兼容的规范、经典初值、局部参考与固定图记录。580入口核同一H⁵目标在固定Einstein背景的一圈标量有效项，并比较553固定Jordan背景的冻结变量差别；已有四组计算和独立预审，尚未计完整轮次。连续一圈局部结果不能直接当原图非微扰细化的重整化证明。\n\n优先接通此前不同量子处方和尺度的对象，不新增容易但无助共同模型的问题。维数、Lorentz结构、手征物质和引力作用来源继续联立追踪，统一目标不变。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第579轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+order.format(p=prefix)+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 225.' not in body
        body+='\n\n## 225. 原完整势与单一谱底减除的限定障碍\n\n'+block+next_steps
    else:
        assert '\n## 130.' not in body
        body+='\n\n## 130. 残余能源限制与维数问题的范围\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[580共同量子处方与有效作用候选](archive_231_/round580_drafts/STATUS.md)已登记，尚未计完成轮次；[最新条件账](archive_231_/unified_physics_condition_ledger_579.md)。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为578／2947','最新科学轮次与检查数为579／2952')
        body=body.replace('完成231—578轮。','完成231—579轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第579轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|579|[单一谱底减除后的残余量子能源](research_note_579.md)|[代码](joint_residual_quantum_energy.py)、'
            '[结果](joint_residual_quantum_energy_results.json)、[核验](research_round_579_checks.json)|\n')
        body+='\n[条件账](unified_physics_condition_ledger_579.md)已更新。**580候选：** [共同量子处方与有效作用](round580_drafts/STATUS.md)，尚未计完成。\n'
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
support_links=0
support_files=('unified_physics_condition_ledger_579.md','round580_drafts/STATUS.md','joint_condition_compression_update_579.md','round580_drafts/entry_progress.md')
for name in support_files:
    path=HERE/name
    for link in core.link_parser()(path.read_text('utf8')):
        assert (path.parent/link).resolve().exists(),(path,link)
        support_links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body:continue
    temp=path.with_name(path.name+'.round579.tmp')
    if temp.exists():assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=579,cumulative_tests=2952,numbered_scientific_files=1049,
    unique_protected_evidence_files=1554,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=580,
    supporting_document_links=support_links,
    supporting_document_hashes={name:core.digest(HERE/name) for name in support_files},
    consolidation_not_counted_as_scientific_round=True,active_goal_unchanged=True,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
