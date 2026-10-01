"""Publish round 584 and the consolidated work order without changing the goal."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round584_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_584_checks.json')
assert checked['round']==584 and checked['all_reported_checks_passed']
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round584_20261001'
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

summary='**第584轮完成：** [共同量子测度、原物质记录与几何响应]({p}research_note_584.md)在原有限图LB分支中，将测度、闭域、Gauss、演化记录、能源及背景权偏导接入同一酉字典；删除非恒定补项会产生同源记录差并改变几何响应。四组复算、十四式及主代理末稿核验通过；已收到独立实质审查，独立最终哈希签名因服务额度限制未取得，详见[审查记录]({p}round584_drafts/final_review.txt)。最新584／2979，1064份编号科学文件，1620份保护证据；[核验]({p}research_round_584_checks.json)、[条件账]({p}unified_physics_condition_ledger_584.md)。限原有限图处方，未完成完整引力测度、动态量子来源或连续匹配。'
order='**执行顺序：** 按[合并清单]({p}joint_condition_compression_update_584.md)，先合并已有共同实现，区分独立输入减少与表示冗余，再核同一量子来源的能源、流、压力与几何约束。复用既有结果，不继续孤立扫描排序项；统一目标不变。'
next_steps='\n\n### 第584轮后：先收束共同实现，再核量子来源与几何约束\n\n原有限图的量子表示已经将测度、算符、闭域、Gauss、记录、能源和背景权偏导共同接通。此处减少的是同一量子模型的表示选择，原排序、维数、作用和固定几何尚未被消去。\n\n下一轮先核同一有限尺度量子态的能源、流、压力能否与同一几何约束兼容。复用572—574的来源、350的均值反馈限制及578—579的连续能源障碍；既有反例不另计轮次。候选若采用半经典近似、约化几何或内部时间，分别列明新增输入和完整约束边界，不直接将平均Hamiltonian称作可用引力源。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第584轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+order.format(p=prefix)+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 230.' not in body
        body+='\n\n## 230. 同一量子表示的记录与背景响应\n\n'+block+next_steps
    else:
        assert '\n## 135.' not in body
        body+='\n\n## 135. 量子表示共同转换仍未选择维数与引力\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[585共同量子来源与几何约束候选](archive_231_/round585_drafts/STATUS.md)已登记，尚未计完成轮次；[最新条件账](archive_231_/unified_physics_condition_ledger_584.md)。\n'
    if path.name=='research_direction.md':
        assert '最新科学轮次与检查数为583／2975' in body
        body=body.replace('最新科学轮次与检查数为583／2975','最新科学轮次与检查数为584／2979')
        assert '完成231—583轮。' in body
        body=body.replace('完成231—583轮。','完成231—584轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第584轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|584|[共同量子测度、原物质记录与几何响应](research_note_584.md)|[代码](joint_quantum_measure_records.py)、'
            '[结果](joint_quantum_measure_records_results.json)、[核验](research_round_584_checks.json)|\n')
        body+='\n[条件账](unified_physics_condition_ledger_584.md)已更新。**585候选：** [共同量子来源与几何约束](round585_drafts/STATUS.md)，尚未计完成。\n'
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
support_links=0
support_files=('unified_physics_condition_ledger_584.md','round585_drafts/STATUS.md','joint_condition_compression_update_584.md')
for name in support_files:
    path=HERE/name
    for link in core.link_parser()(path.read_text('utf8')):
        assert (path.parent/link).resolve().exists(),(path,link)
        support_links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body:continue
    temp=path.with_name(path.name+'.round584.tmp')
    if temp.exists():assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=584,cumulative_tests=2979,numbered_scientific_files=1064,
    unique_protected_evidence_files=1620,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=585,
    supporting_document_links=support_links,
    supporting_document_hashes={name:core.digest(HERE/name) for name in support_files},
    consolidation_not_counted_as_scientific_round=True,active_goal_unchanged=True,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
