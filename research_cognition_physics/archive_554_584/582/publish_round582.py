"""Publish round 582 and the consolidated work order without changing the goal."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round582_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_582_checks.json')
assert checked['round']==582 and checked['all_reported_checks_passed']
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round582_20261001'
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

summary='**第582轮完成：** [同一规范质量、量子连接与共同应力]({p}research_note_582.md)将原规范表示、质量Gram和非零规范应力接入同一标量环：质量Gram约束一部分场依赖规范系数，完整应力固定物质—规范交叉项。六组、十四式及独立终审通过，最新582／2970，1058份编号科学文件，1592份保护证据；[核验]({p}research_round_582_checks.json)、[条件账]({p}unified_physics_condition_ledger_582.md)。限所列标量环及一阶EFT表示；非全部量子环、独立算符基、量子在壳来源或连续图匹配。'
order='**执行顺序：** 按[最新合并清单]({p}joint_condition_compression_update_582.md)先收束已有共同实现，分别登记独立输入减少、表示冗余和待接通部分。下一项直接核两种部分量子化的共同变量与度规—物质混合；不继续孤立规范项清点，统一目标不变。'
next_steps='\n\n### 第582轮后：共同量子变量与背景选择\n\n582把原规范质量、标量环连接及领先规范应力放进同一归一，获得可共同保留的系数限制；局部同应力仍不足以决定所算混合系数。该结果只涉及明确标量内部环和局部有效作用，不能直接替代完整量子来源。\n\n依最新合并清单，已兼容的经典来源、局部参考、有限图量子态和记录功能直接复用。共同实现不自动减少公理，有限阶表示压缩也不直接确定全部独立物理参数。\n\n583转向原共同真空的联合Hessian，核固定Jordan与固定Einstein两个标量子块的关系，以及必须保留的度规—物质混合。入口尚未完成，不能只凭经典框架等价或抽象矩阵恒等式宣称完整量子等价。能源—几何尺度、实际量子来源与物理结构选择继续保留。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第582轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+order.format(p=prefix)+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 228.' not in body
        body+='\n\n## 228. 同一规范质量、量子连接与共同应力\n\n'+block+next_steps
    else:
        assert '\n## 133.' not in body
        body+='\n\n## 133. 共同有效项不等于维数或物理结构选择\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[583共同量子变量与冻结背景候选](archive_231_/round583_drafts/STATUS.md)已登记，尚未计完成轮次；[最新条件账](archive_231_/unified_physics_condition_ledger_582.md)。\n'
    if path.name=='research_direction.md':
        assert '最新科学轮次与检查数为581／2964' in body
        body=body.replace('最新科学轮次与检查数为581／2964','最新科学轮次与检查数为582／2970')
        assert '完成231—581轮。' in body
        body=body.replace('完成231—581轮。','完成231—582轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第582轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|582|[同一规范质量、量子连接与共同应力](research_note_582.md)|[代码](joint_gauged_effective_geometry.py)、'
            '[结果](joint_gauged_effective_geometry_results.json)、[核验](research_round_582_checks.json)|\n')
        body+='\n[条件账](unified_physics_condition_ledger_582.md)已更新。**583候选：** [共同量子变量与冻结背景](round583_drafts/STATUS.md)，尚未计完成。\n'
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
support_links=0
support_files=('unified_physics_condition_ledger_582.md','round583_drafts/STATUS.md','joint_condition_compression_update_582.md')
for name in support_files:
    path=HERE/name
    for link in core.link_parser()(path.read_text('utf8')):
        assert (path.parent/link).resolve().exists(),(path,link)
        support_links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body:continue
    temp=path.with_name(path.name+'.round582.tmp')
    if temp.exists():assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=582,cumulative_tests=2970,numbered_scientific_files=1058,
    unique_protected_evidence_files=1592,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=583,
    supporting_document_links=support_links,
    supporting_document_hashes={name:core.digest(HERE/name) for name in support_files},
    consolidation_not_counted_as_scientific_round=True,active_goal_unchanged=True,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
