"""Publish round 580 and the consolidated work order without changing the goal."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round580_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_580_checks.json')
assert checked['round']==580 and checked['all_reported_checks_passed']
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round580_20261001'
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

summary='**第580轮完成：** [同一曲目标的量子有效项与几何响应]({p}research_note_580.md)给原H⁵标量一圈局部系数及其度规变分。仅改U、K和纯几何项不能吸收所算离壳行列式；RS在平坦背景取值零而几何响应可非零。固定Jordan与固定Einstein的部分量子化亦须映射。七组、十三式及独立终审通过，最新580／2959，1052份编号科学文件，1566份保护证据；[核验]({p}research_round_580_checks.json)、[条件账]({p}unified_physics_condition_ledger_580.md)。限给定四维Euclidean标量环，非全模型净发散、原图重整化或独立物理算符分类。'
order='**执行顺序：** 依[合并清单]({p}joint_condition_compression_update_579.md)先压缩已兼容条件。580接通有效项与几何响应；581检查动态度规和物质场重定义能否减少独立项，同时保留观察量映射。不将更多热核项直接当新认知公理，统一目标不变。'
next_steps='\n\n### 第580轮后：先核有效项冗余与共同变量映射\n\n在固定变量离壳处方中，原二导数标量算符类不足；新增局部项还参与几何变分。该结论不等于每一项都需要独立物理参数，亦不等于原图连续能源已被修复。\n\n581入口沿同一K和U，在明确输入的动态Euclidean Einstein领先作用中，检查曲率项与物质项的微扰场重定义。已有多项式因式分解、原势链式关系及RS块的三组入口计算，尚待完整变量替换、余项和记录源映射核验；不计完成轮次。\n\n同时保留553与574不同冻结变量的量子处方匹配、非零规范来源、物理有限截止及宏观尺度的共同条件。局部算符压缩不能替代这些真实连接，也不独立选择三维或引力作用。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第580轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+order.format(p=prefix)+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 226.' not in body
        body+='\n\n## 226. 曲目标有效项与共同几何响应\n\n'+block+next_steps
    else:
        assert '\n## 131.' not in body
        body+='\n\n## 131. 目标曲率、量子处方与时空维数的边界\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[581动态几何下的算符压缩候选](archive_231_/round581_drafts/STATUS.md)已登记，尚未计完成轮次；[最新条件账](archive_231_/unified_physics_condition_ledger_580.md)。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为579／2952','最新科学轮次与检查数为580／2959')
        body=body.replace('完成231—579轮。','完成231—580轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第580轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|580|[同一曲目标的量子有效项与几何响应](research_note_580.md)|[代码](joint_scalar_effective_geometry.py)、'
            '[结果](joint_scalar_effective_geometry_results.json)、[核验](research_round_580_checks.json)|\n')
        body+='\n[条件账](unified_physics_condition_ledger_580.md)已更新。**581候选：** [动态几何下的算符压缩](round581_drafts/STATUS.md)，尚未计完成。\n'
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
support_links=0
support_files=('unified_physics_condition_ledger_580.md','round581_drafts/STATUS.md')
for name in support_files:
    path=HERE/name
    for link in core.link_parser()(path.read_text('utf8')):
        assert (path.parent/link).resolve().exists(),(path,link)
        support_links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body:continue
    temp=path.with_name(path.name+'.round580.tmp')
    if temp.exists():assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=580,cumulative_tests=2959,numbered_scientific_files=1052,
    unique_protected_evidence_files=1566,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=581,
    supporting_document_links=support_links,
    supporting_document_hashes={name:core.digest(HERE/name) for name in support_files},
    consolidation_not_counted_as_scientific_round=True,active_goal_unchanged=True,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
