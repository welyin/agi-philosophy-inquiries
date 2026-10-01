"""Publish round 583 and the consolidated work order without changing the goal."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round583_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_583_checks.json')
assert checked['round']==583 and checked['all_reported_checks_passed']
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round583_20261001'
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

summary='**第583轮完成：** [共同真空中的联合涨落与冻结背景]({p}research_note_583.md)给原五场与共形度规的联合Hessian：固定Jordan和固定Einstein的标量子块即使在同一在壳真空也不同，标准化后径向谱仍有差。保留混合及拉回冻结条件可接通原作用，非线性积分核验通过。五组、十四式及独立终审通过，最新583／2975，1061份编号科学文件，1601份保护证据；[核验]({p}research_round_583_checks.json)、[条件账]({p}unified_physics_condition_ledger_583.md)。限受限切向与作用映射，未完成完整规范固定、量子测度、引力行列式或连续匹配。'
order='**执行顺序：** 沿[合并清单]({p}joint_condition_compression_update_582.md)核共同对象。583定位两个部分量子化不能直接拼接的混合与背景条件；584检验原有限图的测度、演化算符、Gauss与记录能否同时转换，不继续扫描两个质量参数。统一目标不变。'
next_steps='\n\n### 第583轮后：同一有限图量子处方的实际映射\n\n共同真空并不使两个固定度规的标量算符相同；原作用已给精确冻结条件拉回，完整度规—物质的所列受限二次块保持合同。该结果不替代完整配置空间测度、规范固定或引力环。\n\n584保留574原Laplace–Beltrami量子化，核曲测度到平直坐标测度的酉变换、必须保留的量子补项及577原记录。保留作用的逆度量和全部相互作用，不能把换测度误当平直目标或自由删除排序项。\n\n固定几何、有限图和原排序继续作为输入。若确能保全部算符、态与来源，再连接几何变分；578—579的能源障碍、维数和物理结构选择继续开放。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第583轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+order.format(p=prefix)+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 229.' not in body
        body+='\n\n## 229. 共同真空中的联合涨落与冻结背景\n\n'+block+next_steps
    else:
        assert '\n## 134.' not in body
        body+='\n\n## 134. 共同涨落映射不等于量子引力或维数选择\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[584量子测度与有限图演化候选](archive_231_/round584_drafts/STATUS.md)已登记，尚未计完成轮次；[最新条件账](archive_231_/unified_physics_condition_ledger_583.md)。\n'
    if path.name=='research_direction.md':
        assert '最新科学轮次与检查数为582／2970' in body
        body=body.replace('最新科学轮次与检查数为582／2970','最新科学轮次与检查数为583／2975')
        assert '完成231—582轮。' in body
        body=body.replace('完成231—582轮。','完成231—583轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第583轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|583|[共同真空中的联合涨落与冻结背景](research_note_583.md)|[代码](joint_frame_hessian_matching.py)、'
            '[结果](joint_frame_hessian_matching_results.json)、[核验](research_round_583_checks.json)|\n')
        body+='\n[条件账](unified_physics_condition_ledger_583.md)已更新。**584候选：** [量子测度与有限图演化](round584_drafts/STATUS.md)，尚未计完成。\n'
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
support_links=0
support_files=('unified_physics_condition_ledger_583.md','round584_drafts/STATUS.md','joint_condition_compression_update_582.md')
for name in support_files:
    path=HERE/name
    for link in core.link_parser()(path.read_text('utf8')):
        assert (path.parent/link).resolve().exists(),(path,link)
        support_links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body:continue
    temp=path.with_name(path.name+'.round583.tmp')
    if temp.exists():assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=583,cumulative_tests=2975,numbered_scientific_files=1061,
    unique_protected_evidence_files=1601,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=584,
    supporting_document_links=support_links,
    supporting_document_hashes={name:core.digest(HERE/name) for name in support_files},
    consolidation_not_counted_as_scientific_round=True,active_goal_unchanged=True,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
