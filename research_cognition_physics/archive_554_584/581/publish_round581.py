"""Publish round 581 and the consolidated work order without changing the goal."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round581_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_581_checks.json')
assert checked['round']==581 and checked['all_reported_checks_passed']
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round581_20261001'
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

summary='**第581轮完成：** [共同几何与物质有效项的条件压缩]({p}research_note_581.md)在声明的动态Euclidean Einstein—标量领先作用下，将580完整标量系数分解为约化项、EOM与总导数，并给共同度规／标量一阶重定义。实际原作用替换核二阶余项，长度观察量必须同步拉回。五组、十三式及独立终审通过，最新581／2964，1055份编号科学文件，1580份保护证据；[核验]({p}research_round_581_checks.json)、[条件账]({p}unified_physics_condition_ledger_581.md)。限零规范曲率标量部门及有限EFT阶，非完整独立算符基、全框架量子等价或原图连续证明。'
order='**执行顺序：** 沿[合并清单]({p}joint_condition_compression_update_579.md)压缩共同条件，保留来源与观察量映射。581减少了逐项新增假设的需要；582接回原非零规范背景的目标连接与应力，检查纯标量压缩的适用边界，不继续孤立算符清点。'
next_steps='\n\n### 第581轮后：回到非零规范共同来源\n\n一批局部曲率项可在同一领先Einstein—标量作用中转为物质有效项；实际变量替换保留二阶余项，几何读数须同步转换。该压缩是条件性表示等价，不是引力作用、全部物质或独立物理参数的唯一推导。\n\n582入口接回原U(1)荷3与H⁵目标连接，已有两组局部矩阵计算。规范曲率反号可保持经典规范应力，却改变标量环的梯度—规范交叉系数；该结果仍待完整来源映射与独立核验，不是两份已实现Gauss／Einstein态的比较。\n\n下一步联立非零规范应力和目标连接的有效项，保持原作用与来源一致。旧553／574量子处方匹配、有限尺度和观测继续在主账内；不把纯标量部门的收束当共同模型完成。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第581轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+order.format(p=prefix)+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 227.' not in body
        body+='\n\n## 227. 共同几何与物质有效项的条件压缩\n\n'+block+next_steps
    else:
        assert '\n## 132.' not in body
        body+='\n\n## 132. 表示压缩不等于维数或物质选择\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[582非零规范背景与共同来源候选](archive_231_/round582_drafts/STATUS.md)已登记，尚未计完成轮次；[最新条件账](archive_231_/unified_physics_condition_ledger_581.md)。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为580／2959','最新科学轮次与检查数为581／2964')
        body=body.replace('完成231—580轮。','完成231—581轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第581轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|581|[共同几何与物质有效项的条件压缩](research_note_581.md)|[代码](joint_effective_operator_reduction.py)、'
            '[结果](joint_effective_operator_reduction_results.json)、[核验](research_round_581_checks.json)|\n')
        body+='\n[条件账](unified_physics_condition_ledger_581.md)已更新。**582候选：** [非零规范背景与共同来源](round582_drafts/STATUS.md)，尚未计完成。\n'
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
support_links=0
support_files=('unified_physics_condition_ledger_581.md','round582_drafts/STATUS.md')
for name in support_files:
    path=HERE/name
    for link in core.link_parser()(path.read_text('utf8')):
        assert (path.parent/link).resolve().exists(),(path,link)
        support_links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body:continue
    temp=path.with_name(path.name+'.round581.tmp')
    if temp.exists():assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=581,cumulative_tests=2964,numbered_scientific_files=1055,
    unique_protected_evidence_files=1580,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=582,
    supporting_document_links=support_links,
    supporting_document_hashes={name:core.digest(HERE/name) for name in support_files},
    consolidation_not_counted_as_scientific_round=True,active_goal_unchanged=True,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
