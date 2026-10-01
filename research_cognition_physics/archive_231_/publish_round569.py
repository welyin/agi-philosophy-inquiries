"""Publish round 569; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_quotient_gauge_completion as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round569_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_569_checks.json')
assert checked['round']==569 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round569_20261001'
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

summary='**第569轮完成：** [共同Higgs、残余Gauss与商群磁补全]({p}research_note_569.md)将旧荷与同尺度耦合接入共同图模型，保留电磁约束；最大Z6商不允许直接搬入基本磁迹，正权字符补项可给全局良定义且逐面单位类唯一零的磁势。八组、十六式及独立终审通过，最新569／2891，1019份编号科学文件，1457份保护证据；[核验]({p}research_round_569_checks.json)、[条件账]({p}unified_physics_condition_ledger_569.md)、[合并现状]({p}joint_condition_compression_update_569.md)。磁形状、商群和匹配仍输入，不等于单迹跨尺度成立、量子连续或统一完成。'
next_steps='\n\n### 第569轮后：先去重，再接共同约束与几何\n\n共同Higgs、旧群荷、同尺度耦合和局部传播现在有明确兼容构造；有限正Hamiltonian可在所选商群上定义。该构造增加磁函数输入，不把同模型存在改成条件等价或唯一性，旧单迹/RG限制保留。\n\n548—549已解决无规范场子部门中的物质参考与局部动态几何约束相容，不重复研究。570先回查同一物质增加非零规范应力后，Gauss、几何约束、共同尺度及真空项的接口；简单Abelian延拓若只是成熟工具推论，只作审计，不凑轮次。\n\n优先合并能复用的共同条件，再处理高连接价值缺口；不继续扫描delta、精细斜率或单探针预算。量子连续、手征物质、维数、Einstein动力学来源和实际装置仍开放。目标保持，不创建任务或定时研究。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第569轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 215.' not in body
        body+='\n\n## 215. 残余Gauss与全局商群磁补全\n\n'+block+next_steps
    else:
        assert '\n## 120.' not in body
        body+='\n\n## 120. 共同群与传播仍不选择空间维数\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[570共同规范物质与动态几何候选](archive_231_/round570_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为568／2883','最新科学轮次与检查数为569／2891')
        body=body.replace('完成231—568轮。','完成231—569轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第569轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|569|[共同荷与全局商磁作用](research_note_569.md)|[代码](joint_quotient_gauge_completion.py)、'
            '[结果](joint_quotient_gauge_completion_results.json)、[核验](research_round_569_checks.json)|\n')
        body+='\n**570候选：** [共同规范物质与动态几何](round570_drafts/STATUS.md)，尚未计轮次；先回查549已有约束构造，核新增规范应力及共同尺度，避免重复通用协变化。\n'
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
    temp=path.with_name(path.name+'.round569.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=569,cumulative_tests=2891,numbered_scientific_files=1019,
    unique_protected_evidence_files=1457,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=570,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
