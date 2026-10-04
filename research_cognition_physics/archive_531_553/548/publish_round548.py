"""Publish round 548 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_matter_reference_coordinates as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round548_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_548_checks.json')
assert checked['round']==548 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round548_20260930'
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
        with (folder/key).open('xb') as stream:
            stream.write(raw)
        manifest[str(p.relative_to(ROOT))]=dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    with (folder/'manifest.json').open('x',encoding='utf8') as stream:
        json.dump(manifest,stream,ensure_ascii=False,indent=2)
summary='**第548轮完成：** [同一物质模型的关系参考]({p}research_note_548.md)证明仅Higgs＋s场值的规范不变菜单秩≤2；同一势的局部在壳配置则可由场值与导数不变量形成满秩四参考，免再添加四个基本参考场。7组及独立终审通过，最新548／2756，956份编号科学文件，1255份保护证据；[核验]({p}research_round_548_checks.json)、[条件更新]({p}unified_physics_condition_ledger_548.md)、[条件压缩表]({p}joint_condition_compression_table.md)。给定几何、状态准备、实际量子读取与完整引力反馈仍开放；局部构造未选择三维。'
next_steps='\n\n### 第548轮后：已有物质参考与动态几何约束\n\n先合并条件的首项落实为548：同一Higgs＋s场内容可在指定局部经典状态中兼作导数参考，但零阶规范不变菜单不足。成熟关系坐标方法和局部存在定理不另计新发现。全局条件压缩表已保存；547与548均已核验，早前待冻结描述保留为历史。\n\n下一正式编号549优先核同一物质参考的动态几何初始约束，先恢复非最小曲率项的度规变分，不能因R=0就漏掉改进应力，也不另换不相连的最小耦合模型冒充原谱反馈。对接成熟Einstein—标量／标量—张量初值方法，明确正F、有效阶数、边界数据及状态来源；存在性、有效范围和实际测量逐项分开。\n\n物质本身由旧候选提供不等于制备、探针和记录已实现；局部满秩也不生成背景、维数或Einstein作用。统一目标继续，不改目标、不设定时任务。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第548轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 194.' not in body
        body+='\n\n## 194. 已有物质的关系参考与状态条件\n\n'+block+next_steps
    else:
        assert '\n## 99.' not in body
        body+='\n\n## 99. 导数参考不等于维数生成\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[549进行中的来源与曲率接口准备](archive_231_/round549_drafts/STATUS.md)已保存；入口检查不计为完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为547／2749','最新科学轮次与检查数为548／2756')
        body=body.replace('完成231—547轮。','完成231—548轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第548轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|548|[同一物质模型的关系参考](research_note_548.md)|[代码](joint_matter_reference_coordinates.py)、'
            '[结果](joint_matter_reference_coordinates_results.json)、[核验](research_round_548_checks.json)|\n')
        body+='\n**549准备中：** [同一曲率作用与几何约束入口](round549_drafts/STATUS.md)，尚未完成，不纳入上述累计数。\n'
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
    temp=path.with_name(path.name+'.round548.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=548,cumulative_tests=2756,numbered_scientific_files=956,
    unique_protected_evidence_files=1255,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=549,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
