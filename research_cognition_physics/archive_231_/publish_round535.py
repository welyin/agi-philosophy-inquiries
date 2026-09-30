"""Publish round 535 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import fixed_volume_scalar_stability as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round535_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_535_checks.json')
assert checked['round']==535 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round535_20260930'
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
summary='**第535轮完成：** [固定体积与物质稳定性]({p}research_note_535.md)在同一单代含Majorana模型中证明：几何体积约束禁止半径变化，却保留完整裸谱作用的Higgs下降方向，不能单独产生有限常场局部最小。大半径下几何截断准确仍不保证标量截断驻点有效。8组与独立终审通过，最新535／2663，917份编号科学文件，1157份保护证据；[核验]({p}research_round_535_checks.json)、[条件更新]({p}unified_physics_condition_ledger_535.md)。排除限于所列完整裸作用分支；Wilsonian匹配、量子有效势及整体目标仍开放。'
next_steps='\n\n### 第535轮后：共同规范、Yukawa与Higgs边界\n\n535已关闭仅加几何体积约束的单代完整裸谱常场修补。严格f′<0、固定参数及无额外作用是必要范围；总体积、点态体积形式、有限代数unimodularity和Jacobson空间球体积不得混为同一条件。给定V₀的圆球r已定，乘子随方程求出，但V₀来源和观测有效Λ仍未解释。\n\n下一正式编号536按成熟Wilsonian解释检查同一正权下的规范归一化、物理Yukawa和Higgs四次耦合联合边界。先明确规范化和矩约束，再讨论来源清楚且阈值一致的尺度运行；不能把有限Dirac矩阵中的原始参数直接当实验Yukawa。531规范匹配固定的同一r与尺度须同时约束物质部门，不重新各自调权。\n\n不继续扫描已排除裸球面分支，不把体积乘子当自动消除标量势或解决Λ；保留三维、Lorentz、实际记录、内部资源与量子引力等开放接口。及时形成编号报告与条件账；不改应用目标、不设定时任务、不做图像检验。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第535轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 181.' not in body
        body+='\n\n## 181. 体积约束与完整物质变化\n\n'+block+next_steps
    else:
        assert '\n## 86.' not in body
        body+='\n\n## 86. 固定体积未提供有限标量稳定态\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为534／2655','最新科学轮次与检查数为535／2663')
        body=body.replace('完成231—534轮。','完成231—535轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第535轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|535|[固定体积与物质稳定性](research_note_535.md)|[代码](fixed_volume_scalar_stability.py)、'
            '[结果](fixed_volume_scalar_stability_results.json)、[核验](research_round_535_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round535.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=535,cumulative_tests=2663,numbered_scientific_files=917,
    unique_protected_evidence_files=1157,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=536,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
