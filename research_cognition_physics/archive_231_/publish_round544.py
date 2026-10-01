"""Publish round 544 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_singlet_common_mass_rg as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round544_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_544_checks.json')
assert checked['round']==544 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round544_20260930'
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
summary='**第544轮完成：** [共同质量运行与阈值条件]({p}research_note_544.md)联立实σ、保护中微子对和两种质量的一圈运行；解析证明旧平衡点在明确正号条件下进入严格双凝聚，非零门户例越过旧同尺度比值界。新模起初轻，未取得重阈值或全族排除。8组与独立终审通过，最新544／2732，944份编号科学文件，1223份保护证据；[核验]({p}research_round_544_checks.json)、[条件更新]({p}unified_physics_condition_ledger_544.md)。MS共同高边界属新增匹配输入，完整物理拟合与统一仍开放。'
next_steps='\n\n### 第544轮后：同一混合谱与物理阈值\n\n544证明同一高能负质量边界会在共同一圈运行中分裂；旧q=S零σ曲率点进入严格双凝聚，不能将543静态界跨尺度沿用。公共质量重标仍不能改变无量纲层级。所得是运行参数的固定尺度树势，不是完整Coleman–Weinberg真空、实际阈值或观测拟合。\n\n下一正式编号545把同一径向混合矩阵与完整中性费米子谱联立，核Schur补、可观测Higgs成分、活跃轻子权重及重场条件。缺少分离时不得仅用形式λ_eff签收物理Higgs；先给解析关系或反例，不用加密网格代替全族证明。\n\n明确保留σ来源、裸到MS有限匹配、共同几何、维数、Lorentz、实际记录与资源等开放输入。目标保持进行中，不另设任务或自动化。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第544轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 190.' not in body
        body+='\n\n## 190. 共同质量运行与实际谱条件\n\n'+block+next_steps
    else:
        assert '\n## 95.' not in body
        body+='\n\n## 95. 运行解除边界但未完成重阈值\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为543／2724','最新科学轮次与检查数为544／2732')
        body=body.replace('完成231—543轮。','完成231—544轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第544轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|544|[共同质量运行与阈值条件](research_note_544.md)|[代码](joint_singlet_common_mass_rg.py)、'
            '[结果](joint_singlet_common_mass_rg_results.json)、[核验](research_round_544_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round544.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=544,cumulative_tests=2732,numbered_scientific_files=944,
    unique_protected_evidence_files=1223,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=545,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
