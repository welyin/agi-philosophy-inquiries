"""Publish round 552 with existing navigation/history and concurrent-write checks."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_split_mass_geometry as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round552_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_552_checks.json')
assert checked['round']==552 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round552_20260930'
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

summary='**第552轮完成：** [不等质量与共同几何判据]({p}research_note_552.md)证明同一F₀同时控制常场曲率反求与径向Hessian惯性，给正F且正四次势仍为径向鞍点的反例及退化真平坦谷；并把544运行终点接入共同质量与几何匹配。7组、14式及独立终审通过，最新552／2787，968份编号科学文件，1291份保护证据；[核验]({p}research_round_552_checks.json)、[条件更新]({p}unified_physics_condition_ledger_552.md)。终点为精确有理代理，完整曲率RG、量子真空与统一仍开放。'
next_steps='\n\n### 第552轮后：共同判据与弯曲时空一圈接口\n\n同一非最小h＋s常场模型中，F₀将曲率反求、退化与径向正性连接。544运行终点可在声明的引力系数及真空匹配下给平直正面族；物质相区仍须逐分量核查，有限十进制代理不是真实ODE严格误差包络。\n\n后继553先核同一物质对曲率耦合、Einstein系数、真空项和曲率平方反项的共同要求。对接成熟一圈结果，明确量子物质／经典几何范围；不能把树势的辅助尺度变化误称物理真空变化。只有新的交叉约束或共同构造才计新轮次。\n\n统一目标未改，旧轮次和所有冻结证据保留；阶段尚未完成。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第552轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 198.' not in body
        body+='\n\n## 198. 不等质量的曲率与径向共同判据\n\n'+block+next_steps
    else:
        assert '\n## 103.' not in body
        body+='\n\n## 103. 径向稳定与曲率连接不等于维数生成\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[553弯曲时空一圈接口草稿](archive_231_/round553_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为551／2780','最新科学轮次与检查数为552／2787')
        body=body.replace('完成231—551轮。','完成231—552轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第552轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|552|[不等质量与共同几何判据](research_note_552.md)|[代码](joint_split_mass_geometry.py)、'
            '[结果](joint_split_mass_geometry_results.json)、[核验](research_round_552_checks.json)|\n')
        body+='\n**553候选：** [弯曲时空一圈接口](round553_drafts/STATUS.md)，尚未完成；成熟反项映射不另计轮次。\n'
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
    temp=path.with_name(path.name+'.round552.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=552,cumulative_tests=2787,numbered_scientific_files=968,
    unique_protected_evidence_files=1291,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=553,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
