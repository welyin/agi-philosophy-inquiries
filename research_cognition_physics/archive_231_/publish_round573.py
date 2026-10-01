"""Publish round 573; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_gravity_material_coordinates as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round573_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_573_checks.json')
assert checked['round']==573 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round573_20261001'
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

summary='**第573轮完成：** [同源动态几何中的物质参考]({p}research_note_573.md)在572反流后的同一规范物质及几何中，证明局部四参考满秩、h类时；连续屏障和最大值原理承担证明，未另添参考场或改来源。同一菜单在对称壁退化，不能作全局单一坐标。六组、十六式及独立终审通过，最新573／2921，1031份编号科学文件，1495份保护证据；[核验]({p}research_round_573_checks.json)、[条件账]({p}unified_physics_condition_ledger_573.md)。维数、作用与匹配仍输入，实际量子记录、低曲率与统一仍开放。'
next_steps='\n\n### 第573轮后：复用同一经典解，转接量子记录与应力\n\n572同一约束来源与几何现在同时支持短时经典发展、局部物质参考和类时时钟。证明用原来源的规范能量、连续上屏障和严格单调性，不借另一组局部初值，也未额外调整来源。对称壁及镜像读数限制保留，局部坐标不等于实际量子仪器。\n\n已经合并的经典部分直接复用，停止行列式、网格与源参数精扫。574先回查554—566的有限量子读口、来源矩、Gauss及反作用，和568—573同尺度完整H／几何的实际对象差别。优先解决同一来源的实际记录与共同几何应力接口；不能仅凭相同字段名或写出复合算符名称便视为接通。\n\n原有限平直二标量仪器不自动覆盖非平直目标、规范场、动态几何或连续导数平方。若采用有效场论／半经典方法，明确量子—经典映射、来源与尺度前提，不冒称量子引力。低曲率、完整量子连续、手征物质、维数与引力作用来源继续开放，完整目标保持。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第573轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 219.' not in body
        body+='\n\n## 219. 同一约束几何中的局部物质坐标\n\n'+block+next_steps
    else:
        assert '\n## 124.' not in body
        body+='\n\n## 124. 局部坐标辨识并非维数选择\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[574共同物质、量子记录与几何应力候选](archive_231_/round574_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为572／2915','最新科学轮次与检查数为573／2921')
        body=body.replace('完成231—572轮。','完成231—573轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第573轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|573|[同源动态几何中的物质参考](research_note_573.md)|[代码](joint_gravity_material_coordinates.py)、'
            '[结果](joint_gravity_material_coordinates_results.json)、[核验](research_round_573_checks.json)|\n')
        body+='\n**574候选：** [共同物质、量子记录与几何应力](round574_drafts/STATUS.md)，尚未计轮次；复用经典共同模块，优先核有限量子仪器与新来源、应力和尺度的真实映射。\n'
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
    temp=path.with_name(path.name+'.round573.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=573,cumulative_tests=2921,numbered_scientific_files=1031,
    unique_protected_evidence_files=1495,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=574,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
