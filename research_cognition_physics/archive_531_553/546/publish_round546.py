"""Publish round 546 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_top_boundary_identifiability as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round546_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_546_checks.json')
assert checked['round']==546 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round546_20260930'
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
summary='**第546轮完成：** [top唯一反求与共同物质条件]({p}research_note_546.md)证明固定规范轨道和尺度下，终点top参数在整个共同一圈族内严格随q0递增；可达目标唯一固定高能起点及其余共同参数。反求例同时核真空与谱，唯一性不保证物理可行。5组与独立终审通过，最新546／2744，950份编号科学文件，1239份保护证据；[核验]({p}research_round_546_checks.json)、[条件更新]({p}unified_physics_condition_ledger_546.md)。人造MS目标不是实验拟合，阈值与统一仍开放。'
next_steps='\n\n### 第546轮后：未知味方向下的弱流条件\n\n546在给定一圈、固定规范轨道及共同边界族中证明top映射严格单调，可达目标唯一反求q0，并固定其余归一化参数。Yukawa有限存在性不保证全部标量轨道或稳定真空，α／Higgs份额也未自动单调。\n\n下一正式编号547把同一谱接到与味方向无关的中性弱流测量，区分重态关闭和可上壳情形、费米常数归一化与可观测宽度比。不把开放重态的总中性宽度直接当不可见宽度；先形成解析必要条件，再判断是否需要严格全族积分或结构修订。\n\n有限匹配、极点、σ来源、共同几何、维数、Lorentz、设备与资源仍开放。保持统一目标，不另立任务或自动化。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第546轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 192.' not in body
        body+='\n\n## 192. 同一top条件固定共同参数\n\n'+block+next_steps
    else:
        assert '\n## 97.' not in body
        body+='\n\n## 97. 反求唯一不等于物理可行\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为545／2739','最新科学轮次与检查数为546／2744')
        body=body.replace('完成231—545轮。','完成231—546轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第546轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|546|[top唯一反求与共同物质条件](research_note_546.md)|[代码](joint_top_boundary_identifiability.py)、'
            '[结果](joint_top_boundary_identifiability_results.json)、[核验](research_round_546_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round546.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=546,cumulative_tests=2744,numbered_scientific_files=950,
    unique_protected_evidence_files=1239,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=547,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
