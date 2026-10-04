"""Publish round 541 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import protected_pair_finite_matching as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round541_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_541_checks.json')
assert checked['round']==541 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round541_20260930'
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
summary='**第541轮完成：** [真实有限匹配与质量项联动]({p}research_note_541.md)将成熟重场匹配接入共同边界，证明精确保护配对在声明层级与维四下方运行中，任意允许top／阈值仍有λ>0.16259。负质量目标要求联动二次项输入，不能仅调整四次项。7组与独立终审通过，最新541／2709，935份编号科学文件，1201份保护证据；[核验]({p}research_round_541_checks.json)、[条件更新]({p}unified_physics_condition_ledger_541.md)。维六反馈、一般破缺及完整统一仍开放。'
next_steps='\n\n### 第541轮后：检查同一匹配生成的维六反馈\n\n541把真实一圈有限匹配接入540预算；精确保护配对、匹配前质量比绝对值≤0.01、一次阈值及维四下方运行中，有限跃变绝对值<0.005，低能λ仍>0.16259。三例同时对接历史负质量参数，但该质量为独立反向输入，尚未施加高能谱质量边界。\n\n下一正式编号542核C5=0仍存活的树级Higgs—轻子流算符，利用成熟SEFT方程判断其对λ反馈的符号和预算；不能只因轻质量相消就删除维六作用。区分固定一圈贡献、部分RG重求和和完整SMEFT；一般非简并破缺与高能质量接口仍待核，不再扫描已排除维四分支。\n\n条件账从本轮恢复首版C编号：运行与截断=C20，观测=C27，质量=C17—C18。旧短表标签漂移已明示，旧冻结科学结果不改；不据标签签收宇宙学／紫外问题。维数、Lorentz、共同几何、实际记录和内部资源仍开放，应用目标保持不变。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第541轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 187.' not in body
        body+='\n\n## 187. 有限匹配与同一质量条件\n\n'+block+next_steps
    else:
        assert '\n## 92.' not in body
        body+='\n\n## 92. 保护分支有限匹配的限定结论\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为540／2702','最新科学轮次与检查数为541／2709')
        body=body.replace('完成231—540轮。','完成231—541轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第541轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|541|[真实有限匹配与质量项联动](research_note_541.md)|[代码](protected_pair_finite_matching.py)、'
            '[结果](protected_pair_finite_matching_results.json)、[核验](research_round_541_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round541.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=541,cumulative_tests=2709,numbered_scientific_files=935,
    unique_protected_evidence_files=1201,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=542,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
