"""Publish round 555 with existing navigation/history and concurrent-write checks."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_matter_energy_moment_control as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round555_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_555_checks.json')
assert checked['round']==555 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round555_20260930'
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

summary='**第555轮完成：** [共同能量矩与非Gaussian参考记录]({p}research_note_555.md)证明同一正四次势有限模型中，一个二阶能量矩预算控制四类实际记录，且在未测静态演化下保持；独立正测量密度积分核非Gaussian源，并给无约束势反例。6组、12式及独立终审通过，最新555／2805，977份编号科学文件，1319份保护证据；[核验]({p}research_round_555_checks.json)、[条件更新]({p}unified_physics_condition_ledger_555.md)。来源预算与仪器仍输入，有限风险不等于稳定坐标；统一仍开放。'
next_steps='\n\n### 第555轮后：共同来源矩充分性与规范物质接口\n\n同一有限相互作用H的二阶矩预算，已接到任意非Gaussian源的四个实际参考读数；静态未测演化保留预算。来源准备、标定稳定和仪器物质身份仍未生成，也未取得连续或规模一致极限。此接口的充分性已闭合，不以继续压小常数代替共同模型推进。\n\n后继556回查完整Higgs规范不变观测、径向量子化与实际仪器的映射，复用360—365区域约束及421守恒量审计。不得将全实轴h的正则量子化直接当作严格正半轴径向子部门，不重复一般Gauss或WAY障碍。\n\n全局统一目标未改，历史冻结材料保留，阶段未完成。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第555轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 201.' not in body
        body+='\n\n## 201. 共同能量矩与实际参考来源\n\n'+block+next_steps
    else:
        assert '\n## 106.' not in body
        body+='\n\n## 106. 源风险守恒不等于坐标标定稳定\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[556规范物质与仪器候选](archive_231_/round556_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为554／2799','最新科学轮次与检查数为555／2805')
        body=body.replace('完成231—554轮。','完成231—555轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第555轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|555|[共同能量矩与非Gaussian参考记录](research_note_555.md)|[代码](joint_matter_energy_moment_control.py)、'
            '[结果](joint_matter_energy_moment_control_results.json)、[核验](research_round_555_checks.json)|\n')
        body+='\n**556候选：** [完整规范物质与实际仪器](round556_drafts/STATUS.md)，尚未完成；径向量子化与规范读口须映射，不重复一般约束障碍。\n'
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
    temp=path.with_name(path.name+'.round555.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=555,cumulative_tests=2805,numbered_scientific_files=977,
    unique_protected_evidence_files=1319,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=556,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
