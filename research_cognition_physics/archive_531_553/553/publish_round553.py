"""Publish round 553 with existing navigation/history and concurrent-write checks."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_curved_matter_renormalization as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round553_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_553_checks.json')
assert checked['round']==553 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round553_20260930'
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

summary='**第553轮完成：** [曲率运行与全尺度谱关系]({p}research_note_553.md)将同一物质的一圈曲率反项接入共同模型；证明Weyl／弱耦合谱关系的偏离斜率为(12N²+299)/(120·16π²)>0，代数目和正谱权不能消去。共同共形面本阶保持，单尺度谱边界仍可用。6组、11式及独立终审通过，最新553／2793，971份编号科学文件，1300份保护证据；[核验]({p}research_round_553_checks.json)、[条件更新]({p}unified_physics_condition_ledger_553.md)。仅限物质一圈MS、经典度规；量子真空与统一仍开放。'
next_steps='\n\n### 第553轮后：可保持条件与单尺度匹配\n\n共同共形曲率面在声明的一圈物质系统中保持；未修订的Weyl／弱耦合谱关系不能跨尺度保持，代数目和谱权不能修复。保留单尺度谱边界，再按同一物质规则运行高阶曲率参数。真空树势的辅助尺度变化由有效势显式对数抵消，不能解释为物理真空的尺度依赖。\n\n后继554先核同一物质的有限量子真空匹配、正F与曲率运行能否共同实现，优先减少独立输入；若只是自由V₀调节或非退化驻点微扰延续，记成熟接口，不据此凑轮次。随后回到全局条件账，避免单一部门无限细化。\n\n统一目标未改；所有冻结证据保留，阶段尚未完成。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第553轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 199.' not in body
        body+='\n\n## 199. 共同曲率运行与谱关系的尺度边界\n\n'+block+next_steps
    else:
        assert '\n## 104.' not in body
        body+='\n\n## 104. 共形保持与尺度反项不等于维数生成\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[554共同量子匹配候选](archive_231_/round554_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为552／2787','最新科学轮次与检查数为553／2793')
        body=body.replace('完成231—552轮。','完成231—553轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第553轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|553|[曲率运行与全尺度谱关系](research_note_553.md)|[代码](joint_curved_matter_renormalization.py)、'
            '[结果](joint_curved_matter_renormalization_results.json)、[核验](research_round_553_checks.json)|\n')
        body+='\n**554候选：** [共同量子匹配](round554_drafts/STATUS.md)，尚未完成；成熟接口或自由真空调节不另计轮次。\n'
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
    temp=path.with_name(path.name+'.round553.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=553,cumulative_tests=2793,numbered_scientific_files=971,
    unique_protected_evidence_files=1300,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=554,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
