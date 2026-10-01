"""Publish round 563; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_direct_radial_record as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round563_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_563_checks.json')
assert checked['round']==563 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round563_20260930'
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

summary='**第563轮完成：** [同一物质源的直接径向记录与共同预算]({p}research_note_563.md)在原全H与560同一Gauss来源中，直接R仪器无须角向对齐即可产生有限时动态记录；弱径向域、实际阈值、源反作用、整体半有界和指针／开关预算共同核清。新RP耦合及补偿末读仍输入，充分成本保守。6组、14式及独立终审通过，最新563／2853，1001份编号科学文件，1398份保护证据；[核验]({p}research_round_563_checks.json)、[条件账]({p}unified_physics_condition_ledger_563.md)。未把单径向记录当四坐标、连续极限或统一完成。'
next_steps='\n\n### 第563轮后：从已合并径向记录转向共同几何菜单\n\n同一物质源、规范约束、直接径向动态记录与完整预算现在可共同实现；直接R仪器分支省去角向冻结准备，但没有消除新测量耦合和自治装置输入。继续按先合并已证条件、再攻关键连接的顺序，单探针常数优化后置。\n\n已回查548的在壳导数菜单及554的四参数联合读口：不重复多加指针的参数辨识。后继优先审计原全H中径向预测是否仍依赖内部取向／链路资料，及该约化与同物质几何菜单的映射。若只是旧工具推论，只更新对象审计；真正联合实现或限定反例才编号。\n\n统一目标不变。C04／C21、连续尺度、维数及完整几何反馈仍开放；冻结历史和既有反例全部保留。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第563轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 209.' not in body
        body+='\n\n## 209. 直接径向记录与共同来源\n\n'+block+next_steps
    else:
        assert '\n## 114.' not in body
        body+='\n\n## 114. 径向可读不等于完整空间参考\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[564共同几何菜单与径向动力学候选](archive_231_/round564_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为562／2847','最新科学轮次与检查数为563／2853')
        body=body.replace('完成231—562轮。','完成231—563轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第563轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|563|[同一源直接径向记录](research_note_563.md)|[代码](joint_direct_radial_record.py)、'
            '[结果](joint_direct_radial_record_results.json)、[核验](research_round_563_checks.json)|\n')
        body+='\n**564候选：** [径向动力学与共同几何菜单](round564_drafts/STATUS.md)，尚未完成；先复用548／554和旧记录闭合审计，核同一全H中的资料依赖，不重复参数满秩或优化单探针常数。\n'
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
    temp=path.with_name(path.name+'.round563.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=563,cumulative_tests=2853,numbered_scientific_files=1001,
    unique_protected_evidence_files=1398,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=564,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))

