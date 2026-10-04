"""Publish round 565; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_full_matter_large_b_limit as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round565_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_565_checks.json')
assert checked['round']==565 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round565_20261001'
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

summary='**第565轮完成：** [全物质快链路极限与相互作用损失]({p}research_note_565.md)证明固定a、k、W而b增大时，完整无界物质与同一有限时径向仪器对固定P₀输入有强极限；预算可一致，但极限未测源跨端分解。五组、十四式及独立终审通过，最新565／2863，1007份编号科学文件，1418份保护证据；[核验]({p}research_round_565_checks.json)、[条件账]({p}unified_physics_condition_ledger_565.md)。限制固定参数冻结分支，RP仪器仍可联端，未证明空间、GR或统一完成。'
next_steps='\n\n### 第565轮后：继承可相容的合并，先核共同场论接口\n\n同一规范物质、来源、实际径向记录、预算和受控极限已得到联合证书；固定参数极限会删除跨端物质作用，这一限制不靠继续加强冻结解决。优先保留564完整联合态，并与成熟格点规范—标量理论核对区域、作用、尺度及测量的共同映射。\n\n按用户确认的顺序，已有同模型结论直接继承；先处理映射明确、可减少独立输入的条件，再集中到连续尺度、同物质装置与几何反馈等关键缺口。严格不相容的拼接及时排除，单探针精度和常数优化后置。成熟定理复述与表格整理不计新轮次。\n\n同步耦合／质量匹配保留为替代，须另核完整物质稳定和动力学。给定图、群、维数、势常数及新作用明确列账；统一目标保持，历史证据不覆盖。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第565轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 211.' not in body
        body+='\n\n## 211. 受控约化与源相互作用\n\n'+block+next_steps
    else:
        assert '\n## 116.' not in body
        body+='\n\n## 116. 有效来源与原有传播能否共同保留\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[566完整结构与共同场论接口候选](archive_231_/round566_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为564／2858','最新科学轮次与检查数为565／2863')
        body=body.replace('完成231—564轮。','完成231—565轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第565轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|565|[全物质强极限与源作用损失](research_note_565.md)|[代码](joint_full_matter_large_b_limit.py)、'
            '[结果](joint_full_matter_large_b_limit_results.json)、[核验](research_round_565_checks.json)|\n')
        body+='\n**566候选：** [保留完整结构的共同场论接口](round566_drafts/STATUS.md)，尚未完成；优先对接成熟理论并核共同尺度，保留同步耦合匹配为替代，不重复固定参数冻结或单探针常数优化。\n'
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
    temp=path.with_name(path.name+'.round565.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=565,cumulative_tests=2863,numbered_scientific_files=1007,
    unique_protected_evidence_files=1418,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=566,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))

