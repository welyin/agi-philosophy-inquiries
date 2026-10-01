"""Publish round 556 with existing navigation/history and concurrent-write checks."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_higgs_invariant_reference_readout as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round556_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_556_checks.json')
assert checked['round']==556 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round556_20260930'
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

summary='**第556轮完成：** [Higgs不变读取与后态约束]({p}research_note_556.md)在单胞元全局SU(2)部门中给同势径向域反例及合法不变POVM，共用完整能量矩预算；具体兼容仪器却有11/12后态离开singlet，抽象平方根分箱仪器可保部门。6组、12式及独立终审通过，最新556／2811，980份编号科学文件，1330份保护证据；[核验]({p}research_round_556_checks.json)、[条件更新]({p}unified_physics_condition_ledger_556.md)。未实现局域Gauss、物质探针或时空坐标；统一仍开放。'
next_steps='\n\n### 第556轮后：从单胞元不变读取转向共同规范链路\n\n完整物质的表示、势、能量矩与不变概率读取已在单胞元模型中接通。旧平直径向读口不能直接继承；原Cartesian测量后制备又不保singlet，故后态和实现条件须保留。有限分箱平方根仪器只证明抽象相容，不签收物质来源与测后能源。\n\n后继557优先检验相邻物质、动态SU(2)链路、协变梯度与实际联合读数。复用360—365已有区域约束；不能把未知链路免费当经典平行运输，或把尖锐群元读取的后态代价略去。只有新的同模型连接才增加轮次。\n\n统一目标未改，历史冻结材料保留，阶段未完成。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第556轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 202.' not in body
        body+='\n\n## 202. 完整不变物质、概率与状态更新\n\n'+block+next_steps
    else:
        assert '\n## 107.' not in body
        body+='\n\n## 107. 内部规范不变读数不等于时空坐标\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[557规范链路与相邻参考候选](archive_231_/round557_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为555／2805','最新科学轮次与检查数为556／2811')
        body=body.replace('完成231—555轮。','完成231—556轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第556轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|556|[Higgs不变读取与后态约束](research_note_556.md)|[代码](joint_higgs_invariant_reference_readout.py)、'
            '[结果](joint_higgs_invariant_reference_readout_results.json)、[核验](research_round_556_checks.json)|\n')
        body+='\n**557候选：** [规范链路与相邻参考](round557_drafts/STATUS.md)，尚未完成；保留动态链路与仪器代价，复用既有Gauss拼接。\n'
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
    temp=path.with_name(path.name+'.round556.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=556,cumulative_tests=2811,numbered_scientific_files=980,
    unique_protected_evidence_files=1330,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=557,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
