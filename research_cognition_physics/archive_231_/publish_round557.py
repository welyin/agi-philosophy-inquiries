"""Publish round 557; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_gauge_link_reference as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round557_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_557_checks.json')
assert checked['round']==557 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round557_20260930'
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

summary='**第557轮完成：** [共同规范链路与测量注能]({p}research_note_557.md)在同势双胞元SU(2)模型中，把相邻作用、两端Gauss、配置读口和源反作用接到同一边势；一个E₁预算控制实际二阶读数及测后平均能量。细群仪器有具体部门泄漏，直接不变指针逐分支保约束。6组、14式及独立终审通过，最新557／2817，983份编号科学文件，1340份保护证据；[核验]({p}research_round_557_checks.json)、[条件更新]({p}unified_physics_condition_ledger_557.md)。图、群、指针与理想脉冲仍输入；未生成时空或完成统一。'
next_steps='\n\n### 第557轮后：共同作用和记录的尺度接口\n\n同一动态链路与正四次物质在给定图上已接通局域约束、边势读取和非选择源注能；这个配置读口只需E₁预算，但556完整动量菜单仍保留E₂条件。指针和理想脉冲的来源尚未消除，不把源注能改称全装置耗散。\n\n后继558先核同一链路消去后的作用和记录能否共同闭合。复用363压缩测量矩及434—435有效作用，不重报一般投影恒等式；区分链路常模、真实低能部门及固定背景近似。优先减少共享输入和发现跨部门限制，不细扫单边精度。\n\n保留历史冻结材料及原合并顺序；最新条件状态见557总账，旧压缩表中下一轮548属于历史。统一目标不变，阶段未完成。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第557轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 203.' not in body
        body+='\n\n## 203. 同一规范边势的约束、读口与源反作用\n\n'+block+next_steps
    else:
        assert '\n## 108.' not in body
        body+='\n\n## 108. 相邻规范读数与维数来源仍须区分\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[558同一链路的作用与记录尺度候选](archive_231_/round558_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为556／2811','最新科学轮次与检查数为557／2817')
        body=body.replace('完成231—556轮。','完成231—557轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第557轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|557|[共同规范链路与测量注能](research_note_557.md)|[代码](joint_gauge_link_reference.py)、'
            '[结果](joint_gauge_link_reference_results.json)、[核验](research_round_557_checks.json)|\n')
        body+='\n**558候选：** [共同作用与记录的尺度闭合](round558_drafts/STATUS.md)，尚未完成；先回查363与434—435，不重做通用压缩结论。\n'
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
    temp=path.with_name(path.name+'.round557.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=557,cumulative_tests=2817,numbered_scientific_files=983,
    unique_protected_evidence_files=1340,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=558,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))

