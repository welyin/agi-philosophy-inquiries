"""Publish round 554 with existing navigation/history and concurrent-write checks."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_matter_reference_quantum_readout as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round554_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_554_checks.json')
assert checked['round']==554 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round554_20260930'
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

summary='**第554轮完成：** [同物质参考的量子联合读取]({p}research_note_554.md)在明确有限截断中联立h、s导数菜单与Gaussian共同记录；保留源偏差，给同平均全H能量和同协方差却有无界平方读口风险的反例。6组、12式及独立终审通过，最新554／2799，974份编号科学文件，1310份保护证据；[核验]({p}research_round_554_checks.json)、[条件更新]({p}unified_physics_condition_ledger_554.md)。仪器与来源仍输入，标定不是连续时空定位；完整物质／引力与统一仍开放。'
next_steps='\n\n### 第554轮后：共同读口与来源高阶矩\n\n同一约化有限h、s菜单已有精确对易限制和明确共同POVM，实际平方读取需要高阶来源矩；平均能量和二阶协方差不足。源偏差、仪器噪声、有限宽度、径向截断和真实坐标映射分别保留。前置真空一圈审查属于未编号接口，不以自由V₀匹配签收新理论。\n\n后继555先核同一相互作用H下的来源矩控制：有限〈H²〉或有来源的热态是否足够，能否在非Gaussian演化后仍控制平方读口风险。不得把非对易算符H≥T直接平方，不只继续优化Gaussian精度。\n\n全局统一目标未改，历史冻结材料保留，阶段未完成。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第554轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 200.' not in body
        body+='\n\n## 200. 同物质参考的共同读取与来源条件\n\n'+block+next_steps
    else:
        assert '\n## 105.' not in body
        body+='\n\n## 105. 有限量子参考标定不等于时空坐标生成\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[555来源高阶矩候选](archive_231_/round555_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为553／2793','最新科学轮次与检查数为554／2799')
        body=body.replace('完成231—553轮。','完成231—554轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第554轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|554|[同物质参考的量子联合读取](research_note_554.md)|[代码](joint_matter_reference_quantum_readout.py)、'
            '[结果](joint_matter_reference_quantum_readout_results.json)、[核验](research_round_554_checks.json)|\n')
        body+='\n**555候选：** [同一H下的来源矩控制](round555_drafts/STATUS.md)，尚未完成；来源存在与实际准备分开，不以一般矩公式凑轮次。\n'
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
    temp=path.with_name(path.name+'.round554.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=554,cumulative_tests=2799,numbered_scientific_files=974,
    unique_protected_evidence_files=1310,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=555,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
