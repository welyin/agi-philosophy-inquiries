"""Publish round 566; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_shared_vertex_gluing as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round566_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_566_checks.json')
assert checked['round']==566 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round566_20261001'
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

summary='**第566轮完成：** [共享物质节点、联合记录与交叉反作用]({p}research_note_566.md)精确约化同一有限图全H，证明共享节点角能须保留原动力学固定的交叉项；两源同逐边态与同平均全H，共同Gaussian读数和瞬时源注能仍不同。六组、十四式及独立终审通过，最新566／2869，1010份编号科学文件，1426份保护证据；[核验]({p}research_round_566_checks.json)、[条件账]({p}unified_physics_condition_ledger_566.md)。逐来源存在有限时区分；不是局部层析失效，图、群、装置及连续几何仍有输入，统一未完成。'
next_steps='\n\n### 第566轮后：准确拼接后，核共同物质传播与尺度\n\n共享节点、完整规范源、联合读口和反作用已在同一有限图中连接；交叉项由原a、k固定，不另设任意边—边耦合。完整联合态、指针协方差与实际装置边界继续保留。\n\n后继567先核548—549连续h＋s参考与557—566图源的导数部门是否匹配。原s现场动能和门户势能否承担共同传播主部，须从当前同一势、同一真空的图线性化检验；339已有一般多物种判据只继承，不重报。需要补s边势或改变尺度时明确列为新增条件；量子诱导与引力混合另核。\n\n再审计闭环磁作用及共同连续尺度。先合并对象明确的条件，再攻关键缺口；不把成熟公式复述、单探针常数优化计成新轮次。图、群、维数和完整几何反馈仍开放，统一目标与冻结历史保持。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第566轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 212.' not in body
        body+='\n\n## 212. 共享节点的精确动力学与联合操作\n\n'+block+next_steps
    else:
        assert '\n## 117.' not in body
        body+='\n\n## 117. 给定图上的共同载体与联合资料\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[567共同物质传播与尺度候选](archive_231_/round567_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为565／2863','最新科学轮次与检查数为566／2869')
        body=body.replace('完成231—565轮。','完成231—566轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第566轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|566|[共享节点与联合记录代价](research_note_566.md)|[代码](joint_shared_vertex_gluing.py)、'
            '[结果](joint_shared_vertex_gluing_results.json)、[核验](research_round_566_checks.json)|\n')
        body+='\n**567候选：** [共同物质传播与连续尺度](round567_drafts/STATUS.md)，尚未完成；已有同势线性化入口复算，先核h＋s导数部门映射，再审共同尺度与闭环作用；不把补入的梯度当已生成。\n'
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
    temp=path.with_name(path.name+'.round566.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=566,cumulative_tests=2869,numbered_scientific_files=1010,
    unique_protected_evidence_files=1426,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=567,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))

