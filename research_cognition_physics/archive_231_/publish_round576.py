"""Publish round 576 and the consolidated work order without changing the goal."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round576_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_576_checks.json')
assert checked['round']==576 and checked['all_reported_checks_passed']
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round576_20261001'
resuming=folder.exists()
if resuming:
    manifest=core.read(folder/'manifest.json')
    raws={p:(folder/manifest[str(p.relative_to(ROOT))]['snapshot']).read_bytes() for p in paths}
    for p,raw in raws.items():
        assert hashlib.sha256(raw).hexdigest()==manifest[str(p.relative_to(ROOT))]['sha256']
else:
    raws={p:p.read_bytes() for p in paths};folder.mkdir(exist_ok=False);manifest={}
    for p,raw in raws.items():
        key=('root_' if p.parent==ROOT else 'research_' if p.parent==RESEARCH else 'archive_')+p.name
        with (folder/key).open('xb') as stream:stream.write(raw)
        manifest[str(p.relative_to(ROOT))]=dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    with (folder/'manifest.json').open('x',encoding='utf8') as stream:
        json.dump(manifest,stream,ensure_ascii=False,indent=2)

summary='**第576轮完成：** [自主指针的物质谱与新增通道]({p}research_note_576.md)核575同一作用：原两径向质量保留，新增中性无质量模，并由连接曲率固定非零重→轻＋指针树级通道。原非均匀源的逐点零指针动量不被保持；裸加质量与相关势补全具有不同低能结果。六组、十四式及独立终审通过，最新576／2940，1040份编号科学文件，1523份保护证据；[核验]({p}research_round_576_checks.json)、[条件账]({p}unified_physics_condition_ledger_576.md)。限常真空有效作用与树级谱，未作实验拟合或量子连续证明，统一仍开放。'
order='**执行顺序：** 已完成新增探针的兼容筛查；[576条件账]({p}unified_physics_condition_ledger_576.md)明确额外物理输入与可选补全。接续原物质区域记录入口，优先检验能否减少新增字段及相互作用；不继续精扫χ质量和寿命。'
next_steps='\n\n### 第576轮后：由装置的物质限制转向原区域记录\n\n同一指针连接在记录机制之外，还固定额外物质模及三点作用；二阶旧质量不变不能解释成原物理完全不变。裸pin改变低能动能；相关pin保谱但增加新势，尚未复验非均匀几何或记录。\n\n按先合并兼容条件、再抓联合卡点的顺序，577先核原574完整有限图Hamiltonian能否在不同相邻区域之间传递相位信息到有界读口，不另加χ或新的测量耦合。准备、终读、图和量子化仍计输入；Gauss物理空间不能默认张量分解。\n\n只有完整算符、正规态、实际概率和资源的共同命题通过才计新轮。该新准备的应力不能直接继承572初值；共同几何和有限ℏ自洽仍须接通。保持维数、引力来源、手征物质和共同尺度的全局目标，不把仪器参数微调作为自动接续。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第576轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+order.format(p=prefix)+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 222.' not in body
        body+='\n\n## 222. 探针物质限制与同一候选的相容范围\n\n'+block+next_steps
    else:
        assert '\n## 127.' not in body
        body+='\n\n## 127. 物质谱与装置身份不等于空间维数选择\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[577原物质区域记录候选](archive_231_/round577_drafts/STATUS.md)已登记，尚未计完成轮次；[最新条件账](archive_231_/unified_physics_condition_ledger_576.md)。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为574／2927','最新科学轮次与检查数为576／2940')
        body=body.replace('完成231—574轮。','完成231—576轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第576轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|576|[自主指针的物质谱与新增通道](research_note_576.md)|[代码](joint_pointer_matter_compatibility.py)、'
            '[结果](joint_pointer_matter_compatibility_results.json)、[核验](research_round_576_checks.json)|\n')
        body+='\n[条件账](unified_physics_condition_ledger_576.md)已更新。**577候选：** [原物质区域的自主记录](round577_drafts/STATUS.md)，尚未计完成。\n'
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
support_links=0
support_files=('unified_physics_condition_ledger_576.md','round577_drafts/STATUS.md')
for name in support_files:
    path=HERE/name
    for link in core.link_parser()(path.read_text('utf8')):
        assert (path.parent/link).resolve().exists(),(path,link)
        support_links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body:continue
    temp=path.with_name(path.name+'.round576.tmp')
    if temp.exists():assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=576,cumulative_tests=2940,numbered_scientific_files=1040,
    unique_protected_evidence_files=1523,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=577,
    supporting_document_links=support_links,
    supporting_document_hashes={name:core.digest(HERE/name) for name in support_files},
    consolidation_not_counted_as_scientific_round=True,active_goal_unchanged=True,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
