"""Publish round 575 and the consolidated work order without changing the goal."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round575_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_575_checks.json')
assert checked['round']==575 and checked['all_reported_checks_passed']
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round575_20261001'
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

summary='**第575轮完成：** [内部物质时钟、自主指针与共同几何]({p}research_note_575.md)明示增加圆值指针及局域连接，将自主响应、完整应力交换、共同经典约束和保留的局部参考接通；固定图正规Gauss态有非零短时有界读口响应。七组、十七式及独立终审通过，最新575／2934，1037份编号科学文件，1515份保护证据；[核验]({p}research_round_575_checks.json)、[条件账]({p}unified_physics_condition_ledger_575.md)。探针、耦合、准备及末读仍输入，有限ℏ应力与几何尚未自洽，统一仍开放。'
order='**执行顺序：** [先收束可兼容条件，再处理关键接口]({p}joint_condition_compression_update_575.md)。共同实现与独立假设减少分别记账；已有证明直接复用，先核新增探针与原物质谱／有效尺度的兼容，再选择牵动多个部门的缺口。整理不增轮次，不继续单一仪器参数优化。'
next_steps='\n\n### 第575轮后：收束兼容结果，再筛查跨部门卡点\n\n575将原物质时钟、明确新增的自主探针、共同经典应力和局部参考相连，并给固定图量子有界效果的短时响应。必须分别保留连续经典与固定图量子范围；这尚非有限ℏ的量子应力—几何自洽。\n\n执行顺序按最新合并清单：可直接继承的Gauss、正性、约束和参考结果不重做；共同实现增加不等于独立输入减少。576先对同一作用作新增探针的谱／低能兼容筛查，区分额外基本场与已有物质组成的有效装置。如果只是成熟工具应用，就记入口审计，不新增轮次；发现真正冲突及时检验替代候选。\n\n随后按共同模型的限制力选择有限ℏ几何、共同尺度与区域拼接等缺口。3+1、Lorentz、Einstein作用来源与手征物质仍在联合主账中。停止指针强度、波包、网格和记录寿命精扫，完整目标不改，未设置定时任务。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第575轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+order.format(p=prefix)+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 221.' not in body
        body+='\n\n## 221. 自主记录与共同经典几何的明示扩展\n\n'+block+next_steps
    else:
        assert '\n## 126.' not in body
        body+='\n\n## 126. 自主目标连接不选择时空维数\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[576新增探针与共同物质谱候选](archive_231_/round576_drafts/STATUS.md)已登记，尚未计完成轮次；[最新条件合并清单](archive_231_/joint_condition_compression_update_575.md)。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为574／2927','最新科学轮次与检查数为575／2934')
        body=body.replace('完成231—574轮。','完成231—575轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第575轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|575|[内部物质时钟、自主指针与共同几何](research_note_575.md)|[代码](joint_autonomous_pointer_geometry.py)、'
            '[结果](joint_autonomous_pointer_geometry_results.json)、[核验](research_round_575_checks.json)|\n')
        body+='\n[条件合并清单](joint_condition_compression_update_575.md)为执行整理，不增轮次。**576候选：** [新增探针与原谱／低能兼容](round576_drafts/STATUS.md)，尚未计完成。\n'
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
support_links=0
support_files=('joint_condition_compression_update_575.md','round576_drafts/STATUS.md')
for name in support_files:
    path=HERE/name
    for link in core.link_parser()(path.read_text('utf8')):
        assert (path.parent/link).resolve().exists(),(path,link)
        support_links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body:continue
    temp=path.with_name(path.name+'.round575.tmp')
    if temp.exists():assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=575,cumulative_tests=2934,numbered_scientific_files=1037,
    unique_protected_evidence_files=1515,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=576,
    supporting_document_links=support_links,
    supporting_document_hashes={name:core.digest(HERE/name) for name in support_files},
    consolidation_not_counted_as_scientific_round=True,active_goal_unchanged=True,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
