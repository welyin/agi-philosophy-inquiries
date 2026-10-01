"""Publish round 549 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_reference_gravity_constraints as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round549_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_549_checks.json')
assert checked['round']==549 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round549_20260930'
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
summary='**第549轮完成：** [同一物质参考与动态几何]({p}research_note_549.md)在明确四维二导数分支中，构造同一非最小h＋s及几何的局部解析解：满足引力初始约束，并保留满秩关系参考与完整势常数。8组及独立终审通过，最新549／2764，959份编号科学文件，1264份保护证据；[核验]({p}research_round_549_checks.json)、[条件更新]({p}unified_physics_condition_ledger_549.md)。这是局部共同实现；完整谱近似、状态准备和量子读取仍开放，未生成维数或Einstein作用。'
next_steps='\n\n### 第549轮后：动态参考已相容，审查共同曲率尺度\n\n同一非最小h＋s的局部参考与引力约束已在二导数解析分支共同实现；保留势常数、混合动量与真实场加速度。局部存在不是全局有限资源构造，也不是完整谱作用的受控近似。\n\n下一正式编号550先核同一裸模型的非最小场方程是否强制动态曲率尺度，并检验改变参考振幅、取向或单位能否改变它。继承533静态系数界，不另加未登记真空抵消。若裸分支不受控，限定结论后选择明示的有效匹配分支。\n\n继续按用户要求先合并共同条件，再处理跨部门卡点；既有统一目标、27项审计入口和冻结历史均保留，不设置定时任务。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第549轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 195.' not in body
        body+='\n\n## 195. 物质参考与动态几何的局部约束相容\n\n'+block+next_steps
    else:
        assert '\n## 100.' not in body
        body+='\n\n## 100. 动态参考共同实现仍不选择维数\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[550动态曲率尺度草稿](archive_231_/round550_drafts/STATUS.md)已开启，尚未计入完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为548／2756','最新科学轮次与检查数为549／2764')
        body=body.replace('完成231—548轮。','完成231—549轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第549轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|549|[同一物质参考与动态几何](research_note_549.md)|[代码](joint_reference_gravity_constraints.py)、'
            '[结果](joint_reference_gravity_constraints_results.json)、[核验](research_round_549_checks.json)|\n')
        body+='\n**550准备中：** [同一裸模型的动态曲率尺度](round550_drafts/STATUS.md)，尚未完成，不纳入上述累计数。\n'
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
    temp=path.with_name(path.name+'.round549.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=549,cumulative_tests=2764,numbered_scientific_files=959,
    unique_protected_evidence_files=1264,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=550,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
