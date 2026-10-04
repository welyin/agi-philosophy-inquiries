"""Publish round 550 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_dynamic_curvature_scale as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round550_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_550_checks.json')
assert checked['round']==550 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round550_20260930'
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
summary='**第550轮完成：** [共同裸模型的动态曲率界]({p}research_note_550.md)把完整非最小场方程与同一谱矩、规范归一联立，证明固定非零弱耦合下逐点R/F≥6g_w²；减弱参考梯度或逼近相对F边界不能消除此界。8组及独立终审通过，最新550／2772，962份编号科学文件，1273份保护证据；[核验]({p}research_round_550_checks.json)、[条件更新]({p}unified_physics_condition_ledger_550.md)。非单调截止与独立真空匹配给明确例外；只限制所列裸分支，不是完整谱理论无解或统一完成。'
next_steps='\n\n### 第550轮后：联合曲率界与最少匹配输入\n\n549的局部参考与几何共同实现保留；550证明同一正单调截止裸分支具有动态曲率／系数下界，不能仅靠降低参考振幅进入任意弱曲率极限。完整谱近似仍未受控；R/F不是已经匹配的实测Newton尺度。\n\n下一正式编号551先检验仅增加一个明示真空匹配常数是否已足以同时释放物质／引力尺度。回查543—544物质层级与533规范归一，只核新增引力接口；任何新增独立质量、Einstein或非单调截止输入分支都要单列，不能只在一个部门消掉问题后称统一完成。\n\n研究目标、27项条件入口和所有冻结历史保持；继续按共同模型合并条件与处理跨部门卡点。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第550轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 196.' not in body
        body+='\n\n## 196. 同一裸模型的动态曲率尺度限制\n\n'+block+next_steps
    else:
        assert '\n## 101.' not in body
        body+='\n\n## 101. 动态曲率界不等于维数或引力生成\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[551物质与引力尺度草稿](archive_231_/round551_drafts/STATUS.md)已开启，尚未计入完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为549／2764','最新科学轮次与检查数为550／2772')
        body=body.replace('完成231—549轮。','完成231—550轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第550轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|550|[共同裸模型的动态曲率界](research_note_550.md)|[代码](joint_dynamic_curvature_scale.py)、'
            '[结果](joint_dynamic_curvature_scale_results.json)、[核验](research_round_550_checks.json)|\n')
        body+='\n**551准备中：** [真空匹配后的物质与引力尺度](round551_drafts/STATUS.md)，尚未完成，不纳入上述累计数。\n'
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
    temp=path.with_name(path.name+'.round550.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=550,cumulative_tests=2772,numbered_scientific_files=962,
    unique_protected_evidence_files=1273,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=551,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
