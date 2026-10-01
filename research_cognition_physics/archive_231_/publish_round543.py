"""Publish round 543 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_singlet_origin_and_vacuum as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round543_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_543_checks.json')
assert checked['round']==543 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round543_20260930'
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
summary='**第543轮完成：** [中性标量来源与共同质量障碍]({p}research_note_543.md)核现有内涨落不能生成Majorana标量；显式添加σ后，负Higgs代数修正与同一裸质量条件联立，非零top严格双凝聚分支存在尺度比上界，旧平衡点没有重σ。零top平坦族及独立双质量是明确例外。8组与独立终审通过，最新543／2724，941份编号科学文件，1215份保护证据；[核验]({p}research_round_543_checks.json)、[条件更新]({p}unified_physics_condition_ledger_543.md)。运行分裂、完整物理拟合与统一仍开放。'
next_steps='\n\n### 第543轮后：共同高质量边界是否能运行出重轻分离\n\n543对接实中性标量修订：现有代数普通内涨落不生成σ；作为显式新增候选，同一保护对给λ_eff=3q²/T，但共同二次系数下，q>0的严格双凝聚只能q>S，尺度比受限。旧q=S平衡点σ为零曲率，不能直接当重阈值；q=0平坦族及独立质量是明确例外。全部限制只在同尺度、平直常场a4裸势，未排除量子运行或不同反项。\n\n下一正式编号544核同一高能质量边界在含σ、Higgs与Majorana的共同运行下能否形成重轻分离。先核两标量质量β、Yukawa与场归一化、阈值和轻子数保护，不直接照搬自由指定低能双质量的文献拟合。有限代数／允许涨落的扩展并行保留为结构来源条件；若新增反项或输入必须列账，不能称旧条件已自动产生。\n\n阶段目标仍是共同模型的逻辑相融。冻结历史、条件总账与编号报告持续保存；几何、维数、Lorentz、实际记录和资源仍开放，不更改应用目标或设置定时任务。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第543轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 189.' not in body
        body+='\n\n## 189. 标量修订的来源与共同尺度条件\n\n'+block+next_steps
    else:
        assert '\n## 94.' not in body
        body+='\n\n## 94. 结构修订不能省略真空与尺度\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为542／2716','最新科学轮次与检查数为543／2724')
        body=body.replace('完成231—542轮。','完成231—543轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第543轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|543|[中性标量来源与共同质量障碍](research_note_543.md)|[代码](joint_singlet_origin_and_vacuum.py)、'
            '[结果](joint_singlet_origin_and_vacuum_results.json)、[核验](research_round_543_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round543.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=543,cumulative_tests=2724,numbered_scientific_files=941,
    unique_protected_evidence_files=1215,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=544,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
