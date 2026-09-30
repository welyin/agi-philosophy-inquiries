"""Publish round 534 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import spectral_background_feedback as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round534_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_534_checks.json')
assert checked['round']==534 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round534_20260930'
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
summary=('**第534轮完成：** [共同背景与完整谱尺度检验]({p}research_note_534.md)'
    '分类同一a₄截断的常场几何解：M=0无非零径向解，Majorana分支可有。'
    '进一步以完整Dirac谱证明：严格递减截止、固定Λ／常Φ且自由半径的圆S⁴裸作用没有有限半径驻点，'
    '截断球面解不能升级为完整真空。8组与独立终审通过，最新534／2655，914份编号科学文件，'
    '1151份保护证据；[核验]({p}research_round_534_checks.json)、[条件更新]({p}unified_physics_condition_ledger_534.md)。'
    '只排除所列Euclidean常场分支；Lorentz、动态场、体积约束及整体统一仍开放。')
next_steps='\n\n### 第534轮后：明确允许的尺度变化及有效作用\n\n534已用完整谱排除固定截止、常物质场、无其它源且允许自由半径的圆S⁴驻点，不能再在这个分支扫描参数。M=0与Majorana截断解的区别已分类；完整作用的严格单调来自每一正谱模式。非严格紧支撑截止、固定体积或额外量子有效作用属于不同条件，尚未排除。\n\n下一正式编号535优先对接Chamseddine–Connes 0812.0165引言式(9)的体积约束，以及既有量子有效作用。先区分三个已出现但不等同的条件：532有限代数的unimodularity；谱几何的四维体积约束；306的Jacobson局部空间球固定体积。不得凭同名或都限制体积就直接相认。检验它们是否允许同一物质背景及规范归一化、是否只把宇宙常数移为未定乘子，以及物理体积／谱计数与内部资源之间还欠哪些明确映射。\n\n这一步允许成熟物理反向约束认知定义，但输入与推导分列。533裸系数界不自动限制加了乘子的有效宇宙项；同样，取消径向变分也不自动构造实际测量、选出3+1或解决观测Λ。新增物理条件和竞争分支回填总账，保留历史和冻结证据，不改应用目标或定时任务。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第534轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 180.' not in body
        body+='\n\n## 180. 共同常场背景与完整作用的尺度障碍\n\n'+block+next_steps
    else:
        assert '\n## 85.' not in body
        body+='\n\n## 85. 完整谱排除一个常场几何分支\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为533／2647','最新科学轮次与检查数为534／2655')
        body=body.replace('完成231—533轮。','完成231—534轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第534轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|534|[完整谱与共同几何背景](research_note_534.md)|[代码](spectral_background_feedback.py)、'
            '[结果](spectral_background_feedback_results.json)、[核验](research_round_534_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round534.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=534,cumulative_tests=2655,numbered_scientific_files=914,
    unique_protected_evidence_files=1151,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=535,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
