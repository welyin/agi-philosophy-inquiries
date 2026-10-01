"""Publish round 558; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_gauge_link_scale_matching as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round558_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_558_checks.json')
assert checked['round']==558 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round558_20260930'
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

summary='**第558轮完成：** [共同有效作用与测量尺度]({p}research_note_558.md)证明同一快链路谱面对原κ、b的导数共同给边势与电能均值；瞬时仪器的实际波动和注能却不能仅由有效势替代。构造正规Gauss源族，完整初能量一致有界而固定分辨瞬时电注能随b增长。6组、13式及独立终审通过，最新558／2823，986份编号科学文件，1350份保护证据；[核验]({p}research_round_558_checks.json)、[条件更新]({p}unified_physics_condition_ledger_558.md)。固定背景谱与初态族不等于全物质绝热动力学；慢仪器、时空与统一仍开放。'
next_steps='\n\n### 第558轮后：同一有限时间探针的尺度匹配\n\n同一有效谱面及原耦合导数已经共同固定边势、电能平均，避免另外任意配匹系数。实际瞬时读数仍含高带方差，完整物理源态族显示有限初能量不能免除随能隙增长的测量注能；这只限制给定瞬时仪器，未排除低能有效测量。\n\n后继559检验保留源演化、持续有限时间的同一F探针，核实际读数和源注能是否共同受控。先复用363、421—422、434—435及成熟有限时间工具；给定τ和开关仍须列账，不称自主实现。\n\n给定图、规范群、源准备、全物质低能动力学以及时空／引力来源仍开放。统一目标不变，历史冻结证据保留，阶段未完成。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第558轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 204.' not in body
        body+='\n\n## 204. 同一有效作用、实际记录与仪器尺度\n\n'+block+next_steps
    else:
        assert '\n## 109.' not in body
        body+='\n\n## 109. 共同尺度匹配尚未生成宏观维数\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[559有限时间探针的共同匹配候选](archive_231_/round559_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为557／2817','最新科学轮次与检查数为558／2823')
        body=body.replace('完成231—557轮。','完成231—558轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第558轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|558|[共同有效作用与测量尺度](research_note_558.md)|[代码](joint_gauge_link_scale_matching.py)、'
            '[结果](joint_gauge_link_scale_matching_results.json)、[核验](research_round_558_checks.json)|\n')
        body+='\n**559候选：** [有限时间探针与共同读数](round559_drafts/STATUS.md)，尚未完成；先回查363、421—422与434—435，不重做通用测量结论。\n'
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
    temp=path.with_name(path.name+'.round558.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=558,cumulative_tests=2823,numbered_scientific_files=986,
    unique_protected_evidence_files=1350,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=559,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))

