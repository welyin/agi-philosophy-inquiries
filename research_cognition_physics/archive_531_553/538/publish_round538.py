"""Publish round 538 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import protected_pair_rg as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round538_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_538_checks.json')
assert checked['round']==538 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round538_20260930'
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
summary='**第538轮完成：** [保护复相关分支的共同尺度运行]({p}research_note_538.md)把537存活边界接入规范、top、ν及Higgs一圈运行和同时树级阈值。轻质量相消保留，但Hermitian回路贡献不消失；错误省略ν使低尺度输出偏移。8组与独立终审通过，最新538／2687，926份编号科学文件，1179份保护证据；[核验]({p}research_round_538_checks.json)、[条件更新]({p}unified_physics_condition_ledger_538.md)。有限一圈匹配、非零质量及现实拟合未完成，整体目标开放。'
next_steps='\n\n### 第538轮后：非零轻质量与共同边界的接续\n\n538确认Y=(c,ic,0)/√2的受保护配对可以保持C₅=0，但仍影响top与Higgs运行；重场须按同一运行质量阈值匹配，不能据零轻质量删除它们的回路。四个阈值和五个共同边界样本仅作一圈维四＋树级匹配诊断，不构成全参数域排除或现实拟合。\n\n下一正式编号539先以独立成熟来源核清非零Weinberg项的四次系数，再检验小轻子数破缺能否同时实现非零质量和同一规范／Higgs边界。来源因子核对不单独虚增轮次；需要给出真实新增联合条件、构造或反例。保留有限一圈阈值、完整味谱、实际测量与共同几何等开放条件。不另调规范权，不添加无来源物质，不改应用目标、不设定时任务。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第538轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 184.' not in body
        body+='\n\n## 184. 保护复结构与实际运行的联立\n\n'+block+next_steps
    else:
        assert '\n## 89.' not in body
        body+='\n\n## 89. 尺度接续与几何生成的边界\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为537／2679','最新科学轮次与检查数为538／2687')
        body=body.replace('完成231—537轮。','完成231—538轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第538轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|538|[保护复相关分支的共同运行](research_note_538.md)|[代码](protected_pair_rg.py)、'
            '[结果](protected_pair_rg_results.json)、[核验](research_round_538_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round538.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=538,cumulative_tests=2687,numbered_scientific_files=926,
    unique_protected_evidence_files=1179,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=539,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
