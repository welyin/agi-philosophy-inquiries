"""Publish round 545 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_singlet_spectrum_and_weights as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round545_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_545_checks.json')
assert checked['round']==545 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round545_20260930'
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
summary='**第545轮完成：** [共同真空的质量与测量成分]({p}research_note_545.md)联立完整中性谱、径向质量和Higgs／活跃轻子权重；共同裸边界中轻标量越轻，其Higgs成分越少、重对活跃成分越多。给Schur静态响应与实际谱的精确区别，并复用544轨道证明裸界不能跨尺度沿用。7组与独立终审通过，最新545／2739，947份编号科学文件，1230份保护证据；[核验]({p}research_round_545_checks.json)、[条件更新]({p}unified_physics_condition_ledger_545.md)。物理阈值、真实拟合与统一仍开放。'
next_steps='\n\n### 第545轮后：top参数与共同可测条件\n\n545在543同尺度严格裸分支联立完整中性质量、径向质量与谱权重，得到轻质量／Higgs身份／活跃混合之间的精确限制。544运行后的q0=.5例可改变裸权重界，因此不能直接排除运行分支。Schur补测静态逆响应，单场质量还须动能归一化与动量分离。\n\n下一正式编号546优先证明top边界到低尺度值的顺序和唯一性，再将给定参数目标与同一中性／标量读取条件联立。采用声明的参数目标与误差，不能把SM里提取的MS量不经扩展模型匹配当实验排除；不加密无证书网格凑结论。\n\nσ来源、有限匹配、共同几何、维数、Lorentz、实际设备与资源仍开放。目标保持进行中，及时编号报告并保存历史。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第545轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 191.' not in body
        body+='\n\n## 191. 共同真空的质量与可测成分\n\n'+block+next_steps
    else:
        assert '\n## 96.' not in body
        body+='\n\n## 96. 粒子身份须由同一谱权重检验\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为544／2732','最新科学轮次与检查数为545／2739')
        body=body.replace('完成231—544轮。','完成231—545轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第545轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|545|[共同真空的质量与测量成分](research_note_545.md)|[代码](joint_singlet_spectrum_and_weights.py)、'
            '[结果](joint_singlet_spectrum_and_weights_results.json)、[核验](research_round_545_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round545.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=545,cumulative_tests=2739,numbered_scientific_files=947,
    unique_protected_evidence_files=1230,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=546,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
