"""Publish round 536 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_yukawa_higgs_matching as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round536_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_536_checks.json')
assert checked['round']==536 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round536_20260930'
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
summary='**第536轮完成：** [规范、Yukawa与Higgs共同边界]({p}research_note_536.md)将531同一权重与533规范化联立，证明未知中微子Dirac奇异值补全二／四次矩的精确充要条件并给显式构造。历史规范基准下top＋三ν族有λ≥0.17199279，且构造λ=0.2的共同边界见证；均不是低能实验拟合。8组与独立终审通过，最新536／2671，920份编号科学文件，1164份保护证据；[核验]({p}research_round_536_checks.json)、[条件更新]({p}unified_physics_condition_ledger_536.md)。活跃重场、阈值、seesaw、运行及整体统一仍开放。'
next_steps='\n\n### 第536轮后：将共同边界接到中微子阈值与尺度运行\n\n536给出规范、物理Yukawa和Higgs维四边界的联合可行域。三规范耦合固定同一r、T，已知带电矩K₂、K₄后未知n_g个ν奇异值能补全，当且仅当Rν=T−K₂≥0且K₄+Rν²/(n_g r)≤λT≤K₄+Rν²/r。该充要性不包含seesaw、扰动性、阈值、实际观测或完整稳定真空；原始有限D参数不能直接作为物理Yukawa。\n\n下一正式编号537：先检验同一候选是否同时容许活跃ν_R、轻中微子与536要求的Dirac矩；优先利用成熟type-I seesaw与可能的复相位抵消，不能把无抵消估计当一般反证。阈值之后才可选择一致的RG轨道，将低能输入同共同边界比较。若树级条件已经相冲，先给限定清楚的排除／替代构造，再决定是否进行长程运行；不任意增加未记账的物质。\n\n531冻结规范一圈基准继续作为历史诊断，规范单态ν_R不改变该一圈规范系数，但会影响Yukawa及Higgs运行。1204.0328的λ_src=6λ，三代和全部活跃中微子范围必须保留。四维、Lorentz、参考、几何反馈和观测接口仍开放，不把人工高能λ=0.13与低能数直接比较。每个实质单元继续完整编号报告与条件账；不改应用目标、不设定时任务。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第536轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 182.' not in body
        body+='\n\n## 182. 同一规范与物质边界的可行域\n\n'+block+next_steps
    else:
        assert '\n## 87.' not in body
        body+='\n\n## 87. 共同有效边界仍不选择几何维数\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为535／2663','最新科学轮次与检查数为536／2671')
        body=body.replace('完成231—535轮。','完成231—536轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第536轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|536|[规范与物质共同边界](research_note_536.md)|[代码](joint_yukawa_higgs_matching.py)、'
            '[结果](joint_yukawa_higgs_matching_results.json)、[核验](research_round_536_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round536.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=536,cumulative_tests=2671,numbered_scientific_files=920,
    unique_protected_evidence_files=1164,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=537,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
