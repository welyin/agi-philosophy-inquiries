"""Publish round 533 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import shared_spectral_budget as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round533_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_533_checks.json')
assert checked['round']==533 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round533_20260930'
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
summary=('**第533轮完成：** [同一谱作用的联合系数界]({p}research_note_533.md)'
    '在明确的正非增截止、同一正权及a₄截断下证明V≥6g_w²K²和径向λ≥g_w²/6。'
    '调模块权或复制同类代不能在固定弱耦合下任意压小裸V/M_eff⁴；非单调正截止给必要前提反例。'
    '9组与独立终审通过，最新533／2647，911份编号科学文件，1143份保护证据；'
    '[核验]({p}research_round_533_checks.json)、[条件更新]({p}unified_physics_condition_ledger_533.md)。'
    '这是同尺度裸截断约束，未求出耦合引力真空、观测宇宙学常数或完整统一理论。')
next_steps='\n\n### 第533轮后：联立标量与几何背景\n\n533在同一权重／截止下消去规范和身份迹的共同自由度，得到V≥6g_w²K²及径向λ≥g_w²/6。成熟热核与矩不等式直接继承；非增截止、无独立部门反项和当前有限模块均是明确附加输入。正而非单调截止反例说明不能扩大适用域。\n\n下一正式编号534：在同一K(t)、V(t)及Weyl平方项下求常场／常曲率的耦合驻点，检验533平直径向极小是否保持，以及解是否落入受控低曲率范围。不得用固定平直的Higgs极小替代引力反馈，或把裸V直接当观测Λ；必要时检验明确的新截止／场内容／完整谱作用分支。三维、Lorentz、实际测量、重整化及整体统一仍开放。保留旧研究和冻结证据，及时交付编号报告，不改应用目标、不设置定时任务。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第533轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 179.' not in body
        body+='\n\n## 179. 同一谱作用的规范、标量与真空约束\n\n'+block+next_steps
    else:
        assert '\n## 84.' not in body
        body+='\n\n## 84. 联合谱系数不等于空间或引力真空推导\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为532／2638','最新科学轮次与检查数为533／2647')
        body=body.replace('完成231—532轮。','完成231—533轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第533轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|533|[共同谱系数约束](research_note_533.md)|[代码](shared_spectral_budget.py)、'
            '[结果](shared_spectral_budget_results.json)、[核验](research_round_533_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round533.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=533,cumulative_tests=2647,numbered_scientific_files=911,
    unique_protected_evidence_files=1143,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=534,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
