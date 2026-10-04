"""Publish verified 696, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round696_navigation_checks.json'
checked=core.read(HERE/'research_round_696_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round696_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第696轮完成：** [完整辅助平均后的固定接缝负核]({p}research_note_696.md)'
         '原中心空间历史的全部S9积分精确求出；固定时间接缝的2×2标量核仍有严格负方向。'
         '两组、十六式通过，最新696／3299，1400份编号科学文件、2720份保护证据。'
         '[核验]({p}research_round_696_checks.json)、[全条件账]({p}unified_physics_condition_ledger_696.md)。'
         '时间Haar与H_b尚未纳入，不构成原物理RP或H_F反例；单空间自环不替代实际空间传播。')
order=('**当前执行顺序（696后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[697中心边界的完整规范闭合]({p}round697_drafts/STATUS.md)，'
       '复用673同步输送，保留最终全群holonomy及原H_b；'
       '旧空间接口逐项复用，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第696轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 342.' not in text
        text+='\n\n## 342. 完整辅助平均与尚待完成的规范闭合\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 247.' not in text
        text+='\n\n## 247. 旧空间合同复用，单空间自环不等于空间生成\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—695轮。' in text and '最新科学轮次与检查数为695／3297' in text
        text=text.replace('完成231—695轮。','完成231—696轮。').replace('最新科学轮次与检查数为695／3297','最新科学轮次与检查数为696／3299')
    if p==HERE/'README.md':
        text+='\n\n## 第696轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|696|[完整辅助平均后的固定接缝负核](research_note_696.md)|[代码](joint_dynamic_auxiliary_integral.py)、[结果](joint_dynamic_auxiliary_integral_results.json)、[核验](research_round_696_checks.json)|\n'
    planned[p]=text.replace('\n',newline).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder.mkdir(exist_ok=False)
manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as stream:stream.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as stream:json.dump(manifest,stream,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round696.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=696,cumulative_tests=3299,numbered_scientific_files=1400,
            unique_protected_evidence_files=2720,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=697,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
