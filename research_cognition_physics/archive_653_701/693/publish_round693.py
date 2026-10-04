"""Publish verified 693, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round693_navigation_checks.json'
checked=core.read(HERE/'research_round_693_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round693_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第693轮完成：** [完整规范平均的近零谱与原来源误差]({p}research_note_693.md)'
         '给出对全部原空间规范背景统一的近零谱Haar幂界，接成原完整物理来源及九点调节表达的误差式。'
         '两组、十六式通过，最新693／3293，1391份编号科学文件、2668份保护证据。'
         '[核验]({p}research_round_693_checks.json)、[全条件账]({p}unified_physics_condition_ledger_693.md)。'
         '固定有限盒，界很松、热矩系数未求值；不构成实际符号、RP或连续时空证书。')
order=('**当前执行顺序（693后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[694临界谱附近的原物理来源]({p}round694_drafts/STATUS.md)，'
       '核实际来源的消去或跳变，不做裸谱常数优化；旧空间382—386、425、522—523逐项复用，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第693轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 339.' not in text
        text+='\n\n## 339. 完整规范平均的近零谱与来源控制\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 244.' not in text
        text+='\n\n## 244. 旧空间合同复用，有限谱模量不是空间尺度\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—692轮。' in text and '最新科学轮次与检查数为692／3291' in text
        text=text.replace('完成231—692轮。','完成231—693轮。').replace('最新科学轮次与检查数为692／3291','最新科学轮次与检查数为693／3293')
    if p==HERE/'README.md':
        text+='\n\n## 第693轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|693|[完整规范平均的近零谱与原来源误差](research_note_693.md)|[代码](joint_gauge_sublevel_control.py)、[结果](joint_gauge_sublevel_control_results.json)、[核验](research_round_693_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round693.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=693,cumulative_tests=3293,numbered_scientific_files=1391,
            unique_protected_evidence_files=2668,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=694,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
