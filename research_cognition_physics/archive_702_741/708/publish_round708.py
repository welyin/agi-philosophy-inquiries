"""Publish708 only after scientific verification; preserve navigation snapshots."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round708_navigation_checks.json'
checked=core.read(HERE/'research_round_708_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round708_20261003'
assert not folder.exists() and not TARGET.exists()
summary=('**第708轮完成：** [继续认知所需的共同涨落资料]({p}research_note_708.md)'
         '原同一P₂同时约束集体读口预算、动能主系数及Majorana差分耦合；'
         '仅保中点／幅度不足，加两矩也不自动完成未来闭合。'
         '三组、十七式通过，最新708／3326，1436份编号科学文件、2915份保护证据。'
         '[核验]({p}research_round_708_checks.json)、[全条件账]({p}unified_physics_condition_ledger_708.md)。'
         '结论限于声明任务和原候选，旧空间与统一目标不变。')
order=('**当前执行顺序（708后，优先于下方历史安排）：** 按2026-10-03用户确认，'
       '用少量认知假说共同约束缺口，具体工程设计后置。'
       '接[709内部条件态与有效记忆]({p}round709_drafts/STATUS.md)，'
       '比较完整保留、有限任务和保记忆方案，核原参考、实际预算及质量响应能否同源；'
       '停止逐矩反例，保留四分支及699范围。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第708轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 354.' not in text
        text+='\n\n## 354. 原内部涨落连接读口预算与质量差分\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 259.' not in text
        text+='\n\n## 259. 旧空间接口保留，认知任务约束不冒充维数生成\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—707轮。' in text and '最新科学轮次与检查数为707／3323' in text
        text=text.replace('完成231—707轮。','完成231—708轮。').replace('最新科学轮次与检查数为707／3323','最新科学轮次与检查数为708／3326')
    if p==HERE/'README.md':
        text+='\n\n## 第708轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|708|[继续认知所需的共同涨落资料](research_note_708.md)|[代码](joint_block_fluctuation_contract.py)、[结果](joint_block_fluctuation_contract_results.json)、[核验](research_round_708_checks.json)|\n'
    planned[p]=text.replace('\n',nl).encode(enc)
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
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round708.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-03',latest_round=708,cumulative_tests=3326,numbered_scientific_files=1436,
    unique_protected_evidence_files=2915,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=709,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
