"""Publish verified 613, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round613_navigation_checks.json'
checked=core.read(HERE/'research_round_613_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round613_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第613轮完成：** [同一Hamiltonian的手征、规范与几何来源条件]({p}research_note_613.md)'
         '成熟连续时间候选可保有限热态，但朴素手征筛选破坏原弱代数和U(1)周期；'
         '同谱海能与压力不能仅由单点Λ抵消共同匹配。'
         '三组、十二式通过，最新613／3084，1151份编号科学文件、1884份保护证据。'
         '[核验]({p}research_round_613_checks.json)、[C01—C27更新总账]({p}unified_physics_condition_ledger_613.md)。'
         '仅限指定向量样候选和接法；真实手征、连续及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（613后，优先于下方历史安排）：** 先合并物理成立条件，认知系统设计后置。'
       '接[614原一代物质与成熟手征构造]({p}round614_drafts/STATUS.md)，'
       '核原表示、商群与Yukawa／Majorana能否共用候选输入，逐项记录新增场及来源；'
       '不重复精度优化，目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第613轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 259.' not in text
        text+='\n\n## 259. 同一Hamiltonian、规范及几何来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 164.' not in text
        text+='\n\n## 164. 手征、规范和来源不能分别签收\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—612轮。' in text and '最新科学轮次与检查数为612／3081' in text
        text=text.replace('完成231—612轮。','完成231—613轮。').replace('最新科学轮次与检查数为612／3081','最新科学轮次与检查数为613／3084')
    if p==HERE/'README.md':
        text+='\n\n## 第613轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|613|[同一Hamiltonian的手征、规范与几何来源条件](research_note_613.md)|[代码](joint_hamiltonian_gauge_contract.py)、[结果](joint_hamiltonian_gauge_contract_results.json)、[核验](research_round_613_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round613.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=613,cumulative_tests=3084,numbered_scientific_files=1151,
            unique_protected_evidence_files=1884,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=614,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
