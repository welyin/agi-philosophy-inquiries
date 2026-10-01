"""Publish verified 606, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round606_navigation_checks.json'
checked=core.read(HERE/'research_round_606_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round606_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第606轮完成：** [同一过程映射连接手征子空间、电动能与来源]({p}research_note_606.md)'
         '用早期完整过程要求组织接口：同一编码确定联络、正横向项、Fock提升及几何来源；'
         '非恒定手征投影通常不交织原裸电演化，压缩Hamiltonian不能冒充过程等价。'
         '三组、十四式通过，最新606／3063，1130份编号科学文件、1835份保护证据。'
         '[核验]({p}research_round_606_checks.json)、[条件账]({p}unified_physics_condition_ledger_606.md)。'
         'GW核为明确Euclidean诊断；未完成实际权限、手征Hamiltonian或物理连续，无新增独立代理审查。')
order=('**当前执行顺序（606后，优先于下方历史安排）：** 以早期操作结构寻找共同连接，认知实现后置。'
       '接[607同一编码的动力学交织与多部门条件合并]({p}round607_drafts/STATUS.md)，'
       '联查完整记录、参考、演化与来源映射；不继续单独参数扫描，目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第606轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 252.' not in text
        text+='\n\n## 252. 同一过程映射确定物质与来源补项\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 157.' not in text
        text+='\n\n## 157. 等距子空间不自动保持自然演化\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—605轮。' in text and '最新科学轮次与检查数为605／3060' in text
        text=text.replace('完成231—605轮。','完成231—606轮。').replace('最新科学轮次与检查数为605／3060','最新科学轮次与检查数为606／3063')
    if p==HERE/'README.md':
        text+='\n\n## 第606轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|606|[手征子空间、电动能与共同过程映射](research_note_606.md)|[代码](joint_chiral_fibre_source.py)、[结果](joint_chiral_fibre_source_results.json)、[核验](research_round_606_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round606.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=606,cumulative_tests=3063,numbered_scientific_files=1130,
            unique_protected_evidence_files=1835,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=607,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
