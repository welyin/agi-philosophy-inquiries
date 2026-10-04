"""Publish verified 588, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round588_navigation_checks.json'
checked=core.read(HERE/'research_round_588_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round588_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第588轮完成：** [局部组合、有限生成元障碍与几何补全]({p}research_note_588.md)'
         '同一双lapse符号接原切向来源，但原有限分块菜单产生不可吸收的跨节点动量项；'
         '酉表示扩展保完整过程，完整空间度规补单一ψ遗漏的剪切。五组、十四式及主代理核验通过，'
         '最新588／2999，1076份编号科学文件、1668份保护证据。'
         '[核验]({p}research_round_588_checks.json)、[条件账]({p}unified_physics_condition_ledger_588.md)。'
         '未建立全量子几何、固定ℏ连续或全部HDA；未取得新的独立代理审查。')
order=('**当前执行顺序（588后，优先于下方历史安排）：** 接续'
       '[589完整度规与同一有限量子物质]({p}round589_drafts/STATUS.md)，'
       '核一般正度规的非对角项、Gauss、正性、真实记录及原共形分支；不再重复菜单不闭合反例，统一目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第588轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 234.' not in text
        text+='\n\n## 234. 局部组合的条件与完整几何接口\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 139.' not in text
        text+='\n\n## 139. 单一共形几何不足以容纳一般切向作用\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—587轮。' in text and '最新科学轮次与检查数为587／2994' in text
        text=text.replace('完成231—587轮。','完成231—588轮。').replace('最新科学轮次与检查数为587／2994','最新科学轮次与检查数为588／2999')
    if p==HERE/'README.md':
        text+='\n\n## 第588轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|588|[局部组合与几何补全](research_note_588.md)|[代码](joint_local_geometry_generators.py)、[结果](joint_local_geometry_generators_results.json)、[核验](research_round_588_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round588.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=588,cumulative_tests=2999,numbered_scientific_files=1076,
            unique_protected_evidence_files=1668,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=589,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
