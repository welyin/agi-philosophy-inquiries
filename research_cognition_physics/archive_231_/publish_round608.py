"""Publish verified 608, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round608_navigation_checks.json'
checked=core.read(HERE/'research_round_608_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round608_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第608轮完成：** [投影补全的独立组合、局域性与共同过程]({p}research_note_608.md)'
         '整体P／Q规则一般不保独立组合，全空间可出现距离无关影响；'
         '实际在位约束的逐项补全保局部支持，并保原P内全部历史、参考和来源。'
         '三组、十二式通过，最新608／3069，1136份编号科学文件、1849份保护证据。'
         '[核验]({p}research_round_608_checks.json)、[条件账]({p}unified_physics_condition_ledger_608.md)。'
         '反例操作出P；实际GW、无界传播及完整物理连接仍开放，无新增独立代理审查。')
order=('**当前执行顺序（608后，优先于下方历史安排）：** 继续由共同过程与组合约束合并条件，认知实现后置。'
       '接[609一粒子协变补全与Fock共同来源]({p}round609_drafts/STATUS.md)，'
       '返回真实投影及原质量接口，不把在位乘积分解输入GW核，目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第608轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 254.' not in text
        text+='\n\n## 254. 独立组合检验约束全局投影补全\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 159.' not in text
        text+='\n\n## 159. 共同过程与局域支持的条件性整合\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—607轮。' in text and '最新科学轮次与检查数为607／3066' in text
        text=text.replace('完成231—607轮。','完成231—608轮。').replace('最新科学轮次与检查数为607／3066','最新科学轮次与检查数为608／3069')
    if p==HERE/'README.md':
        text+='\n\n## 第608轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|608|[投影补全的独立组合、局域性与共同过程](research_note_608.md)|[代码](joint_projection_local_composition.py)、[结果](joint_projection_local_composition_results.json)、[核验](research_round_608_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round608.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=608,cumulative_tests=3069,numbered_scientific_files=1136,
            unique_protected_evidence_files=1849,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=609,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
