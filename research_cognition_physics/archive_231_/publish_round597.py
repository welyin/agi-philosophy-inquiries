"""Publish verified 597, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round597_navigation_checks.json'
checked=core.read(HERE/'research_round_597_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round597_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第597轮完成：** [原物质符号记忆、几何响应与内部参照]({p}research_note_597.md)'
         '原势的双符号最低点、全局对称及关系记录的原边能源已接通；全局翻转不改几何过程，局部相对符号可改变来源。'
         '长期保持仍须真实谱隙与读出对比度，准备／装置未完成。三组、九式及主代理核验通过，'
         '最新597／3037，1103份编号科学文件、1764份保护证据。'
         '[核验]({p}research_round_597_checks.json)、[条件账]({p}unified_physics_condition_ledger_597.md)。'
         '本轮是原候选实现与相容性检验，非认知设计的必然性证明；未取得新的独立代理审查。')
order=('**当前执行顺序（597后，优先于下方历史安排）：** 接续'
       '[598原交互的关系记录]({p}round598_drafts/STATUS.md)，'
       '用577旧结果核关系信号和共同来源；接通后返回联合条件总账，不自动扩展成长程装置设计。统一目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第597轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 243.' not in text
        text+='\n\n## 243. 符号记忆与共同几何的适用边界\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 148.' not in text
        text+='\n\n## 148. 关系记录的原来源与准备边界\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—596轮。' in text and '最新科学轮次与检查数为596／3034' in text
        text=text.replace('完成231—596轮。','完成231—597轮。').replace('最新科学轮次与检查数为596／3034','最新科学轮次与检查数为597／3037')
    if p==HERE/'README.md':
        text+='\n\n## 第597轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|597|[原物质符号记忆与共同来源](research_note_597.md)|[代码](joint_singlet_memory_symmetry.py)、[结果](joint_singlet_memory_symmetry_results.json)、[核验](research_round_597_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round597.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=597,cumulative_tests=3037,numbered_scientific_files=1103,
            unique_protected_evidence_files=1764,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=598,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
