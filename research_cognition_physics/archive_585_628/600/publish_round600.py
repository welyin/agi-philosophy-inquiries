"""Publish verified 600, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round600_navigation_checks.json'
checked=core.read(HERE/'research_round_600_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round600_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第600轮完成：** [完整物质应力、Yukawa力与共同有效作用约化]({p}research_note_600.md)'
         '将非零应力迹、标量质量力和质量函数变换接入同一离壳EFT映射，'
         '原无迹规范捷径在费米物质加入后必须补交叉项；原Dirac与Majorana变化已复算。'
         '三组、十二式核验通过，最新600／3045，1112份编号科学文件、1789份保护证据。'
         '[核验]({p}research_round_600_checks.json)、[条件账]({p}unified_physics_condition_ledger_600.md)。'
         '限给定Einstein领先作用的一阶EFT；完整量子来源与约束未完成，未取得独立代理审查。')
order=('**当前执行顺序（600后，优先于下方历史安排）：** 继续先整合条件，后置认知实现设计。'
       '接[601同一EFT阶的约束与物理分支]({p}round601_drafts/STATUS.md)，'
       '核lapse、shift来源、完整物质方程及约束传播，保留低能范围和旧反例；统一目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第600轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 246.' not in text
        text+='\n\n## 246. 完整物质应力与共同有效作用约化\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 151.' not in text
        text+='\n\n## 151. 完整物质约化仍依赖给定的几何与维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—599轮。' in text and '最新科学轮次与检查数为599／3042' in text
        text=text.replace('完成231—599轮。','完成231—600轮。').replace('最新科学轮次与检查数为599／3042','最新科学轮次与检查数为600／3045')
    if p==HERE/'README.md':
        text+='\n\n## 第600轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|600|[完整物质应力与Yukawa力的共同约化](research_note_600.md)|[代码](joint_full_matter_operator_reduction.py)、[结果](joint_full_matter_operator_reduction_results.json)、[核验](research_round_600_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round600.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=600,cumulative_tests=3045,numbered_scientific_files=1112,
            unique_protected_evidence_files=1789,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=601,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
