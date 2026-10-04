"""Publish verified 614, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round614_navigation_checks.json'
checked=core.read(HERE/'research_round_614_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round614_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第614轮完成：** [原物质、全局商群与成熟手征表示的共同字典]({p}research_note_614.md)'
         '原16内部表示、Z₆商、32模式CAR质量及来源可共同映入候选容器的子群，'
         '无需扩大物理规范群；若扩大为全Spin(10)，原实singlet Majorana接法失败。'
         '三组、十二式通过，最新614／3087，1154份编号科学文件、1891份保护证据。'
         '[核验]({p}research_round_614_checks.json)、[条件账]({p}unified_physics_condition_ledger_614.md)。'
         '有限质量字典不是完整手征测度或时间重建，局域性、连续及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（614后，优先于下方历史安排）：** 继续整合原群、量子过程和共同来源，认知设计后置。'
       '接[615原子群的Weyl测度与共同来源]({p}round615_drafts/STATUS.md)，'
       '核配置限制、测度变化与同一来源的真实连接，保留成熟候选的局域性和重建边界；目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第614轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 260.' not in text
        text+='\n\n## 260. 原群与原质量的共同表示字典\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 165.' not in text
        text+='\n\n## 165. 表示容器不等于物理规范群扩大\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—613轮。' in text and '最新科学轮次与检查数为613／3084' in text
        text=text.replace('完成231—613轮。','完成231—614轮。').replace('最新科学轮次与检查数为613／3084','最新科学轮次与检查数为614／3087')
    if p==HERE/'README.md':
        text+='\n\n## 第614轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|614|[原物质、全局商群与成熟手征表示的共同字典](research_note_614.md)|[代码](joint_spinor_subgroup_mass.py)、[结果](joint_spinor_subgroup_mass_results.json)、[核验](research_round_614_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round614.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=614,cumulative_tests=3087,numbered_scientific_files=1154,
            unique_protected_evidence_files=1891,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=615,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
