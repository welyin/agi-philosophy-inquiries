"""Publish verified 590, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round590_navigation_checks.json'
checked=core.read(HERE/'research_round_590_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round590_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第590轮完成：** [共同几何功、记录过程与来源涨落]({p}research_note_590.md)'
         '原物质的含时传播接几何功与实际注能；严格Gauss有限平均能源态仍可有无限来源二阶量。'
         '显式内部量子几何控制分支实现同源能量交换，惯量与势仍为输入。六组、十四式及主代理核验通过，'
         '最新590／3011，1082份编号科学文件、1691份保护证据。'
         '[核验]({p}research_round_590_checks.json)、[条件账]({p}unified_physics_condition_ledger_590.md)。'
         '未得到Einstein动力学、固定ℏ连续或自主装置；未取得新的独立代理审查。')
order=('**当前执行顺序（590后，优先于下方历史安排）：** 接续'
       '[591低能几何来源与引力条件]({p}round591_drafts/STATUS.md)，'
       '先核同一完整物质的低能谱、来源域与有效几何响应，再比较已有引力条件；不继续孤立控制器参数扫描，统一目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第590轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 236.' not in text
        text+='\n\n## 236. 几何功、来源涨落与内部能量账\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 141.' not in text
        text+='\n\n## 141. 内部几何控制仍需物理动力学条件\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—589轮。' in text and '最新科学轮次与检查数为589／3005' in text
        text=text.replace('完成231—589轮。','完成231—590轮。').replace('最新科学轮次与检查数为589／3005','最新科学轮次与检查数为590／3011')
    if p==HERE/'README.md':
        text+='\n\n## 第590轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|590|[共同几何功与来源涨落](research_note_590.md)|[代码](joint_geometry_work_noise.py)、[结果](joint_geometry_work_noise_results.json)、[核验](research_round_590_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round590.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=590,cumulative_tests=3011,numbered_scientific_files=1082,
            unique_protected_evidence_files=1691,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=591,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
