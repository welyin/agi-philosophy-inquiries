"""Publish verified 592, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round592_navigation_checks.json'
checked=core.read(HERE/'research_round_592_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round592_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第592轮完成：** [原记录的有限带近似与共同几何来源]({p}research_note_592.md)'
         '原连续仪器无非零有限维精确封闭带；同一初始H²预算经实际读取递推，控制完整记录后态及几何平均来源的共同近似。'
         '仅有迹距离与平均能源不足，已给原模型反例。四组、九式及主代理核验通过，'
         '最新592／3021，1088份编号科学文件、1713份保护证据。'
         '[核验]({p}research_round_592_checks.json)、[条件账]({p}unified_physics_condition_ledger_592.md)。'
         '谱重置为数学比较，非免费物理操作；变化几何、装置与引力约束未完成，未取得新的独立代理审查。')
order=('**当前执行顺序（592后，优先于下方历史安排）：** 接续'
       '[593联合量子几何的记录反作用与约束能源账]({p}round593_drafts/STATUS.md)，'
       '核原读取在内部几何分支中的能源、力和约束来源，落实尚未计入的装置交换；不继续优化截断常数，统一目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第592轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 238.' not in text
        text+='\n\n## 238. 原记录与几何来源的共同有限带近似\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 143.' not in text
        text+='\n\n## 143. 共同来源预算不代替内部装置\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—591轮。' in text and '最新科学轮次与检查数为591／3017' in text
        text=text.replace('完成231—591轮。','完成231—592轮。').replace('最新科学轮次与检查数为591／3017','最新科学轮次与检查数为592／3021')
    if p==HERE/'README.md':
        text+='\n\n## 第592轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|592|[有限带记录与共同来源](research_note_592.md)|[代码](joint_record_band_control.py)、[结果](joint_record_band_control_results.json)、[核验](research_round_592_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round592.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=592,cumulative_tests=3021,numbered_scientific_files=1088,
            unique_protected_evidence_files=1713,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=593,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
