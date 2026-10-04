"""Publish verified 634, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round634_navigation_checks.json'
checked=core.read(HERE/'research_round_634_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round634_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第634轮完成：** [记录条件化、守恒来源与几何约束的共同边界]({p}research_note_634.md)'
         '原连续记录的概率共同固定分支能源／动量及交叉噪声；'
         '无补偿地拼接前后来源产生非零表面散度，临时接触项不能抵消净荷。'
         '三组、十四式通过，最新634／3143，1214份编号科学文件、2036份保护证据。'
         '[核验]({p}research_round_634_checks.json)、[条件账]({p}unified_physics_condition_ledger_634.md)。'
         '仅排除指定接法；真实总过程、图映射及GR仍开放。')
order=('**当前执行顺序（634后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[635区域熵、原物种与共同几何有效项]({p}round635_drafts/STATUS.md)，'
       '先回查旧面积熵与边界结果，核同一处方的真实连接；不继续记录硬件设计，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第634轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 280.' not in text
        text+='\n\n## 280. 原记录四动量与几何守恒的共同边界\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 185.' not in text
        text+='\n\n## 185. 记录来源匹配不完成空间或GR推导\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—633轮。' in text and '最新科学轮次与检查数为633／3140' in text
        text=text.replace('完成231—633轮。','完成231—634轮。').replace('最新科学轮次与检查数为633／3140','最新科学轮次与检查数为634／3143')
    if p==HERE/'README.md':
        text+='\n\n## 第634轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|634|[记录条件化、守恒来源与几何约束的共同边界](research_note_634.md)|[代码](joint_record_ward_boundary.py)、[结果](joint_record_ward_boundary_results.json)、[核验](research_round_634_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round634.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=634,cumulative_tests=3143,numbered_scientific_files=1214,
            unique_protected_evidence_files=2036,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=635,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
