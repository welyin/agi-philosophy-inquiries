"""Publish verified 650, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round650_navigation_checks.json'
checked=core.read(HERE/'research_round_650_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round650_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第650轮完成：** [原引力约束、面积匹配与物质边界通量的联合条件]({p}research_note_650.md)'
         '原物质模有满足全部经典初始约束的小数据解族；两侧面积匹配仍有单侧通量，'
         '完整物质—几何动量拼接使两侧通量相消。'
         '三组、十六式通过，最新650／3188，1262份编号科学文件、2165份保护证据。'
         '[核验]({p}research_round_650_checks.json)、[条件账]({p}unified_physics_condition_ledger_650.md)。'
         '限给定领先Einstein作用、固定类时切口与经典分支；关系区域及全量子统一仍开放。')
order=('**当前执行顺序（650后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[651物质关系区域、动态嵌入与共同边界配对]({p}round651_drafts/STATUS.md)，'
       '核原参考的实际输送及完整边界条件，不继续面积或控制器参数扫描；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第650轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 296.' not in text
        text+='\n\n## 296. 面积匹配与实际物质—引力通量\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 201.' not in text
        text+='\n\n## 201. 区域演化还须保持完整边界配对\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—649轮。' in text and '最新科学轮次与检查数为649／3185' in text
        text=text.replace('完成231—649轮。','完成231—650轮。').replace('最新科学轮次与检查数为649／3185','最新科学轮次与检查数为650／3188')
    if p==HERE/'README.md':
        text+='\n\n## 第650轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|650|[原引力约束、面积匹配与物质边界通量的联合条件](research_note_650.md)|[代码](joint_gravity_boundary_flux.py)、[结果](joint_gravity_boundary_flux_results.json)、[核验](research_round_650_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round650.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=650,cumulative_tests=3188,numbered_scientific_files=1262,
            unique_protected_evidence_files=2165,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=651,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
