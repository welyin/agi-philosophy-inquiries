"""Publish verified 636, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round636_navigation_checks.json'
checked=core.read(HERE/'research_round_636_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round636_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第636轮完成：** [原商群的联合边界记录、区域熵与电几何来源]({p}research_note_636.md)'
         '原Z₆边界标签独立化产生13/18禁标签权；'
         '原非对角电来源同均值同协方差的两个合法态仍有不同区域熵。'
         '三组、十四式通过，最新636／3149，1220份编号科学文件、2050份保护证据。'
         '[核验]({p}research_round_636_checks.json)、[条件账]({p}unified_physics_condition_ledger_636.md)。'
         '限原固定图的明示态及电来源；全热态熵、连续匹配和GR仍开放。')
order=('**当前执行顺序（636后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[637原完整热态、区域能源与熵]({p}round637_drafts/STATUS.md)，'
       '复用603完整H与617切分，核共同态的区域熵条件；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第636轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 282.' not in text
        text+='\n\n## 282. 原商群边界记录与电来源的共同限制\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 187.' not in text
        text+='\n\n## 187. 边界熵约束不等于连续几何已涌现\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—635轮。' in text and '最新科学轮次与检查数为635／3146' in text
        text=text.replace('完成231—635轮。','完成231—636轮。').replace('最新科学轮次与检查数为635／3146','最新科学轮次与检查数为636／3149')
    if p==HERE/'README.md':
        text+='\n\n## 第636轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|636|[原商群的联合边界记录、区域熵与电几何来源](research_note_636.md)|[代码](joint_quotient_boundary_entropy.py)、[结果](joint_quotient_boundary_entropy_results.json)、[核验](research_round_636_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round636.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=636,cumulative_tests=3149,numbered_scientific_files=1220,
            unique_protected_evidence_files=2050,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=637,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
