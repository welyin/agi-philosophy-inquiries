"""Publish verified 635, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round635_navigation_checks.json'
checked=core.read(HERE/'research_round_635_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round635_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第635轮完成：** [原物种、曲目标梯度与区域熵的共同几何系数]({p}research_note_635.md)'
         '原五标量和整代费米质量共同匹配面积与R项；'
         '原RS项强制给梯度几何熵，同场值同面积不再由一个常数系数覆盖。'
         '三组、十二式通过，最新635／3146，1217份编号科学文件、2043份保护证据。'
         '[核验]({p}research_round_635_checks.json)、[条件账]({p}unified_physics_condition_ledger_635.md)。'
         '限部分一圈及局部UV；规范边界、有限熵、Newton值和GR仍开放。')
order=('**当前执行顺序（635后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[636原Gauss区域、边界通量与规范熵]({p}round636_drafts/STATUS.md)，'
       '回查617原商群切分及旧熵结果，对接共同态与边界模式；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第635轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 281.' not in text
        text+='\n\n## 281. 原物种和曲目标的共同几何熵\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 186.' not in text
        text+='\n\n## 186. 几何熵匹配不选择维数或Newton常数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—634轮。' in text and '最新科学轮次与检查数为634／3143' in text
        text=text.replace('完成231—634轮。','完成231—635轮。').replace('最新科学轮次与检查数为634／3143','最新科学轮次与检查数为635／3146')
    if p==HERE/'README.md':
        text+='\n\n## 第635轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|635|[原物种、曲目标梯度与区域熵的共同几何系数](research_note_635.md)|[代码](joint_entropy_geometric_matching.py)、[结果](joint_entropy_geometric_matching_results.json)、[核验](research_round_635_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round635.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=635,cumulative_tests=3146,numbered_scientific_files=1217,
            unique_protected_evidence_files=2043,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=636,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
