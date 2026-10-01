"""Publish verified 660, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round660_navigation_checks.json'
checked=core.read(HERE/'research_round_660_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round660_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第660轮完成：** [原手征权重、局部配对与观测变换的共同表示]({p}research_note_660.md)'
         '同一局部配对精确保留原带相位权重及变换后的观测；自由有限盒候选的反射正性与严格归一成立。'
         '原观测的时间支撑尚未接通，帧Jacobian及来源必须同时输送。'
         '四组、十六式通过，最新660／3225，1292份编号科学文件、2269份保护证据。'
         '[核验]({p}research_round_660_checks.json)、[条件账]({p}unified_physics_condition_ledger_660.md)。'
         '原共同Hamiltonian、一般规范场及量子GR仍开放。')
order=('**当前执行顺序（660后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[661原观测代数与共同物理时间]({p}round661_drafts/STATUS.md)，'
       '核实际观测的时间支撑、反射及演化映射，限定接口后回到其它共同条件；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第660轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 306.' not in text
        text+='\n\n## 306. 原权重与局部正表示的观测接口\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 211.' not in text
        text+='\n\n## 211. 权重相同仍须保留观测与来源映射\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—659轮。' in text and '最新科学轮次与检查数为659／3221' in text
        text=text.replace('完成231—659轮。','完成231—660轮。').replace('最新科学轮次与检查数为659／3221','最新科学轮次与检查数为660／3225')
    if p==HERE/'README.md':
        text+='\n\n## 第660轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|660|[原手征权重、局部配对与观测变换的共同表示](research_note_660.md)|[代码](joint_local_mirror_process.py)、[结果](joint_local_mirror_process_results.json)、[核验](research_round_660_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round660.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=660,cumulative_tests=3225,numbered_scientific_files=1292,
            unique_protected_evidence_files=2269,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=661,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
