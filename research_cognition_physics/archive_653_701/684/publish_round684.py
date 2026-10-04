"""Publish verified 684, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round684_navigation_checks.json'
checked=core.read(HERE/'research_round_684_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round684_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第684轮完成：** [辅助反射的有限复制障碍与原物理子代数]({p}research_note_684.md)'
         '排除原补偿的一整类有限实复制反射；完整补偿可精确保全部物理来源，同时保留辅助负方向。'
         '完整辅助RP仅为充分路线。两组、十六式通过，最新684／3275，1364份编号科学文件、2567份保护证据。'
         '[核验]({p}research_round_684_checks.json)、[全条件账]({p}unified_physics_condition_ledger_684.md)。'
         '原动态Q0正性、H_F身份、共同连续及量子GR仍开放，旧空间接口复用。')
order=('**当前执行顺序（684后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[685原共同极限中的体权重与物理来源]({p}round685_drafts/STATUS.md)，'
       '停止同类复制设计，核真实权重及来源，保留直接物理子代数路线；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第684轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 330.' not in text
        text+='\n\n## 330. 辅助全代数正性并非原物理必要条件\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 235.' not in text
        text+='\n\n## 235. 旧空间合同复用，回到真正物理子代数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—683轮。' in text and '最新科学轮次与检查数为683／3273' in text
        text=text.replace('完成231—683轮。','完成231—684轮。').replace('最新科学轮次与检查数为683／3273','最新科学轮次与检查数为684／3275')
    if p==HERE/'README.md':
        text+='\n\n## 第684轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|684|[辅助反射的有限复制障碍与原物理子代数](research_note_684.md)|[代码](joint_auxiliary_physical_positivity.py)、[结果](joint_auxiliary_physical_positivity_results.json)、[核验](research_round_684_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round684.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=684,cumulative_tests=3275,numbered_scientific_files=1364,
            unique_protected_evidence_files=2567,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=685,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
