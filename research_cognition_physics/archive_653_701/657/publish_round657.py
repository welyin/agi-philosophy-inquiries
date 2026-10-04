"""Publish verified 657, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round657_navigation_checks.json'
checked=core.read(HERE/'research_round_657_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round657_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第657轮完成：** [原空间辅助测度的反射内积与完整球面积分]({p}research_note_657.md)'
         '原自由偶数时间盒的带相位辅助Pfaffian给同一BCS反射Gram核；全S⁹积分为显式平方范数。'
         '两时间片严格归一，任意偶数长度未归一正性成立。'
         '四组、二十一式通过，最新657／3214，1283份编号科学文件、2239份保护证据。'
         '[核验]({p}research_round_657_checks.json)、[条件账]({p}unified_physics_condition_ledger_657.md)。'
         '一般长度严格归一、共同单步过程、物理CAR与量子GR仍开放。')
order=('**当前执行顺序（657后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[658严格归一与跨时间共同过程]({p}round658_drafts/STATUS.md)，'
       '核原平均向量支撑及不同时间盒的实际映射；不把每盒正内积当同一演化，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第657轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 303.' not in text
        text+='\n\n## 303. 原空间测度的反射内积\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 208.' not in text
        text+='\n\n## 208. 原测度、反射内积与完整球面平均须共用对象\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—656轮。' in text and '最新科学轮次与检查数为656／3210' in text
        text=text.replace('完成231—656轮。','完成231—657轮。').replace('最新科学轮次与检查数为656／3210','最新科学轮次与检查数为657／3214')
    if p==HERE/'README.md':
        text+='\n\n## 第657轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|657|[原空间辅助测度的反射内积与完整球面积分](research_note_657.md)|[代码](joint_auxiliary_reflection_gluing.py)、[结果](joint_auxiliary_reflection_gluing_results.json)、[核验](research_round_657_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round657.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=657,cumulative_tests=3214,numbered_scientific_files=1283,
            unique_protected_evidence_files=2239,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=658,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
