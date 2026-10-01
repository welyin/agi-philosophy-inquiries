"""Publish verified 639, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round639_navigation_checks.json'
checked=core.read(HERE/'research_round_639_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round639_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第639轮完成：** [原区域读量、物理体积与尺度统一的相对信息预算]({p}research_note_639.md)'
         '原完整H的区域注能范数由有效物理体积给出，固定支撑可保统一信息上界；'
         '非均匀几何变分还须保读取权重的变化。'
         '三组、十五式通过，最新639／3157，1229份编号科学文件、2074份保护证据。'
         '[核验]({p}research_round_639_checks.json)、[条件账]({p}unified_physics_condition_ledger_639.md)。'
         '集体读取为候选操作；因果实现、跨图态极限及GR仍开放。')
order=('**当前执行顺序（639后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[640实际尺度映射、共同参考与联合观测]({p}round640_drafts/STATUS.md)，'
       '核原群、参考、记录及物质／几何来源能否共同输送；不继续优化剖面或装置，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第639轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 285.' not in text
        text+='\n\n## 285. 原区域读量与物理体积的统一预算\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 190.' not in text
        text+='\n\n## 190. 区域信息上界不替代跨图连续状态\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—638轮。' in text and '最新科学轮次与检查数为638／3154' in text
        text=text.replace('完成231—638轮。','完成231—639轮。').replace('最新科学轮次与检查数为638／3154','最新科学轮次与检查数为639／3157')
    if p==HERE/'README.md':
        text+='\n\n## 第639轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|639|[原区域读量、物理体积与尺度统一的相对信息预算](research_note_639.md)|[代码](joint_regional_read_scale.py)、[结果](joint_regional_read_scale_results.json)、[核验](research_round_639_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round639.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=639,cumulative_tests=3157,numbered_scientific_files=1229,
            unique_protected_evidence_files=2074,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=640,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
