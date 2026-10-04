"""Publish verified 617, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round617_navigation_checks.json'
checked=core.read(HERE/'research_round_617_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round617_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第617轮完成：** [原规范区域拼接、完整演化与共同能源来源]({p}research_note_617.md)'
         '原紧商群的匹配切分保持原Gauss过程、CAR、热态和几何来源；'
         '非对角电动能须共同分配，重复计能或删剪切均有明确差额。'
         '三组、十二式通过，最新617／3095，1163份编号科学文件、1915份保护证据。'
         '[核验]({p}research_round_617_checks.json)、[条件账]({p}unified_physics_condition_ledger_617.md)。'
         '仅表示切分，非独立区域或真实细化；连续与GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（617后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[618区域拼接、共同作用与几何边界]({p}round618_drafts/STATUS.md)，'
       '先回查旧边界和引力变分结果，再核同一来源的缺口；目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第617轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 263.' not in text
        text+='\n\n## 263. 规范匹配、非对角能源与共同来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 168.' not in text
        text+='\n\n## 168. 表示切分不等于真实空间细化\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—616轮。' in text and '最新科学轮次与检查数为616／3092' in text
        text=text.replace('完成231—616轮。','完成231—617轮。').replace('最新科学轮次与检查数为616／3092','最新科学轮次与检查数为617／3095')
    if p==HERE/'README.md':
        text+='\n\n## 第617轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|617|[原规范区域拼接、完整演化与共同能源来源](research_note_617.md)|[代码](joint_region_energy_gluing.py)、[结果](joint_region_energy_gluing_results.json)、[核验](research_round_617_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round617.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=617,cumulative_tests=3095,numbered_scientific_files=1163,
            unique_protected_evidence_files=1915,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=618,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
