"""Publish verified 642, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round642_navigation_checks.json'
checked=core.read(HERE/'research_round_642_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round642_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第642轮完成：** [完整规范物质的共同表示、空间消元与几何来源]({p}research_note_642.md)'
         '原完整商群、曲目标及CAR共用树表示；合法Gauss态在同一内部CAR消元后相同，'
         '完整能源和几何来源仍不同，限定了状态单独闭合的空间合并。'
         '三组、十六式通过，最新642／3165，1238份编号科学文件、2097份保护证据。'
         '[核验]({p}research_round_642_checks.json)、[条件账]({p}unified_physics_condition_ledger_642.md)。'
         '见证非Gibbs／低能态；连续统一与GR仍未完成。')
order=('**当前执行顺序（642后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[643共同模型与剩余跨分支接口]({p}round643_drafts/STATUS.md)，'
       '回填全条件账，核固定图、连续手征、参考与动态几何的共同资料；不继续树控制器设计，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第642轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 288.' not in text
        text+='\n\n## 288. 完整物质的空间表示与共同来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 193.' not in text
        text+='\n\n## 193. 完整Gauss表示不等于物质来源可删除\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—641轮。' in text and '最新科学轮次与检查数为641／3162' in text
        text=text.replace('完成231—641轮。','完成231—642轮。').replace('最新科学轮次与检查数为641／3162','最新科学轮次与检查数为642／3165')
    if p==HERE/'README.md':
        text+='\n\n## 第642轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|642|[完整规范物质的共同表示、空间消元与几何来源](research_note_642.md)|[代码](joint_tree_matter_source.py)、[结果](joint_tree_matter_source_results.json)、[核验](research_round_642_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round642.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=642,cumulative_tests=3165,numbered_scientific_files=1238,
            unique_protected_evidence_files=2097,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=643,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
