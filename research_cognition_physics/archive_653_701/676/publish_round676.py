"""Publish verified 676, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round676_navigation_checks.json'
checked=core.read(HERE/'research_round_676_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round676_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第676轮完成：** [非平坦边界初始空间与原物理泛函的零权重支撑]({p}research_note_676.md)'
         '原整数通量精确给72／56负谱分区，659固定64／64边界像误差始终为1；'
         '同一分区同时使辅助权重及全部纯物理来源为零。'
         '因此是指定表示障碍，不是完整物理边缘反例；第五方向正性保留。'
         '两组、十二式通过，最新676／3259，1340份编号科学文件、2460份保护证据。'
         '[核验]({p}research_round_676_checks.json)、[全条件账]({p}unified_physics_condition_ledger_676.md)。'
         '限当前两传播方向盒，完整物理时间、连续及量子GR仍开放。')
order=('**当前执行顺序（676后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[677全谱有理符号与共同正则表示]({p}round677_drafts/STATUS.md)，'
       '保留原观测、零支撑及完整平均，不把第五方向当物理时间；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第676轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 322.' not in text
        text+='\n\n## 322. 非平坦边界表示与原物理来源支撑\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 227.' not in text
        text+='\n\n## 227. 固定初始空间失效与物理边缘反例的区别\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—675轮。' in text and '最新科学轮次与检查数为675／3257' in text
        text=text.replace('完成231—675轮。','完成231—676轮。').replace('最新科学轮次与检查数为675／3257','最新科学轮次与检查数为676／3259')
    if p==HERE/'README.md':
        text+='\n\n## 第676轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|676|[非平坦边界初始空间与原物理泛函的零权重支撑](research_note_676.md)|[代码](joint_nonflat_seed_support.py)、[结果](joint_nonflat_seed_support_results.json)、[核验](research_round_676_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round676.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=676,cumulative_tests=3259,numbered_scientific_files=1340,
            unique_protected_evidence_files=2460,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=677,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
