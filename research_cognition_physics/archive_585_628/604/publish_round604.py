"""Publish verified 604, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round604_navigation_checks.json'
checked=core.read(HERE/'research_round_604_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round604_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第604轮完成：** [手征动能与共同量子来源]({p}research_note_604.md)'
         '原质量在明确的中心差分分支中产生八套同质量低能激发，热量、标量力与尺度来源共同复制；'
         '独立连续量子系数不匹配，而旧单项组合对此失明。'
         '三组、十四式通过，最新604／3056，1124份编号科学文件、1817份保护证据。'
         '[核验]({p}research_round_604_checks.json)、[条件账]({p}unified_physics_condition_ledger_604.md)。'
         '仅排除具体接法，未反证其他手征正则化或统一目标；无新增独立代理审查。')
order=('**当前执行顺序（604后，优先于下方历史安排）：** 先整合条件，认知实现后置。'
       '接[605成熟手征正则化与原量子模型的共同条件]({p}round605_drafts/STATUS.md)，'
       '核配置域、量子测度、Gauss与正性接口；不继续节点复制旧实验或主体装置设计，目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第604轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 250.' not in text
        text+='\n\n## 250. 手征物种与共同量子来源的匹配\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 155.' not in text
        text+='\n\n## 155. 手征动能仍需同态同源的连续接口\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—603轮。' in text and '最新科学轮次与检查数为603／3053' in text
        text=text.replace('完成231—603轮。','完成231—604轮。').replace('最新科学轮次与检查数为603／3053','最新科学轮次与检查数为604／3056')
    if p==HERE/'README.md':
        text+='\n\n## 第604轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|604|[手征动能与共同量子来源](research_note_604.md)|[代码](joint_chiral_source_matching.py)、[结果](joint_chiral_source_matching_results.json)、[核验](research_round_604_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round604.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=604,cumulative_tests=3056,numbered_scientific_files=1124,
            unique_protected_evidence_files=1817,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=605,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
