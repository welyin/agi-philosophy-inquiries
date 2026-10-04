"""Publish verified 665, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round665_navigation_checks.json'
checked=core.read(HERE/'research_round_665_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round665_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第665轮完成：** [原Gauss量子过程、接触作用与辅助历史的共同条件]({p}research_note_665.md)'
         '原固定图的声明接触扩展保共同域、Gauss热态和原记录；辅助提升、归一及空间能源流须同步。'
         '三组、十八式通过，最新665／3237，1307份编号科学文件、2316份保护证据。'
         '[核验]({p}research_round_665_checks.json)、[全条件账]({p}unified_physics_condition_ledger_665.md)。'
         '限保原K的正几何有限图；连续测度及量子GR仍开放。')
order=('**当前执行顺序（665后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[666原CAR到全左手表示的排序与共同测度]({p}round666_drafts/STATUS.md)，'
       '核真实粒子—空穴字典、真空与原辅助对象；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第665轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 311.' not in text
        text+='\n\n## 311. 原完整Gauss过程与接触辅助表示\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 216.' not in text
        text+='\n\n## 216. 同一辅助表示须保真空提升、归一和空间来源\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—664轮。' in text and '最新科学轮次与检查数为664／3234' in text
        text=text.replace('完成231—664轮。','完成231—665轮。').replace('最新科学轮次与检查数为664／3234','最新科学轮次与检查数为665／3237')
    if p==HERE/'README.md':
        text+='\n\n## 第665轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|665|[原Gauss量子过程、接触作用与辅助历史的共同条件](research_note_665.md)|[代码](joint_contact_gauss_history.py)、[结果](joint_contact_gauss_history_results.json)、[核验](research_round_665_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round665.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=665,cumulative_tests=3237,numbered_scientific_files=1307,
            unique_protected_evidence_files=2316,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=666,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
