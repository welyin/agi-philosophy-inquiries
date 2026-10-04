"""Publish verified 641, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round641_navigation_checks.json'
checked=core.read(HERE/'research_round_641_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round641_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第641轮完成：** [同一空间合并过程的非线性记录与几何响应]({p}research_note_641.md)'
         '原二次场的同一时间关联与对易核确定正弦二值读取的两时记录及几何得分；'
         '等时热态替代和删除对易核均给实际记录差异。'
         '两组、十二式通过，最新641／3162，1235份编号科学文件、2089份保护证据。'
         '[核验]({p}research_round_641_checks.json)、[条件账]({p}unified_physics_condition_ledger_641.md)。'
         '字段线性化仍为输入；完整Gauss、连续与GR未完成。')
order=('**当前执行顺序（641后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[642完整非线性物质与Gauss的共同尺度接口]({p}round642_drafts/STATUS.md)，'
       '回查原商群、曲目标和非对角电动能，不继续Gaussian记录精度优化；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第641轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 287.' not in text
        text+='\n\n## 287. 原空间过程的非线性记录与几何响应\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 192.' not in text
        text+='\n\n## 192. 二次记录接口不替代完整Gauss尺度匹配\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—640轮。' in text and '最新科学轮次与检查数为640／3160' in text
        text=text.replace('完成231—640轮。','完成231—641轮。').replace('最新科学轮次与检查数为640／3160','最新科学轮次与检查数为641／3162')
    if p==HERE/'README.md':
        text+='\n\n## 第641轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|641|[同一空间合并过程的非线性记录与几何响应](research_note_641.md)|[代码](joint_blocked_record_history.py)、[结果](joint_blocked_record_history_results.json)、[核验](research_round_641_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round641.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=641,cumulative_tests=3162,numbered_scientific_files=1235,
            unique_protected_evidence_files=2089,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=642,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
