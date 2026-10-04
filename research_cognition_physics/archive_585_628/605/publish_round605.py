"""Publish verified 605, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round605_navigation_checks.json'
checked=core.read(HERE/'research_round_605_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round605_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第605轮完成：** [手征配置限制、Gauss热态与共同来源域]({p}research_note_605.md)'
         '原Gauss热态不自动满足硬允许域；规范不变的空间零延拓形式可保热迹与读取来源，'
         '软惩罚给受控热态极限，锐筛选与时间单面截断分别存在能源和正性障碍。'
         '四组、十四式通过，最新605／3060，1127份编号科学文件、1824份保护证据。'
         '[核验]({p}research_round_605_checks.json)、[条件账]({p}unified_physics_condition_ledger_605.md)。'
         '空间限制未完成四维手征测度或物理连续；无新增独立代理审查。')
order=('**当前执行顺序（605后，优先于下方历史安排）：** 先整合条件，认知实现后置。'
       '接[606变化手征子空间与共同来源]({p}round606_drafts/STATUS.md)，'
       '核场依赖投影与原电动能；保留早期认知对象到当前模型的映射缺口，目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第605轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 251.' not in text
        text+='\n\n## 251. 配置域、Gauss热态与来源的共同条件\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 156.' not in text
        text+='\n\n## 156. 空间配置限制不等同四维手征完成\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—604轮。' in text and '最新科学轮次与检查数为604／3056' in text
        text=text.replace('完成231—604轮。','完成231—605轮。').replace('最新科学轮次与检查数为604／3056','最新科学轮次与检查数为605／3060')
    if p==HERE/'README.md':
        text+='\n\n## 第605轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|605|[手征配置限制、Gauss热态与共同来源域](research_note_605.md)|[代码](joint_admissible_gauss_domain.py)、[结果](joint_admissible_gauss_domain_results.json)、[核验](research_round_605_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round605.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=605,cumulative_tests=3060,numbered_scientific_files=1127,
            unique_protected_evidence_files=1824,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=606,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
