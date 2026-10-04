"""Publish verified 680, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round680_navigation_checks.json'
checked=core.read(HERE/'research_round_680_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round680_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第680轮完成：** [完整物理代数中的质量约化与严格归一边界]({p}research_note_680.md)'
         '固定原H_b下，完整Gauss物理代数的质量插入可逆，全部允许局部质量的RP等价于无质量RP；'
         '精确Grassmann反例保留指定归一的独立门槛。两组、十四式通过，最新680／3267，1352份编号科学文件、2513份保护证据。'
         '[核验]({p}research_round_680_checks.json)、[全条件账]({p}unified_physics_condition_ledger_680.md)。'
         '原无质量动态正性、指定归一、H_F身份及共同连续／量子GR仍开放；旧空间接口复用。')
order=('**当前执行顺序（680后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[681原无质量泛函与动态规范构造]({p}round681_drafts/STATUS.md)，'
       '核成熟构造的实际背景、来源与正时间映射，不再质量调参；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第680轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 326.' not in text
        text+='\n\n## 326. 质量约化不替代原过程与空间接口\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 231.' not in text
        text+='\n\n## 231. 旧方向与热参考合同保持，归一门槛独立\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—679轮。' in text and '最新科学轮次与检查数为679／3265' in text
        text=text.replace('完成231—679轮。','完成231—680轮。').replace('最新科学轮次与检查数为679／3265','最新科学轮次与检查数为680／3267')
    if p==HERE/'README.md':
        text+='\n\n## 第680轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|680|[完整物理代数中的质量约化与严格归一边界](research_note_680.md)|[代码](joint_mass_reflection_congruence.py)、[结果](joint_mass_reflection_congruence_results.json)、[核验](research_round_680_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round680.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=680,cumulative_tests=3267,numbered_scientific_files=1352,
            unique_protected_evidence_files=2513,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=681,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
