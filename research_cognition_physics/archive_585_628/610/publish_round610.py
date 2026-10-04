"""Publish verified 610, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round610_navigation_checks.json'
checked=core.read(HERE/'research_round_610_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round610_20261101'
assert not folder.exists() and not TARGET.exists()
summary=('**第610轮完成：** [共同谱隙、局部链路扰动与Fock来源系数]({p}research_note_610.md)'
         '同一谱隙与局部有限秩导数同时控制投影、联络、正能源及Fock尾部；'
         '并核远端改场依赖，界不乘总费米占据数。'
         '三组、十二式通过，最新610／3075，1142份编号科学文件、1863份保护证据。'
         '[核验]({p}research_round_610_checks.json)、[条件账]({p}unified_physics_condition_ledger_610.md)。'
         '硬域非原全热态，系数衰减非无界传播；真实质量、测度及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（610后，优先于下方历史安排）：** 继续整合同一物理对象，认知实现后置。'
       '接[611原物种质量与共同投影的真实字典]({p}round611_drafts/STATUS.md)，'
       '核原32模式与候选旋量投影的质量、规范和来源条件，不继续单独优化衰减率，目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第610轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 256.' not in text
        text+='\n\n## 256. 谱隙、链路导数与共同局部系数\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 161.' not in text
        text+='\n\n## 161. 局部系数不替代无界电动能传播\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—609轮。' in text and '最新科学轮次与检查数为609／3072' in text
        text=text.replace('完成231—609轮。','完成231—610轮。').replace('最新科学轮次与检查数为609／3072','最新科学轮次与检查数为610／3075')
    if p==HERE/'README.md':
        text+='\n\n## 第610轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|610|[共同谱隙、局部链路扰动与Fock来源系数](research_note_610.md)|[代码](joint_gapped_link_locality.py)、[结果](joint_gapped_link_locality_results.json)、[核验](research_round_610_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round610.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=610,cumulative_tests=3075,numbered_scientific_files=1142,
            unique_protected_evidence_files=1863,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=611,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
