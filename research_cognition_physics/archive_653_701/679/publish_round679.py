"""Publish verified 679, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round679_navigation_checks.json'
checked=core.read(HERE/'research_round_679_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round679_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第679轮完成：** [原物理来源的反射、接触抵消与共同内积]({p}research_note_679.md)'
         '方程项变换与来源接触抵消给全部原物理Weyl来源反射，含有限软调节、奇异质量及零权重；'
         '完整Gauss物理型为Hermitian。两组、十六式通过，最新679／3265，1349份编号科学文件、2501份保护证据。'
         '[核验]({p}research_round_679_checks.json)、[全条件账]({p}unified_physics_condition_ledger_679.md)。'
         '正性、非零归一、原H_F时间及共同连续／量子GR仍开放；旧空间接口逐项复用。')
order=('**当前执行顺序（679后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[680完整质量插入与正性约化]({p}round680_drafts/STATUS.md)，'
       '回查已有半区插入，核完整Gauss物理代数的可逆合同；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第679轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 325.' not in text
        text+='\n\n## 325. 来源反射闭合与旧空间条件的逐项复用\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 230.' not in text
        text+='\n\n## 230. 不恢复旧下界已消去条件，不叠加替代坐标路线\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—678轮。' in text and '最新科学轮次与检查数为678／3263' in text
        text=text.replace('完成231—678轮。','完成231—679轮。').replace('最新科学轮次与检查数为678／3263','最新科学轮次与检查数为679／3265')
    if p==HERE/'README.md':
        text+='\n\n## 第679轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|679|[原物理来源的反射、接触抵消与共同内积](research_note_679.md)|[代码](joint_physical_source_reflection.py)、[结果](joint_physical_source_reflection_results.json)、[核验](research_round_679_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round679.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=679,cumulative_tests=3265,numbered_scientific_files=1349,
            unique_protected_evidence_files=2501,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=680,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
