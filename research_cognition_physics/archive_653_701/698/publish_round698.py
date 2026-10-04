"""Publish verified 698, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round698_navigation_checks.json'
checked=core.read(HERE/'research_round_698_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round698_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第698轮完成：** [原静态中心项的完整积分严格上界]({p}research_note_698.md)'
         '整个原群和两个S9平均后|B11|≤U<2.46×10⁻¹⁰，整数常数项已认证。'
         '负核判据仍缺B10的严格下界，未宣称完整RP失败。'
         '两组、十八式通过，最新698／3303，1406份编号科学文件、2761份保护证据。'
         '[核验]({p}research_round_698_checks.json)、[全条件账]({p}unified_physics_condition_ledger_698.md)。'
         '旧空间接口直接复用，目标及四分支不变。')
order=('**当前执行顺序（698后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[699原非对角完整积分的下界]({p}round699_drafts/STATUS.md)，'
       '核B10≥727/10¹⁶是否成立，保留原群、球面与H_b；'
       '数值估计不代替证书，旧空间合同继续复用。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第698轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 344.' not in text
        text+='\n\n## 344. 原静态项的完整积分上界及剩余非对角条件\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 249.' not in text
        text+='\n\n## 249. 旧空间合同复用，积分界不新增空间条件\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—697轮。' in text and '最新科学轮次与检查数为697／3301' in text
        text=text.replace('完成231—697轮。','完成231—698轮。').replace('最新科学轮次与检查数为697／3301','最新科学轮次与检查数为698／3303')
    if p==HERE/'README.md':
        text+='\n\n## 第698轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|698|[原静态中心项的完整积分严格上界](research_note_698.md)|[代码](joint_static_haar_majorant.py)、[结果](joint_static_haar_majorant_results.json)、[核验](research_round_698_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round698.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=698,cumulative_tests=3303,numbered_scientific_files=1406,
            unique_protected_evidence_files=2761,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=699,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
