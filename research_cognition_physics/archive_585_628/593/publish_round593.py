"""Publish verified 593, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round593_navigation_checks.json'
checked=core.read(HERE/'research_round_593_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round593_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第593轮完成：** [联合量子几何的实际记录与封闭约束能源账]({p}research_note_593.md)'
         '原读取保持非选择几何边缘，却改变同一总能源与几何力；联合H²预算仍成立。'
         '给C_g＋H_m的条件性仪器障碍、能源保持替换误差及重复记录装置账。三组、九式及主代理核验通过，'
         '最新593／3024，1091份编号科学文件、1724份保护证据。'
         '[核验]({p}research_round_593_checks.json)、[条件账]({p}unified_physics_condition_ledger_593.md)。'
         '未构造GR约束或物理态，未实现自主装置；未取得新的独立代理审查。')
order=('**当前执行顺序（593后，优先于下方历史安排）：** 接续'
       '[594内部装置、能量交换与可保存指针]({p}round594_drafts/STATUS.md)，'
       '核联合实现及末端记录是否真正保住能源账，避免把相互作用或指针相干当成免费资源；统一目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第593轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 239.' not in text
        text+='\n\n## 239. 联合几何记录与闭合约束能源账\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 144.' not in text
        text+='\n\n## 144. 记录装置必须进入同一约束账\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—592轮。' in text and '最新科学轮次与检查数为592／3021' in text
        text=text.replace('完成231—592轮。','完成231—593轮。').replace('最新科学轮次与检查数为592／3021','最新科学轮次与检查数为593／3024')
    if p==HERE/'README.md':
        text+='\n\n## 第593轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|593|[联合记录与封闭能源账](research_note_593.md)|[代码](joint_record_constraint_balance.py)、[结果](joint_record_constraint_balance_results.json)、[核验](research_round_593_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round593.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=593,cumulative_tests=3024,numbered_scientific_files=1091,
            unique_protected_evidence_files=1724,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=594,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
