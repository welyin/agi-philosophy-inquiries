"""Publish verified 702, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round702_navigation_checks.json'
checked=core.read(HERE/'research_round_702_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round702_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第702轮完成：** [原完整过程的共同热态、演化与记录]({p}research_note_702.md)'
         '保留原全部相互作用的新正分割有统一热尾控制；自身Gauss Gibbs态、真实时间、'
         '记录联合后态、全部能量矩及热熵共同收敛，补齐655变化热态缺口。'
         '两组、十六式通过，最新702／3311，1418份编号科学文件、2823份保护证据。'
         '[核验]({p}research_round_702_checks.json)、[全条件账]({p}unified_physics_condition_ledger_702.md)。'
         '固定图，不修复699或证明空间连续；旧空间接口与目标不变。')
order=('**当前执行顺序（702后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[703共同热过程的实际来源响应]({p}round703_drafts/STATUS.md)，'
       '区分有界物理来源和改变动能的一般几何来源，核同一分割、态及来源导数；'
       '停止有限格点优化，保留四分支及699反例范围。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第702轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 348.' not in text
        text+='\n\n## 348. 原完整过程的热参考、记录和演化共用一个极限\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 253.' not in text
        text+='\n\n## 253. 旧空间接口继承，热尾紧性不替代空间细化\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—701轮。' in text and '最新科学轮次与检查数为701／3309' in text
        text=text.replace('完成231—701轮。','完成231—702轮。').replace('最新科学轮次与检查数为701／3309','最新科学轮次与检查数为702／3311')
    if p==HERE/'README.md':
        text+='\n\n## 第702轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|702|[原完整过程的共同热态、演化与记录](research_note_702.md)|[代码](joint_gibbs_preserving_transfer.py)、[结果](joint_gibbs_preserving_transfer_results.json)、[核验](research_round_702_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round702.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=702,cumulative_tests=3311,numbered_scientific_files=1418,
            unique_protected_evidence_files=2823,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=703,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
