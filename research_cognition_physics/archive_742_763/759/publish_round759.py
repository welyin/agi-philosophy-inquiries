"""Publish759 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round759_navigation_checks.json'
checked=core.read(HERE/'research_round_759_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round759_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第759轮完成（限定反例）：** [确定背景动态闭合的联合来源检验]({p}research_note_759.md)原合法混合态的质量源协方差严格正，排除对全部输入统一保指定联合矩的单确定物质相轨道合同；保多背景、记忆及量子慢变量。三组、十六式通过，最新759／3464，1589份编号科学文件、3613份保护证据。[核验]({p}research_round_759_checks.json)、[条件账]({p}unified_physics_condition_ledger_759.md)。不是对GR或全部宏观近似的反证。'
order='**当前执行顺序（759后，优先于下方历史安排）：** 接[760背景分布与共同过程]({p}round760_drafts/STATUS.md)，让原量子过程决定必要背景分布/相干及来源，先核成熟输运工具的真实前提；不继续反例常数或单项精度。[范围审计]({p}round759_drafts/scope_and_dedup_review.md)。认知共同候选、旧空间、604、649／699及目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第759轮完成（限定反例）：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 405.' not in text
        text+='\n\n## 405. 确定背景动态闭合的联合来源检验\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 310.' not in text
        text+='\n\n## 310. 旧空间合同保持，确定物质背景限制\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—758轮。' in text and '最新科学轮次与检查数为758／3461' in text
        text=text.replace('完成231—758轮。','完成231—759轮。').replace('最新科学轮次与检查数为758／3461','最新科学轮次与检查数为759／3464')
    if p==HERE/'README.md':
        text+='\n\n## 第759轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|759|[确定背景动态闭合的联合来源检验](research_note_759.md)|[代码](joint_background_source_closure.py)、[结果](joint_background_source_closure_results.json)、[核验](research_round_759_checks.json)|\n'
    planned[p]=text.replace('\n',nl).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder.mkdir(exist_ok=False);manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round759.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=759,cumulative_tests=3464,numbered_scientific_files=1589,
    unique_protected_evidence_files=3613,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=760,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
