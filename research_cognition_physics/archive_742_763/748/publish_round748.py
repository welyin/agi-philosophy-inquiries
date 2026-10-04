"""Publish748 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round748_navigation_checks.json'
checked=core.read(HERE/'research_round_748_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round748_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第748轮完成（合成公式）：** [原标量边与异地总系数]({p}research_note_748.md)旧局部密度与原测地边力固定标量通道，和原跳跃共同给总响应；校准实例数值同号约−0.12565，严格总符号尚缺。无需独立通信耦合。两组、十二式通过，最新748／3434，1556份编号科学文件、3480份保护证据。[核验]({p}research_round_748_checks.json)、[条件账]({p}unified_physics_condition_ledger_748.md)。'
order='**当前执行顺序（748后，优先于下方历史安排）：** 接[749共同指认与跨部门合同]({p}round749_drafts/STATUS.md)，先复用区域／尺度／来源旧结果，寻找联合净约束；保留局部总信号证书缺口，不继续特殊读口精度优化。统一目标及旧空间、604、649／699保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第748轮完成（合成公式）：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 394.' not in text
        text+='\n\n## 394. 原标量边与异地总系数\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 299.' not in text
        text+='\n\n## 299. 旧空间合同保持，原总响应接回联合合同\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—747轮。' in text and '最新科学轮次与检查数为747／3432' in text
        text=text.replace('完成231—747轮。','完成231—748轮。').replace('最新科学轮次与检查数为747／3432','最新科学轮次与检查数为748／3434')
    if p==HERE/'README.md':
        text+='\n\n## 第748轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|748|[原标量边与异地总系数](research_note_748.md)|[代码](joint_remote_total_coefficient.py)、[结果](joint_remote_total_coefficient_results.json)、[核验](research_round_748_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round748.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=748,cumulative_tests=3434,numbered_scientific_files=1556,
    unique_protected_evidence_files=3480,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=749,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
