"""Publish the user-requested integration priority without changing the goal."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent
ROOT=RESEARCH.parent
TARGET=HERE/'condition_integration_597_navigation_checks.json'
checked=core.read(HERE/'joint_condition_integration_597_audit_results.json')
assert checked['all_document_checks_passed'] and checked['new_scientific_round_count']==0
for name,value in checked['source_hashes'].items():
    assert core.digest(HERE/name)==value,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_condition_integration_597_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**当前优先：先整合条件，后置认知实现设计。** 按用户最新指示，'
         '[597基线条件整合]({p}joint_condition_compression_update_597.md)已纳入585—597新增连接并逐项覆盖C01—C27；'
         '区分同模型联合实现、表示冗余与尚未连接的部门。'
         '[证据审计]({p}joint_condition_integration_597_audit_results.json)通过，科学基线仍597／3037，本次不新增轮次。')
order=('**当前执行顺序（优先于下方历史安排）：** 先核已有物质／手征表示与当前Gauss量子候选的连接，'
       '再按共同场内容核量子处方、尺度和引力接口；见[更新入口]({p}round598_drafts/condition_integration_status.md)。'
       '原598关系记录候选保留但不继续执行，记忆、控制及组织设计后置；统一目标不变。')
planned={}
for path,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**当前优先：先整合条件，后置认知实现设计。**' not in text
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    head,rest=text.split('\n\n',1)
    text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    if path==HERE/'README.md':
        text+='\n\n## 597后条件整合与执行入口更新（不新增轮次）\n\n[综合表](joint_condition_compression_update_597.md)、[证据审计](joint_condition_integration_597_audit_results.json)、[可复算审计](joint_condition_integration_597_audit.py)、[新入口](round598_drafts/condition_integration_status.md)。\n'
    planned[path]=text.replace('\n',nl).encode(enc)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
assert all(path.read_bytes()==raw for path,raw in before.items())
folder.mkdir(exist_ok=False)
manifest={}
for i,(path,raw) in enumerate(before.items()):
    name=f'{i}_{path.name}'
    with (folder/name).open('xb') as stream:stream.write(raw)
    manifest[path.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as stream:
    json.dump(manifest,stream,ensure_ascii=False,indent=2)
for path,raw in planned.items():
    assert path.read_bytes()==before[path]
    temp=path.with_name(path.name+'.integration597.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,path)
result=dict(date='2026-10-01',kind='non-numbered user priority update',scientific_baseline=597,
            cumulative_scientific_checks_unchanged=3037,navigation_files=5,navigation_links=links,
            broken_links=0,history_snapshots_preserved=True,original_round598_status_unchanged=True,
            active_goal_changed=False,new_tasks_or_automations=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result))
