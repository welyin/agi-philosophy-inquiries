"""Preserve executed741 source-order evidence; formal round remains740/3419."""
import ast
import hashlib
import json
import os
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
RESEARCH=ARCHIVE.parent;ROOT=RESEARCH.parent
sys.path.insert(0,str(ARCHIVE))
import verify_interaction_rounds as core
import absolute_source_order_entry as entry
TARGET=HERE/'entry_checks.json'
assert not TARGET.exists(),'Published entry is immutable.'
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',ARCHIVE/'README.md']
before={p:p.read_bytes() for p in paths}
result=entry.run();assert result==core.read(entry.TARGET)
for name,digest in result['dependencies'].items():assert core.digest(ARCHIVE/name)==digest,name
history=dict(core.read(ARCHIVE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for n in range(584,741):
    receipt=core.read(ARCHIVE/f'research_round_{n}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt[key].items():
            assert name not in history or history[name]==digest
            history[name]=digest
history.update(core.read(ARCHIVE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==3359
for name,digest in history.items():assert core.digest(ARCHIVE/name)==digest,name
note=HERE/'research_note_741_working.md';checks=core.text_checks(note)
assert checks['display_formulas']==7
for name in ('absolute_source_order_entry.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
report_links=0
for link in core.link_parser()(note.read_text('utf8')):
    resolved=(HERE/link).resolve();assert resolved.exists() or resolved==TARGET.resolve(),link
    report_links+=1
heading='**741入口已执行，正式仍740／3419：**'
summary=heading+' [实际绝对来源与共同圈阶]( {p}round741_drafts/research_note_741_working.md)原同一准备与记录的完整历史响应按统一阶次在二阶进入反作用；首阶只需领先背景上的绝对源。原对象校准通过，下一项接完整初始约束、背景族与真实参考的共同残差。'
summary=summary.replace(']( ','](')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    content=raw.decode(enc).replace('\r\n','\n');assert heading not in content
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    first,rest=content.split('\n\n',1)
    planned[p]=(first+'\n\n'+summary.format(p=prefix)+'\n\n'+rest).replace('\n',nl).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items()),'Concurrent navigation edit'
snapshot=ARCHIVE/'navigation_before_round741_entry_20261004';snapshot.mkdir(exist_ok=False)
manifest={}
for index,(p,raw) in enumerate(before.items()):
    name=f'{index}_{p.name}'
    with (snapshot/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (snapshot/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round741entry.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
assert all(p.read_bytes()==raw for p,raw in planned.items())
for name,digest in history.items():assert core.digest(ARCHIVE/name)==digest,name
artifacts=('absolute_source_order_entry.py','absolute_source_order_entry_results.json','research_note_741_working.md',Path(__file__).name)
receipt=dict(date='2026-10-04',entry_round=741,latest_completed_round=740,
             formal_test_count_unchanged=3419,completed_new_round=False,
             protected_artifacts=3359,previous_hashes_unchanged=True,saved_results_reproduced=True,
             text_checks=checks,local_report_links=report_links,navigation_files=5,
             navigation_links=links,broken_links=0,navigation_snapshots_preserved=True,
             artifact_hashes={name:core.digest(HERE/name) for name in artifacts},
             active_goal_unchanged=True,no_scheduled_task_created=True,no_visual_check=True,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ('entry_round','latest_completed_round','navigation_links','all_checks_passed')}))
