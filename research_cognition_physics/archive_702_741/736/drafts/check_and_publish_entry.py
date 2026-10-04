"""Preserve736's executed regularity entry; formal count stays735/3409."""
import ast
import hashlib
import json
import os
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ARCHIVE=HERE.parent
RESEARCH=ARCHIVE.parent
ROOT=RESEARCH.parent
sys.path.insert(0,str(ARCHIVE))
import verify_interaction_rounds as core
import response_regularity_entry as entry
TARGET=HERE/'entry_checks.json'
assert not TARGET.exists(), 'Published entry is immutable.'
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',ARCHIVE/'README.md']
before={p:p.read_bytes() for p in paths}
result=entry.run()
assert result==core.read(entry.TARGET)
for name,digest in result['dependencies'].items():
    assert core.digest(ARCHIVE/name)==digest,name
history=dict(core.read(ARCHIVE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for number in range(584,736):
    receipt=core.read(ARCHIVE/f'research_round_{number}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt[key].items():
            assert name not in history or history[name]==digest
            history[name]=digest
history.update(core.read(ARCHIVE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==3304
for name,digest in history.items():
    assert core.digest(ARCHIVE/name)==digest,name
note=HERE/'response_regularity_entry.md'
text_check=core.text_checks(note)
assert text_check['display_formulas']==5
for name in ('response_regularity_entry.py',Path(__file__).name):
    ast.parse((HERE/name).read_text('utf8'))
report_links=0
for link in core.link_parser()(note.read_text('utf8')):
    resolved=(HERE/link).resolve()
    assert resolved.exists() or resolved==TARGET.resolve(),link
    report_links+=1
heading='**736入口已执行，正式仍735／3409：**'
summary=(heading+' [原有谱与反馈正则性]( {p}round736_drafts/response_regularity_entry.md)'
         '复用630—631完整物质谱，原比较分支的响应含二阶／四阶乘对数；'
         '普通恰好两阶／四阶Sobolev界不成立，局部有限项不能消掉该对数。'
         '此为估计路线限制，不是自洽解不存在；下一项核因果逆核与共同有效阶。').replace(']( ','](')
planned={}
for path,raw in before.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    content=raw.decode(encoding).replace('\r\n','\n')
    assert heading not in content
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    first,rest=content.split('\n\n',1)
    planned[path]=(first+'\n\n'+summary.format(p=prefix)+'\n\n'+rest).replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items()),'Concurrent navigation edit'
snapshot=ARCHIVE/'navigation_before_round736_entry_20261004'
snapshot.mkdir(exist_ok=False)
manifest={}
for index,(path,raw) in enumerate(before.items()):
    name=f'{index}_{path.name}'
    with (snapshot/name).open('xb') as handle:handle.write(raw)
    manifest[path.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (snapshot/'manifest.json').open('x',encoding='utf8') as handle:
    json.dump(manifest,handle,ensure_ascii=False,indent=2)
for path,raw in planned.items():
    assert path.read_bytes()==before[path]
    temp=path.with_name(path.name+'.round736entry.tmp')
    with temp.open('xb') as handle:handle.write(raw)
    os.replace(temp,path)
assert all(p.read_bytes()==raw for p,raw in planned.items())
for name,digest in history.items():
    assert core.digest(ARCHIVE/name)==digest,name
artifacts=('response_regularity_entry.py','response_regularity_entry_results.json',
           'response_regularity_entry.md',Path(__file__).name)
receipt=dict(date='2026-10-04',entry_round=736,latest_completed_round=735,
             formal_test_count_unchanged=3409,completed_new_round=False,
             protected_artifacts=3304,previous_hashes_unchanged=True,saved_results_reproduced=True,
             text_checks=text_check,local_report_links=report_links,navigation_files=5,
             navigation_links=links,broken_links=0,navigation_snapshots_preserved=True,
             artifact_hashes={name:core.digest(HERE/name) for name in artifacts},
             active_goal_unchanged=True,no_scheduled_task_created=True,no_visual_check=True,
             all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as handle:
    handle.write(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ('entry_round','latest_completed_round','navigation_links','all_checks_passed')}))
