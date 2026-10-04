"""Preserve executed747 native remote-population entry evidence; formal round remains746/3430."""
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
import remote_population_entry as entry
TARGET=HERE/'entry_checks.json'
assert not TARGET.exists(),'Published entry is immutable.'
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',ARCHIVE/'README.md']
before={p:p.read_bytes() for p in paths}
result=entry.run();assert result==core.read(entry.TARGET)
for name,digest in result['dependencies'].items():assert core.digest(ARCHIVE/name)==digest,name
history=dict(core.read(ARCHIVE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for n in range(584,747):
    receipt=core.read(ARCHIVE/f'research_round_{n}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt[key].items():
            assert name not in history or history[name]==digest
            history[name]=digest
history.update(core.read(ARCHIVE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==3449
for name,digest in history.items():assert core.digest(ARCHIVE/name)==digest,name
note=HERE/'research_note_747_working.md';checks=core.text_checks(note)
assert checks['display_formulas']==4
for name in ('remote_population_entry.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
report_links=0
for link in core.link_parser()(note.read_text('utf8')):
    resolved=(HERE/link).resolve();assert resolved.exists() or resolved==TARGET.resolve(),link
    report_links+=1
heading='**747入口已执行，正式仍746／3430：**'
summary=heading+' [原占据转移与实际异地读口]({p}round747_drafts/research_note_747_working.md)原完整H首步使邻端sterile占据在二阶出现严格输入差，原质量背景同时保留；这尚非原T读口的实际记录。下一项核真实异地读口及两条传播通道的共同时间次序，不插入理想中间占据测量。'
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
snapshot=ARCHIVE/'navigation_before_round747_entry_20261004';snapshot.mkdir(exist_ok=False)
manifest={}
for index,(p,raw) in enumerate(before.items()):
    name=f'{index}_{p.name}'
    with (snapshot/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (snapshot/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round747entry.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
assert all(p.read_bytes()==raw for p,raw in planned.items())
for name,digest in history.items():assert core.digest(ARCHIVE/name)==digest,name
artifacts=('remote_population_entry.py','remote_population_entry_results.json','research_note_747_working.md',Path(__file__).name)
receipt=dict(date='2026-10-04',entry_round=747,latest_completed_round=746,
             formal_test_count_unchanged=3430,completed_new_round=False,
             protected_artifacts=3449,previous_hashes_unchanged=True,saved_results_reproduced=True,
             text_checks=checks,local_report_links=report_links,navigation_files=5,
             navigation_links=links,broken_links=0,navigation_snapshots_preserved=True,
             artifact_hashes={name:core.digest(HERE/name) for name in artifacts},
             active_goal_unchanged=True,no_scheduled_task_created=True,no_visual_check=True,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ('entry_round','latest_completed_round','navigation_links','all_checks_passed')}))
