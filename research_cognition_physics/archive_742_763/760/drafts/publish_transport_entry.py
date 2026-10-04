"""Preserve and index the unnumbered760 applicability review."""
import hashlib
import json
import os
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent;RESEARCH=ARCHIVE.parent;ROOT=RESEARCH.parent
sys.path.insert(0,str(ARCHIVE))
import verify_interaction_rounds as core
import transport_source_entry as entry
assert entry.run()==core.read(entry.TARGET)
target=HERE/'transport_entry_checks.json';assert not target.exists()
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',ARCHIVE/'README.md']
before={p:p.read_bytes() for p in paths};planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第759轮完成（限定反例）：**' in text and '**760输运入口已核：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    head,rest=text.split('\n\n',1)
    item='**760输运入口已核：** [保背景分支时的完整来源](%sround760_drafts/transport_applicability_working_report.md)原753周期轻子来源有非零跨能源块部分，不能直接套只处理分块观测的输运定理；接保相干状态近似及同一来源误差。正式仍759／3464，不重报旧谱压缩反例，目标保持。'%prefix
    planned[p]=(head+'\n\n'+item+'\n\n'+rest).replace('\n',nl).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
review=HERE/'transport_applicability_working_report.md';core.text_checks(review)
for link in core.link_parser()(review.read_text('utf8')):assert (HERE/link).resolve().exists(),link
assert all(p.read_bytes()==raw for p,raw in before.items())
folder=ARCHIVE/'navigation_before_round760_entry_20261004';folder.mkdir(exist_ok=False)
manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    tmp=p.with_name(p.name+'.round760-entry.tmp')
    with tmp.open('xb') as f:f.write(raw)
    os.replace(tmp,p)
history=dict(core.read(ARCHIVE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for n in range(584,760):
    receipt=core.read(ARCHIVE/f'research_round_{n}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):history.update(receipt[key])
history.update(core.read(ARCHIVE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==3613
for name,digest in history.items():assert core.digest(ARCHIVE/name)==digest,name
for p,raw in planned.items():assert p.read_bytes()==raw,p
artifacts=('transport_applicability_working_report.md','transport_source_entry.py','transport_source_entry_results.json')
report=dict(latest_scientific_round=759,cumulative_numbered_tests=3464,new_numbered_tests=0,
    protected_evidence_unchanged=len(history),navigation_files=len(paths),navigation_links_checked=links,
    broken_links=0,artifact_hashes={n:core.digest(HERE/n) for n in artifacts},
    navigation_snapshot=folder.name,active_goal_unchanged=True,all_checks_passed=True)
with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
