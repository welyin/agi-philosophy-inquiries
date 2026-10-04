"""Index unnumbered common-object review, preserving prior navigation."""
import hashlib
import json
import os
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent;RESEARCH=ARCHIVE.parent;ROOT=RESEARCH.parent
sys.path.insert(0,str(ARCHIVE))
import verify_interaction_rounds as core
import common_instrument_certificate_check as review
assert review.verify()==core.read(review.TARGET)
target=HERE/'common_instrument_entry_navigation_checks.json'
folder=ARCHIVE/'navigation_before_757_common_instrument_entry_20261004'
assert not target.exists() and not folder.exists()
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',ARCHIVE/'README.md']
before={p:p.read_bytes() for p in paths};planned={}
entry='**757入口已执行，正式仍756／3455：** [同一条件过程证书]({p}round757_drafts/common_instrument_certificate.md)复用634／730／750，已核同一中性记录的概率、条件后态、完整来源及条件线性响应。未把旧CP/全方差公式重计轮次；下一项接原实际相互作用、记录形成与共同来源。[复核]({p}round757_drafts/common_instrument_certificate_checks.json)。当前尚非内部自主仪器或全图连续映射，目标保持。'
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**757入口已执行，正式仍756／3455：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    head,rest=text.split('\n\n',1)
    text=head+'\n\n'+entry.format(p=prefix)+'\n\n'+rest
    planned[p]=text.replace('\n',nl).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder.mkdir();manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.757entry.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
assert all(p.read_bytes()==raw for p,raw in planned.items())
for name,entry in manifest.items():assert core.digest(folder/entry['snapshot'])==entry['sha256']
result=dict(date='2026-10-04',latest_scientific_round=756,cumulative_numbered_tests=3455,
    numbered_scientific_files=1580,protected_historical_files=3574,
    new_numbered_scientific_round=False,navigation_files=5,local_links_checked=links,broken_links=0,
    original_navigation_snapshots_preserved=True,working_report_indexed=True,active_goal_unchanged=True,
    all_checks_passed=True)
with target.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
print(json.dumps(result))
