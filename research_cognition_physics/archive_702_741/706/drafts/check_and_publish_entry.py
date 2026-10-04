"""Reproduce706 graded local Stinespring and two-sided comparison-energy domain, protect705, and register entry without a new round."""
import hashlib
import json
import os
from pathlib import Path
import local_domain_entry as probe

HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
RESEARCH=ARCHIVE.parent;ROOT=RESEARCH.parent
import sys
sys.path.insert(0,str(ARCHIVE))
import verify_interaction_rounds as core
TARGET=HERE/'entry_checks.json'
assert not TARGET.exists()
assert core.text_checks(HERE/'local_domain_entry.md')['display_formulas']==6
saved=core.read(probe.TARGET);assert probe.run()==saved
for name,digest in saved['dependencies'].items():assert core.digest(ARCHIVE/name)==digest
receipt=core.read(ARCHIVE/'research_round_705_checks.json')
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in receipt[key].items():assert core.digest(ARCHIVE/name)==digest
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',ARCHIVE/'README.md']
folder=ARCHIVE/'navigation_before_round706_entry_20261002'
assert not folder.exists()
before={p:p.read_bytes() for p in paths};planned={};links=0
marker='**706入口已执行，未计完成轮次：**'
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert marker not in text and '**第705轮完成：**' in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    entry=(marker+' [保边界通道的双向能量域]('+prefix+'round706_drafts/local_domain_entry.md)'
        '以固定辅助能量标签接通705局部通道的正向及伴随图范数，623原共同域直接复用。'
        '[705的CAR宇称分块补充]('+prefix+'round706_drafts/parity_completion_705.md)'
        '已核验，旧冻结文件保留；完整真实来源历史下一步检验。'
        '正式科学轮次保持705，目标不变。')
    head,rest=text.split('\n\n',1);text=head+'\n\n'+entry+'\n\n'+rest
    for link in core.link_parser()(text):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
    planned[p]=text.replace('\n',newline).encode(enc)
for link in core.link_parser()((HERE/'local_domain_entry.md').read_text('utf8')):
    assert (HERE/link).resolve().exists() or (HERE/link).resolve()==TARGET
assert all(p.read_bytes()==raw for p,raw in before.items())
folder.mkdir();manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    tmp=p.with_name(p.name+'.round706_entry.tmp')
    with tmp.open('xb') as f:f.write(raw)
    os.replace(tmp,p)
names=('local_domain_entry.py','local_domain_entry_results.json',
       'local_domain_entry.md','check_and_publish_entry.py','prepare_entry_publication.py','parity_completion_705.md')
out=dict(date='2026-10-02',entry_round=706,latest_formal_round=705,new_formal_round=False,
    saved_results_reproduced=True,previous_round_evidence_unchanged=True,dependency_hashes_unchanged=True,
    navigation_snapshots_preserved=True,navigation_links=links,broken_links=0,
    active_goal_unchanged=True,all_checks_passed=True,
    artifact_hashes={name:core.digest(HERE/name) for name in names})
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(out,ensure_ascii=False))
