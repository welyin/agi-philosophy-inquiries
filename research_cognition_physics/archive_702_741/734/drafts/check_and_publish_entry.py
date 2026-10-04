"""Verify inherited evidence and publish a non-numbered734 entry."""
import ast
import hashlib
import json
import os
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent;RESEARCH=ARCHIVE.parent;ROOT=RESEARCH.parent
sys.path.insert(0,str(ARCHIVE))
import verify_interaction_rounds as core
import causal_reference_entry as model

result=json.loads(json.dumps(model.run()));assert result==core.read(model.TARGET)
for name,digest in result['dependencies'].items():assert core.digest(ARCHIVE/name)==digest,name
history=dict(core.read(ARCHIVE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for n in range(584,734):
    receipt=core.read(ARCHIVE/f'research_round_{n}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt[key].items():
            assert name not in history or history[name]==digest
            history[name]=digest
history.update(core.read(ARCHIVE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==3270
for name,digest in history.items():assert core.digest(ARCHIVE/name)==digest,name
assert core.text_checks(HERE/'causal_reference_entry.md')['display_formulas']==3
for name in ('causal_reference_entry.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
target=HERE/'entry_checks.json';assert not target.exists()
for link in core.link_parser()((HERE/'causal_reference_entry.md').read_text('utf8')):
    p=(HERE/link).resolve();assert p.exists() or p==target.resolve(),link
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',ARCHIVE/'README.md']
before={p:p.read_bytes() for p in paths};planned={}
summary='**734入口已执行，未计完成轮次：** [同一过去的因果参考响应]( {p}round734_drafts/causal_reference_entry.md)原完整背景脉冲终点恢复同一H，实际参考和记录仍留下响应；保原记录的新参考另需准备修正。下一项接共同减除与完整来源，正式保持733／3404。'
summary=summary.replace(']( ','](')
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**734入口已执行，未计完成轮次：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    head,rest=text.split('\n\n',1)
    planned[p]=(head+'\n\n'+summary.format(p=prefix)+'\n\n'+rest).replace('\n',nl).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder=ARCHIVE/'navigation_before_round734_entry_20261004';folder.mkdir(exist_ok=False)
manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round734entry.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
for name,digest in history.items():assert core.digest(ARCHIVE/name)==digest,name
artifact_names=('causal_reference_entry.py','causal_reference_entry_results.json',
                'causal_reference_entry.md','check_and_publish_entry.py')
checks=dict(date='2026-10-04',entry_round=734,latest_formal_round=733,new_formal_round=False,
    saved_results_reproduced=True,previous_round_evidence_unchanged=True,protected_evidence_files=3270,
    dependency_hashes_unchanged=True,navigation_snapshots_preserved=True,navigation_links=links,broken_links=0,
    active_goal_unchanged=True,all_checks_passed=True,
    artifact_hashes={n:core.digest(HERE/n) for n in artifact_names})
with target.open('x',encoding='utf8') as f:json.dump(checks,f,ensure_ascii=False,indent=2)
print(json.dumps({k:checks[k] for k in ('entry_round','latest_formal_round','navigation_links','all_checks_passed')}))
