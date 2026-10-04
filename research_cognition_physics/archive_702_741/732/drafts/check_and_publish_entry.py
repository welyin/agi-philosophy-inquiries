"""Verify inherited evidence and publish a non-numbered732 entry."""
import ast
import hashlib
import json
import os
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent;RESEARCH=ARCHIVE.parent;ROOT=RESEARCH.parent
sys.path.insert(0,str(ARCHIVE))
import verify_interaction_rounds as core
import source_feedback_entry as model

result=json.loads(json.dumps(model.run()));assert result==core.read(model.TARGET)
for name,digest in result['dependencies'].items():assert core.digest(ARCHIVE/name)==digest,name
history=dict(core.read(ARCHIVE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for n in range(584,732):
    receipt=core.read(ARCHIVE/f'research_round_{n}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt[key].items():
            assert name not in history or history[name]==digest
            history[name]=digest
history.update(core.read(ARCHIVE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==3242
for name,digest in history.items():assert core.digest(ARCHIVE/name)==digest,name
assert core.text_checks(HERE/'source_feedback_entry.md')['display_formulas']==3
for name in ('source_feedback_entry.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
target=HERE/'entry_checks.json';assert not target.exists()
for link in core.link_parser()((HERE/'source_feedback_entry.md').read_text('utf8')):
    p=(HERE/link).resolve();assert p.exists() or p==target.resolve(),link
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',ARCHIVE/'README.md']
before={p:p.read_bytes() for p in paths};planned={}
summary='**732入口已执行，未计完成轮次：** [实际来源与准备能量响应]({p}round732_drafts/source_feedback_entry.md)原完整动态参考的几何／径向来源交叉响应不对称，不能直接换成准备能量的普通Hessian。下一项接共同一阶发展与因果来源，正式保持731／3398。'
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**732入口已执行，未计完成轮次：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    head,rest=text.split('\n\n',1)
    planned[p]=(head+'\n\n'+summary.format(p=prefix)+'\n\n'+rest).replace('\n',nl).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder=ARCHIVE/'navigation_before_round732_entry_20261004';folder.mkdir(exist_ok=False)
manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round732entry.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
for name,digest in history.items():assert core.digest(ARCHIVE/name)==digest,name
artifact_names=('source_feedback_entry.py','source_feedback_entry_results.json',
                'source_feedback_entry.md','check_and_publish_entry.py')
checks=dict(date='2026-10-04',entry_round=732,latest_formal_round=731,new_formal_round=False,
    saved_results_reproduced=True,previous_round_evidence_unchanged=True,protected_evidence_files=3242,
    dependency_hashes_unchanged=True,navigation_snapshots_preserved=True,navigation_links=links,broken_links=0,
    active_goal_unchanged=True,all_checks_passed=True,
    artifact_hashes={n:core.digest(HERE/n) for n in artifact_names})
with target.open('x',encoding='utf8') as f:json.dump(checks,f,ensure_ascii=False,indent=2)
print(json.dumps({k:checks[k] for k in ('entry_round','latest_formal_round','navigation_links','all_checks_passed')}))
