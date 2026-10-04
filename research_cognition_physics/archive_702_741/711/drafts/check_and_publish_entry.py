"""Verify inherited evidence and publish a non-numbered711 entry."""
import ast
import hashlib
import json
import os
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent;RESEARCH=ARCHIVE.parent;ROOT=RESEARCH.parent
sys.path.insert(0,str(ARCHIVE))
import verify_interaction_rounds as core
import topology_sampling_entry as model

result=model.run();assert result==core.read(model.TARGET)
for name,digest in result['dependencies'].items():assert core.digest(ARCHIVE/name)==digest,name
history=dict(core.read(ARCHIVE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for n in range(584,711):
    receipt=core.read(ARCHIVE/f'research_round_{n}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt[key].items():
            assert name not in history or history[name]==digest
            history[name]=digest
history.update(core.read(ARCHIVE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==2943
for name,digest in history.items():assert core.digest(ARCHIVE/name)==digest,name
assert core.text_checks(HERE/'topology_sampling_entry.md')['display_formulas']==3
for name in ('topology_sampling_entry.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
target=HERE/'entry_checks.json';assert not target.exists()
for link in core.link_parser()((HERE/'topology_sampling_entry.md').read_text('utf8')):
    p=(HERE/link).resolve();assert p.exists() or p==target.resolve(),link
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',ARCHIVE/'README.md']
before={p:p.read_bytes() for p in paths};planned={}
summary='**711入口已执行，未计完成轮次：** [有限规范记录与连续拓扑资料]({p}round711_drafts/topology_sampling_entry.md)原完整构形的连通性限制整数标签；显式光滑弱SU(2)映射在全部采样节点和原链路上与身份相同，却有非零三维绕数。此为接口筛查，不当四维瞬子或连续路线反证；正式科学轮次保持710。下一步核实际历史、拼接和同一来源。'
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**711入口已执行，未计完成轮次：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    head,rest=text.split('\n\n',1)
    planned[p]=(head+'\n\n'+summary.format(p=prefix)+'\n\n'+rest).replace('\n',nl).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder=ARCHIVE/'navigation_before_round711_entry_20261003';folder.mkdir(exist_ok=False)
manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round711entry.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
for name,digest in history.items():assert core.digest(ARCHIVE/name)==digest,name
artifact_names=('topology_sampling_entry.py','topology_sampling_entry_results.json',
                'topology_sampling_entry.md','check_and_publish_entry.py')
checks=dict(date='2026-10-03',entry_round=711,latest_formal_round=710,new_formal_round=False,
    saved_results_reproduced=True,previous_round_evidence_unchanged=True,protected_evidence_files=2943,
    dependency_hashes_unchanged=True,navigation_snapshots_preserved=True,navigation_links=links,broken_links=0,
    active_goal_unchanged=True,all_checks_passed=True,
    artifact_hashes={n:core.digest(HERE/n) for n in artifact_names})
with target.open('x',encoding='utf8') as f:json.dump(checks,f,ensure_ascii=False,indent=2)
print(json.dumps({k:checks[k] for k in ('entry_round','latest_formal_round','navigation_links','all_checks_passed')}))
