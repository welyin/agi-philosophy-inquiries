"""Verify inherited evidence and publish a non-numbered 714 entry."""
import ast
import hashlib
import json
import os
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent;RESEARCH=ARCHIVE.parent;ROOT=RESEARCH.parent
sys.path.insert(0,str(ARCHIVE))
import verify_interaction_rounds as core
import center_reference_entry as model
import center_cross_gradient_check as supplement

result=json.loads(json.dumps(model.run()));assert result==core.read(model.TARGET)
extra=json.loads(json.dumps(supplement.run()));assert extra==core.read(supplement.TARGET)
for data in (result,extra):
    for name,digest in data['dependencies'].items():assert core.digest(ARCHIVE/name)==digest,name
history=dict(core.read(ARCHIVE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for n in range(584,714):
    receipt=core.read(ARCHIVE/f'research_round_{n}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt[key].items():
            assert name not in history or history[name]==digest
            history[name]=digest
history.update(core.read(ARCHIVE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==2986
for name,digest in history.items():assert core.digest(ARCHIVE/name)==digest,name
assert core.text_checks(HERE/'center_reference_entry.md')['display_formulas']==3
for name in ('center_reference_entry.py','center_cross_gradient_check.py',Path(__file__).name):
    ast.parse((HERE/name).read_text('utf8'))
target=HERE/'entry_checks.json';assert not target.exists()
for link in core.link_parser()((HERE/'center_reference_entry.md').read_text('utf8')):
    p=(HERE/link).resolve();assert p.exists() or p==target.resolve(),link
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',ARCHIVE/'README.md']
before={p:p.read_bytes() for p in paths};planned={}
summary='**714入口已执行，未计完成轮次：** [相反记录、相同总能量与参考边界]({p}round714_drafts/center_reference_entry.md)原整面中心变换给指定正常Gauss态相反环路信号、相同完整H均值及平均几何来源，含589非对角梯度；不等于全H对称或免费制备。正式保持713／3343。下一步核原Gibbs参考上的相对资源与带荷物质贡献，不把裸电能下界当状态转换代价。'
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**714入口已执行，未计完成轮次：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    head,rest=text.split('\n\n',1)
    planned[p]=(head+'\n\n'+summary.format(p=prefix)+'\n\n'+rest).replace('\n',nl).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder=ARCHIVE/'navigation_before_round714_entry_20261003';folder.mkdir(exist_ok=False)
manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round714entry.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
for name,digest in history.items():assert core.digest(ARCHIVE/name)==digest,name
for p,raw in planned.items():assert p.read_bytes()==raw,p
artifact_names=('center_reference_entry.py','center_reference_entry_results.json',
                'center_reference_entry.md','center_cross_gradient_check.py',
                'center_cross_gradient_results.json','check_and_publish_entry.py')
checks=dict(date='2026-10-03',entry_round=714,latest_formal_round=713,new_formal_round=False,
    saved_results_reproduced=True,previous_round_evidence_unchanged=True,protected_evidence_files=2986,
    full_cross_gradient_supplement_reproduced=True,
    dependency_hashes_unchanged=True,navigation_snapshots_preserved=True,navigation_links=links,broken_links=0,
    active_goal_unchanged=True,all_checks_passed=True,
    artifact_hashes={n:core.digest(HERE/n) for n in artifact_names})
with target.open('x',encoding='utf8') as f:json.dump(checks,f,ensure_ascii=False,indent=2)
print(json.dumps({k:checks[k] for k in ('entry_round','latest_formal_round','navigation_links','all_checks_passed')}))
