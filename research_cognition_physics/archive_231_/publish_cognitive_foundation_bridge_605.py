"""Publish non-numbered early-foundation inheritance evidence; preserve science."""
import hashlib,json,os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'cognitive_foundation_bridge_605_navigation.json'
audit=core.read(HERE/'cognitive_foundation_bridge_605_results.json')
assert audit['all_document_checks_passed']
for name,digest in audit['evidence_hashes'].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_foundation_bridge605_20261001'
assert not folder.exists() and not TARGET.exists()
entry=('**001—230轮继承接口复核（不增轮次）：** [认知结构到当前模型的十项映射]({p}cognitive_foundation_bridge_605.md)'
       '区分完整物理系统、有限谱代码、规范区域与经典记录，逐项保留FUCP、组合、操作权限、Time及有限实验的量词。'
       '复Hilbert表示不等于实际合同已实现；早期成果必须接到对象、过程与来源。'
       '[文档审计]({p}cognitive_foundation_bridge_605_results.json)通过，科学基线仍605／3060。'
       '[606补充入口]({p}round606_drafts/foundation_inheritance_entry.md)已登记；先整合条件，认知设计后置，目标不改。')
planned={}
for path,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**001—230轮继承接口复核（不增轮次）：**' not in text
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    if path in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+entry.format(p=prefix)+'\n\n'+rest
    else:text+='\n\n'+entry.format(p=prefix)+'\n'
    planned[path]=text.replace('\n',nl).encode(enc)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
assert all(path.read_bytes()==raw for path,raw in before.items())
folder.mkdir()
manifest={}
for i,(path,raw) in enumerate(before.items()):
    name=f'{i}_{path.name}'
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[path.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for path,raw in planned.items():
    assert path.read_bytes()==before[path]
    tmp=path.with_name(path.name+'.foundation605.tmp')
    with tmp.open('xb') as f:f.write(raw)
    os.replace(tmp,path)
supplemental=('cognitive_foundation_bridge_605.md','cognitive_foundation_bridge_605_audit.py',
'cognitive_foundation_bridge_605_results.json','round606_drafts/foundation_inheritance_entry.md')
result=dict(date='2026-10-01',scientific_baseline=605,scientific_checks=3060,
supplementary_artifact_hashes={n:core.digest(HERE/n) for n in supplemental},
navigation_hashes={p.relative_to(ROOT).as_posix():core.digest(p) for p in paths},
navigation_files=len(paths),local_links_checked=links,broken_links=0,
new_scientific_rounds=0,active_goal_changed=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(science=605,inheritance_rows=10,navigation_files=7,local_links=links,all_passed=True)))
