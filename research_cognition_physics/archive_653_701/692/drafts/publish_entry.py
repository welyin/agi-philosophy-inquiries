"""Register692 entry only; keep691 as the latest completed scientific round."""
import hashlib
import json
import os
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
ARCHIVE=HERE.parent
RESEARCH=ARCHIVE.parent
ROOT=RESEARCH.parent
sys.path.insert(0,str(ARCHIVE))
import verify_interaction_rounds as core
check=core.read(HERE/'entry_checks.json')
assert check['all_checks_passed'] and not check['new_formal_round_claimed']
for name,digest in check['artifact_hashes'].items():
    assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',ARCHIVE/'README.md']
snapshot=ARCHIVE/'navigation_before_round692_entry_20261002'
target=HERE/'navigation_checks.json'
assert not snapshot.exists() and not target.exists()
before={p:p.read_bytes() for p in paths}
planned={}
links=0
for p,raw in before.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(encoding).replace('\r\n','\n')
    marker='**692入口已执行，未计完成轮次：**'
    assert marker not in text and '**第691轮完成：**' in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    entry=(marker+' [原奇Gauss复合来源]('+prefix+'round692_drafts/odd_gauss_source_entry.md)'
        '已接入原非平坦完整矩阵；局部规范、反射和中性偶配对核验通过。'
        '完整S9／双Haar／H_b平均的正性仍未判定，继承675—680，停止无控制抽样。'
        '[入口核验]('+prefix+'round692_drafts/entry_checks.json)。最新正式科学轮次仍为691，旧空间接口及目标不变。')
    head,rest=text.split('\n\n',1)
    text=head+'\n\n'+entry+'\n\n'+rest
    for link in core.link_parser()(text):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
    planned[p]=text.replace('\n',newline).encode(encoding)
assert all(p.read_bytes()==raw for p,raw in before.items())
snapshot.mkdir()
manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (snapshot/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (snapshot/'manifest.json').open('x',encoding='utf8') as f:
    json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round692_entry.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
out=dict(date='2026-10-02',entry_round=692,latest_formal_round=691,
    navigation_files=len(paths),navigation_links=links,broken_links=0,
    prior_navigation_snapshots_preserved=True,scientific_counters_unchanged=True,
    goal_unchanged=True,all_checks_passed=True)
with target.open('x',encoding='utf8') as f:json.dump(out,f,ensure_ascii=False,indent=2)
print(json.dumps(out))
