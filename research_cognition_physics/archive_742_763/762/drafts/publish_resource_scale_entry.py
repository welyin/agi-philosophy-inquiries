"""Publish the working entry; preserve all historical navigation and science."""
import json
import os
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
ARCHIVE=HERE.parent
RESEARCH=ARCHIVE.parent
ROOT=RESEARCH.parent
sys.path.insert(0,str(ARCHIVE))
import verify_interaction_rounds as core
import verify_resource_scale_entry as check

assert check.verify()==core.read(check.TARGET)
target=HERE/'resource_scale_entry_navigation_checks.json'
snapshot=HERE/'navigation_before_resource_scale_entry'
assert not target.exists() and not snapshot.exists()
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',ARCHIVE/'README.md']
before={p:p.read_bytes() for p in paths}
message='**762入口已核，正式仍761／3470：** [资源尺度与同一背景基线]({p}round762_drafts/research_note_762_working.md)给出新增大作用量分支的精确计数，并核原准备宽度对首阶报告的影响；仅认同参数不能把费米圈当成完整绝对首阶。下一项保原Gauss准备，接玻色协方差、均值修正及761同源响应。[入口核验]({p}round762_drafts/resource_scale_entry_checks.json)。这是工作报告，不增加完成轮次，目标保持。'
planned={}
links=0
for path,raw in before.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(encoding).replace('\r\n','\n')
    assert '**762入口已核，正式仍761／3470：**' not in text
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else \
           'archive_231_/' if path.parent==RESEARCH else ''
    heading,rest=text.split('\n\n',1)
    text=heading+'\n\n'+message.format(p=prefix)+'\n\n'+rest
    for link in core.link_parser()(text):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
    planned[path]=text.replace('\n',newline).encode(encoding)
assert all(p.read_bytes()==raw for p,raw in before.items())
snapshot.mkdir()
manifest={}
for i,(path,raw) in enumerate(before.items()):
    name=f'{i}_{path.name}'
    with (snapshot/name).open('xb') as stream:stream.write(raw)
    manifest[path.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=core.digest(path))
with (snapshot/'manifest.json').open('x',encoding='utf8') as stream:
    json.dump(manifest,stream,ensure_ascii=False,indent=2)
for path,raw in planned.items():
    assert path.read_bytes()==before[path]
    temp=path.with_name(path.name+'.resource_scale_entry.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,path)
assert check.verify()==core.read(check.TARGET)
for name,item in manifest.items():
    assert core.digest(snapshot/item['snapshot'])==item['sha256']
result=dict(date='2026-10-04',kind='round762_working_entry_publication',
    latest_completed_scientific_round=761,cumulative_numbered_tests=3470,
    historical_protected_files=3637,navigation_files=5,navigation_links=links,
    navigation_snapshots_preserved=True,broken_links=0,all_checks_passed=True,
    active_goal_unchanged=True,new_automation=False,stage_complete=False)
with target.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result))

