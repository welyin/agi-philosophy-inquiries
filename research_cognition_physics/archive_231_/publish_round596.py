"""Publish verified 596, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round596_navigation_checks.json'
checked=core.read(HERE/'research_round_596_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round596_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第596轮完成：** [原物质的谱边界与有限精度记录容量]({p}research_note_596.md)'
         '原固定图不能精确全时承载595连续谱硬件；有限能源下，允许读错的记录仍有有限容量界。'
         '原几何模式及诊断编码两组、七式经主代理核验，最新596／3034，1100份编号科学文件、1757份保护证据。'
         '[核验]({p}research_round_596_checks.json)、[条件账]({p}unified_physics_condition_ledger_596.md)。'
         '有限窗口与大系统分支保留；未算原全谱或证明记忆装置，未取得新的独立代理审查。')
order=('**当前执行顺序（596后，优先于下方历史安排）：** 接续'
       '[597原物质的有限窗口记忆]({p}round597_drafts/STATUS.md)，'
       '核原势对称结构、量子保持时间、读写与共同几何来源；复用路径积分接口，统一目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第596轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 242.' not in text
        text+='\n\n## 242. 有限能源记录的谱边界\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 147.' not in text
        text+='\n\n## 147. 谱兼容性不替代空间与引力生成\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—595轮。' in text and '最新科学轮次与检查数为595／3032' in text
        text=text.replace('完成231—595轮。','完成231—596轮。').replace('最新科学轮次与检查数为595／3032','最新科学轮次与检查数为596／3034')
    if p==HERE/'README.md':
        text+='\n\n## 第596轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|596|[原物质谱边界与记录容量](research_note_596.md)|[代码](joint_spectral_record_capacity.py)、[结果](joint_spectral_record_capacity_results.json)、[核验](research_round_596_checks.json)|\n'
    planned[p]=text.replace('\n',newline).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder.mkdir(exist_ok=False)
manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as stream:stream.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as stream:json.dump(manifest,stream,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round596.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=596,cumulative_tests=3034,numbered_scientific_files=1100,
            unique_protected_evidence_files=1757,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=597,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
