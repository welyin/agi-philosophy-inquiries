"""Publish verified 630, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round630_navigation_checks.json'
checked=core.read(HERE/'research_round_630_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round630_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第630轮完成：** [原物质谱、真实时间响应与几何来源的共同连接]({p}research_note_630.md)'
         '原整代径向质量的正实时谱复现两个旧UV系数，并共同固定真空噪声与共形几何交叉来源；'
         '最低成对阈值限制同一局部展开窗口。'
         '三组、十九式通过，最新630／3130，1202份编号科学文件、2008份保护证据。'
         '[核验]({p}research_round_630_checks.json)、[全账]({p}unified_physics_condition_ledger_630.md)。'
         '限明确连续真空分支；一般张量应力、图映射、实际态和GR仍开放。')
order=('**当前执行顺序（630后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[631完整物质与张量几何响应]({p}round631_drafts/STATUS.md)，'
       '核其它度规通道能否与同一物质、态及旧几何作用共用响应；不继续脉冲精度优化，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第630轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 276.' not in text
        text+='\n\n## 276. 原物质谱与几何来源的共同连接\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 181.' not in text
        text+='\n\n## 181. 共形响应不是一般度规或空间维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—629轮。' in text and '最新科学轮次与检查数为629／3127' in text
        text=text.replace('完成231—629轮。','完成231—630轮。').replace('最新科学轮次与检查数为629／3127','最新科学轮次与检查数为630／3130')
    if p==HERE/'README.md':
        text+='\n\n## 第630轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|630|[原物质谱、真实时间响应与几何来源的共同连接](research_note_630.md)|[代码](joint_continuum_source_spectrum.py)、[结果](joint_continuum_source_spectrum_results.json)、[核验](research_round_630_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round630.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=630,cumulative_tests=3130,numbered_scientific_files=1202,
            unique_protected_evidence_files=2008,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=631,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
