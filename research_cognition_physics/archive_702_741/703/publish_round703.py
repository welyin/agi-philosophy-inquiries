"""Publish verified 703, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round703_navigation_checks.json'
checked=core.read(HERE/'research_round_703_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round703_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第703轮完成：** [共同几何来源、热响应与原记录]({p}research_note_703.md)'
         '原完整固定图的几何一二阶热来源、Hessian接触、归一态和立即原记录有同一正近似极限；'
         '固定参照只改表示，C²路径不需新增解析公理。'
         '两组、十六式通过，最新703／3313，1421份编号科学文件、2839份保护证据。'
         '[核验]({p}research_round_703_checks.json)、[全条件账]({p}unified_physics_condition_ledger_703.md)。'
         '含等待的近似来源导数及空间连续仍开放；旧空间与目标不变。')
order=('**当前执行顺序（703后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[704共同热源与真实记录响应]({p}round704_drafts/STATUS.md)，'
       '核同一近似族的真实等待、来源及热规整，不重复623—624已存在的原响应；'
       '保留四分支和699反例范围，停止静态热导数优化。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第703轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 349.' not in text
        text+='\n\n## 349. 原几何来源及热记录共用一个正极限\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 254.' not in text
        text+='\n\n## 254. 旧空间接口继承，几何热来源不替代时空生成\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—702轮。' in text and '最新科学轮次与检查数为702／3311' in text
        text=text.replace('完成231—702轮。','完成231—703轮。').replace('最新科学轮次与检查数为702／3311','最新科学轮次与检查数为703／3313')
    if p==HERE/'README.md':
        text+='\n\n## 第703轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|703|[共同几何来源、热响应与原记录](research_note_703.md)|[代码](joint_geometry_thermal_limit.py)、[结果](joint_geometry_thermal_limit_results.json)、[核验](research_round_703_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round703.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=703,cumulative_tests=3313,numbered_scientific_files=1421,
            unique_protected_evidence_files=2839,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=704,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
