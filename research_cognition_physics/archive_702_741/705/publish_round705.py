"""Publish verified 705, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round705_navigation_checks.json'
checked=core.read(HERE/'research_round_705_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round705_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第705轮完成：** [区域压缩、边界电荷与Gauss支持]({p}research_note_705.md)'
         '原完整热态的无限边界支持使严格单侧有限输出与精确Gauss冲突，给出三误差权衡；'
         '保边界、压缩重数的局部通道可行，但整体仍无限维。'
         '三组、二十式通过，最新705／3318，1427份编号科学文件、2870份保护证据。'
         '[核验]({p}research_round_705_checks.json)、[全条件账]({p}unified_physics_condition_ledger_705.md)。'
         '新族动态来源及空间连续仍开放；入口交换子归属363已补正，旧空间与目标不变。')
order=('**当前执行顺序（705后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[706保边界局部近似与原真实过程]({p}round706_drafts/STATUS.md)，'
       '核比较能源到实际来源的共同域，复用592／623—625／637；'
       '保留四分支及699范围，停止有限容量反例与精度扫描。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第705轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 351.' not in text
        text+='\n\n## 351. 区域、边界支持及有限容量共同验收\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 256.' not in text
        text+='\n\n## 256. 旧空间接口继承，内部边界表示不当空间维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—704轮。' in text and '最新科学轮次与检查数为704／3315' in text
        text=text.replace('完成231—704轮。','完成231—705轮。').replace('最新科学轮次与检查数为704／3315','最新科学轮次与检查数为705／3318')
    if p==HERE/'README.md':
        text+='\n\n## 第705轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|705|[区域压缩、边界电荷与Gauss支持](research_note_705.md)|[代码](joint_regional_charge_compression.py)、[结果](joint_regional_charge_compression_results.json)、[核验](research_round_705_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round705.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=705,cumulative_tests=3318,numbered_scientific_files=1427,
            unique_protected_evidence_files=2870,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=706,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
