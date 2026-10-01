"""Publish verified 603, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round603_navigation_checks.json'
checked=core.read(HERE/'research_round_603_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round603_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第603轮完成：** [完整非线性Gauss物质的热态与有限零频来源]({p}research_note_603.md)'
         '原势尾、正动能与有限CAR质量共同保证固定图每个正β的Gauss热态及全部能量矩；'
         '形式来源的守恒涨落与Abel零频响应有限，602匹配身份已接入完整非线性模型。'
         '两组、十二式通过，最新603／3053，1121份编号科学文件、1810份保护证据。'
         '[核验]({p}research_round_603_checks.json)、[条件账]({p}unified_physics_condition_ledger_603.md)。'
         '温度、准备与有限截止仍输入；未证明瞬时噪声或量子连续，无新增独立代理审查。')
order=('**当前执行顺序（603后，优先于下方历史安排）：** 先整合条件，认知实现后置。'
       '接[604同一量子来源与连续有效理论的尺度匹配]({p}round604_drafts/STATUS.md)，'
       '保持相同背景、态、源及误差窗口，不把热态存在等同连续统一；目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第603轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 249.' not in text
        text+='\n\n## 249. 完整Gauss热态与零频来源的共同存在\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 154.' not in text
        text+='\n\n## 154. 全Gauss热来源仍保留空间与引力输入\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—602轮。' in text and '最新科学轮次与检查数为602／3051' in text
        text=text.replace('完成231—602轮。','完成231—603轮。').replace('最新科学轮次与检查数为602／3051','最新科学轮次与检查数为603／3053')
    if p==HERE/'README.md':
        text+='\n\n## 第603轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|603|[完整Gauss热态与有限零频来源](research_note_603.md)|[代码](joint_thermal_gauss_source.py)、[结果](joint_thermal_gauss_source_results.json)、[核验](research_round_603_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round603.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=603,cumulative_tests=3053,numbered_scientific_files=1121,
            unique_protected_evidence_files=1810,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=604,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
