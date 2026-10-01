"""Publish verified 602, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round602_navigation_checks.json'
checked=core.read(HERE/'research_round_602_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round602_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第602轮完成：** [量子来源、参考态与实时响应的共同匹配]({p}research_note_602.md)'
         '原32模式Dirac/Majorana质量部门的静态与隔离响应差额由同一守恒来源涨落固定；'
         'lapse、温度参考、接触项和噪声须共同匹配，独立Fock酉演化已核验。'
         '三组、十式通过，最新602／3051，1118份编号科学文件、1803份保护证据。'
         '[核验]({p}research_round_602_checks.json)、[条件账]({p}unified_physics_condition_ledger_602.md)。'
         '具体计算限冻结物质背景，完整Gauss热态与量子连续未完成；无新增独立代理审查。')
order=('**当前执行顺序（602后，优先于下方历史安排）：** 先整合条件，认知实现后置。'
       '接[603完整非线性Gauss热态与来源域]({p}round603_drafts/STATUS.md)，'
       '核原势是否控制热迹及响应，不把有限质量因子替代完整模型；统一目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第602轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 248.' not in text
        text+='\n\n## 248. 共同量子态与静态及实时来源匹配\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 153.' not in text
        text+='\n\n## 153. 量子来源匹配没有消除维数与几何输入\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—601轮。' in text and '最新科学轮次与检查数为601／3048' in text
        text=text.replace('完成231—601轮。','完成231—602轮。').replace('最新科学轮次与检查数为601／3048','最新科学轮次与检查数为602／3051')
    if p==HERE/'README.md':
        text+='\n\n## 第602轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|602|[量子来源、参考态与实时响应](research_note_602.md)|[代码](joint_quantum_response_matching.py)、[结果](joint_quantum_response_matching_results.json)、[核验](research_round_602_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round602.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=602,cumulative_tests=3051,numbered_scientific_files=1118,
            unique_protected_evidence_files=1803,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=603,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
