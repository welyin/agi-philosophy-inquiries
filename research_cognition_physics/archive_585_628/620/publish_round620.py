"""Publish verified 620, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round620_navigation_checks.json'
checked=core.read(HERE/'research_round_620_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round620_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第620轮完成：** [共同量子来源的交换约束、混合涨落与框架匹配]({p}research_note_620.md)'
         '原32模式同态来源变换与实际能量交换共同固定交叉关联；'
         '删掉混合噪声会破坏框架方差和交换身份。'
         '三组、十四式通过，最新620／3104，1172份编号科学文件、1937份保护证据。'
         '[核验]({p}research_round_620_checks.json)、[条件账]({p}unified_physics_condition_ledger_620.md)。'
         '限条件质量部门及明示Ward前提；全Gauss连续来源、动态几何与GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（620后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[621全Gauss过程、涂抹来源与联合响应域]({p}round621_drafts/STATUS.md)，'
       '复用有限能源／零频旧界，核实际联合来源域与连续接口；目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第620轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 266.' not in text
        text+='\n\n## 266. 共同来源、混合涨落与交换约束\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 171.' not in text
        text+='\n\n## 171. 联合量子来源条件不选择维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—619轮。' in text and '最新科学轮次与检查数为619／3101' in text
        text=text.replace('完成231—619轮。','完成231—620轮。').replace('最新科学轮次与检查数为619／3101','最新科学轮次与检查数为620／3104')
    if p==HERE/'README.md':
        text+='\n\n## 第620轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|620|[共同量子来源的交换约束、混合涨落与框架匹配](research_note_620.md)|[代码](joint_quantum_exchange_noise.py)、[结果](joint_quantum_exchange_noise_results.json)、[核验](research_round_620_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round620.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=620,cumulative_tests=3104,numbered_scientific_files=1172,
            unique_protected_evidence_files=1937,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=621,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
