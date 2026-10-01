"""Publish verified 649, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round649_navigation_checks.json'
checked=core.read(HERE/'research_round_649_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round649_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第649轮完成：** [面积匹配的物理内积、动力学下降与共同量子几何]({p}research_note_649.md)'
         '原面积匹配有局部正物理内积，但复制的原几何动能不保持零范数类；'
         '共享唯一原量子几何可扩充规范切口并保持完整原过程。'
         '三组、十四式通过，最新649／3185，1259份编号科学文件、2155份保护证据。'
         '[核验]({p}research_round_649_checks.json)、[条件账]({p}unified_physics_condition_ledger_649.md)、'
         '[648式14勘误]({p}round649_drafts/round648_formula14_erratum.md)。'
         '限声明的内部几何模型；完整引力约束、连续映射及统一目标仍开放。')
order=('**当前执行顺序（649后，优先于下方历史安排）：** 继续整合共同条件，认知设计后置。'
       '接[650真实引力约束、动态嵌入与区域生成元]({p}round650_drafts/STATUS.md)，'
       '返回原引力作用的实际边界条件，不以控制模型失败代替GR检验，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第649轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 295.' not in text
        text+='\n\n## 295. 面积约化与同一量子几何的共同过程\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 200.' not in text
        text+='\n\n## 200. 正物理内积还须容纳原演化\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—648轮。' in text and '最新科学轮次与检查数为648／3182' in text
        text=text.replace('完成231—648轮。','完成231—649轮。').replace('最新科学轮次与检查数为648／3182','最新科学轮次与检查数为649／3185')
    if p==HERE/'README.md':
        text+='\n\n## 第649轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|649|[面积匹配的物理内积、动力学下降与共同量子几何](research_note_649.md)|[代码](joint_area_reduction_process.py)、[结果](joint_area_reduction_process_results.json)、[核验](research_round_649_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round649.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=649,cumulative_tests=3185,numbered_scientific_files=1259,
            unique_protected_evidence_files=2155,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=650,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
