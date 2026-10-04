"""Publish verified 598, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round598_navigation_checks.json'
checked=core.read(HERE/'research_round_598_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round598_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第598轮完成：** [费米物质、Gauss约束与共同记录来源]({p}research_note_598.md)'
         '已将旧物质表示接入原有限图量子过程：原束缚势控制新增质量耦合，'
         '合并H半有界、自伴并保Gauss；原读取的注能与平均来源身份保留。'
         '三组、十二式核验通过，最新598／3040，1106份编号科学文件、1775份保护证据。'
         '[核验]({p}research_round_598_checks.json)、[条件账]({p}unified_physics_condition_ledger_598.md)。'
         'CAR、物质谱及跳跃仍明示输入；手征连续、共同传播和GR约束未完成。未取得独立代理审查。')
order=('**当前执行顺序（598后，优先于下方历史安排）：** 继续先整合物理条件，后置认知实现设计。'
       '接[599共同量子处方与尺度]({p}round599_drafts/STATUS.md)，核同一场内容、正则化、'
       '来源和尺度接口，再检验引力连接；不恢复旧598关系记忆支线，统一目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第598轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 244.' not in text
        text+='\n\n## 244. 费米物质与共同量子来源的有限截止连接\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 149.' not in text
        text+='\n\n## 149. 有限费米合并不等于维数或手征连续的推导\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—597轮。' in text and '最新科学轮次与检查数为597／3037' in text
        text=text.replace('完成231—597轮。','完成231—598轮。').replace('最新科学轮次与检查数为597／3037','最新科学轮次与检查数为598／3040')
    if p==HERE/'README.md':
        text+='\n\n## 第598轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|598|[费米物质、Gauss与共同记录来源](research_note_598.md)|[代码](joint_fermion_gauss_completion.py)、[结果](joint_fermion_gauss_completion_results.json)、[核验](research_round_598_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round598.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=598,cumulative_tests=3040,numbered_scientific_files=1106,
            unique_protected_evidence_files=1775,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=599,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
