"""Publish verified 694, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round694_navigation_checks.json'
checked=core.read(HERE/'research_round_694_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round694_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第694轮完成：** [原物理来源临界谱跳变的严格证书]({p}research_note_694.md)'
         '整数惯性与原512维N的逆残差证明固定合法E的标量来源存在非零单侧跳变，'
         '排除来源自动平滑消去近零谱的普遍接法。两组、十四式通过，最新694／3295，'
         '1394份编号科学文件、2681份保护证据。'
         '[核验]({p}research_round_694_checks.json)、[全条件账]({p}unified_physics_condition_ledger_694.md)。'
         '尚非完整辅助／规范平均结论，不构成RP、H_F或连续时空证明。')
order=('**当前执行顺序（694后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[695完整辅助平均与临界支撑]({p}round695_drafts/STATUS.md)，'
       '核原自旋配对及相位，不能从固定E非零直接推积分非零；旧空间接口逐项复用，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第694轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 340.' not in text
        text+='\n\n## 340. 原物理来源临界支撑与精确证书\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 245.' not in text
        text+='\n\n## 245. 旧空间合同复用，规范路径不是空间坐标\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—693轮。' in text and '最新科学轮次与检查数为693／3293' in text
        text=text.replace('完成231—693轮。','完成231—694轮。').replace('最新科学轮次与检查数为693／3293','最新科学轮次与检查数为694／3295')
    if p==HERE/'README.md':
        text+='\n\n## 第694轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|694|[原物理来源临界谱跳变的严格证书](research_note_694.md)|[代码](joint_critical_source_certificate.py)、[结果](joint_critical_source_certificate_results.json)、[核验](research_round_694_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round694.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=694,cumulative_tests=3295,numbered_scientific_files=1394,
            unique_protected_evidence_files=2681,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=695,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
