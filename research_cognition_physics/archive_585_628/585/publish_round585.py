"""Publish round 585 and the consolidated work order without changing the goal."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round585_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_585_checks.json')
assert checked['round']==585 and checked['all_reported_checks_passed']
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round585_20261001'
resuming=folder.exists()
if resuming:
    manifest=core.read(folder/'manifest.json')
    raws={p:(folder/manifest[str(p.relative_to(ROOT))]['snapshot']).read_bytes() for p in paths}
    for p,raw in raws.items():
        assert hashlib.sha256(raw).hexdigest()==manifest[str(p.relative_to(ROOT))]['sha256']
else:
    raws={p:p.read_bytes() for p in paths};folder.mkdir(exist_ok=False);manifest={}
    for p,raw in raws.items():
        key=('root_' if p.parent==ROOT else 'research_' if p.parent==RESEARCH else 'archive_')+p.name
        with (folder/key).open('xb') as stream:stream.write(raw)
        manifest[str(p.relative_to(ROOT))]=dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    with (folder/'manifest.json').open('x',encoding='utf8') as stream:
        json.dump(manifest,stream,ensure_ascii=False,indent=2)

summary='**第585轮完成：** [同一量子演化与局部几何来源]({p}research_note_585.md)构造原有限图H的正局部lapse扩展：单位lapse的演化、Gauss、记录与空间权响应全部相同，原Gauss相位准备却给不同局部来源。光滑经典采样差为二阶；无一致能源界的量子来源可保留或放大差。继承的经典法向括号排除所列固定校准下的单独平均接法。五组、十四式及主代理核验通过，最新585／2984，1067份编号科学文件，1630份保护证据；[核验]({p}research_round_585_checks.json)、[条件账]({p}unified_physics_condition_ledger_585.md)。未完成量子约束、shift、协变应力或连续匹配；未取得新的独立代理审查。'
order='**执行顺序：** 沿[合并清单]({p}joint_condition_compression_update_584.md)复用已有成果；585把共同来源缺口定位到完整局部时间与切向生成元。下一步从原H及局部正分配取得能源流，再核与几何动量的匹配，不继续扫描平均参数。统一目标不变。'
next_steps='\n\n### 第585轮后：同一Hamiltonian的能源流与切向几何\n\n原单位lapse模型连同空间权响应还没有决定完整H[N,beta;gamma]。585给正局部扩展的原Gauss见证，并将经典采样精度与量子能源预算一同登记；任意局部时间推进的匹配不能由总能源替代。\n\n586先从原局部正形式推导能源流，核原曲目标测地边作用、Gauss和共同定义域，再对照切向物质生成元。一般连续性恒等式直接复用，不将单边诊断或经典括号当全图量子引力。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第585轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+order.format(p=prefix)+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 231.' not in body
        body+='\n\n## 231. 同一量子演化尚未决定局部几何来源\n\n'+block+next_steps
    else:
        assert '\n## 136.' not in body
        body+='\n\n## 136. 局部法向生成元与共同来源仍需匹配\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[586能源流与切向几何候选](archive_231_/round586_drafts/STATUS.md)已登记，尚未计完成轮次；[最新条件账](archive_231_/unified_physics_condition_ledger_585.md)。\n'
    if path.name=='research_direction.md':
        assert '最新科学轮次与检查数为584／2979' in body
        body=body.replace('最新科学轮次与检查数为584／2979','最新科学轮次与检查数为585／2984')
        assert '完成231—584轮。' in body
        body=body.replace('完成231—584轮。','完成231—585轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第585轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|585|[同一量子演化与局部几何来源](research_note_585.md)|[代码](joint_local_lapse_source.py)、'
            '[结果](joint_local_lapse_source_results.json)、[核验](research_round_585_checks.json)|\n')
        body+='\n[条件账](unified_physics_condition_ledger_585.md)已更新。**586候选：** [能源流与切向几何](round586_drafts/STATUS.md)，尚未计完成。\n'
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
support_links=0
support_files=('unified_physics_condition_ledger_585.md','round586_drafts/STATUS.md','joint_condition_compression_update_584.md')
for name in support_files:
    path=HERE/name
    for link in core.link_parser()(path.read_text('utf8')):
        assert (path.parent/link).resolve().exists(),(path,link)
        support_links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body:continue
    temp=path.with_name(path.name+'.round585.tmp')
    if temp.exists():assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=585,cumulative_tests=2984,numbered_scientific_files=1067,
    unique_protected_evidence_files=1630,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=586,
    supporting_document_links=support_links,
    supporting_document_hashes={name:core.digest(HERE/name) for name in support_files},
    consolidation_not_counted_as_scientific_round=True,active_goal_unchanged=True,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
