"""Publish round 547 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_neutral_current_flavour_invariant as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round547_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_547_checks.json')
assert checked['round']==547 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round547_20260930'
resuming=folder.exists()
if resuming:
    manifest=core.read(folder/'manifest.json')
    raws={p:(folder/manifest[str(p.relative_to(ROOT))]['snapshot']).read_bytes() for p in paths}
    for p,raw in raws.items():
        assert hashlib.sha256(raw).hexdigest()==manifest[str(p.relative_to(ROOT))]['sha256']
else:
    raws={p:p.read_bytes() for p in paths}; folder.mkdir(exist_ok=False); manifest={}
    for p,raw in raws.items():
        key=('root_' if p.parent==ROOT else 'research_' if p.parent==RESEARCH else 'archive_')+p.name
        with (folder/key).open('xb') as stream:
            stream.write(raw)
        manifest[str(p.relative_to(ROOT))]=dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    with (folder/'manifest.json').open('x',encoding='utf8') as stream:
        json.dump(manifest,stream,ensure_ascii=False,indent=2)
summary='**第547轮完成：** [共同谱的中性弱流条件]({p}research_note_547.md)在重态通道关闭的精确保护分支中，证明总轻中性流仅依赖总活跃混合量，味重分配不能消除偏差；同一谱同时固定重阈值，并给消去树级归一的宽度／不对称量组合。5组及独立终审通过，最新547／2749，953份编号科学文件，1246份保护证据；[核验]({p}research_round_547_checks.json)、[条件更新]({p}unified_physics_condition_ledger_547.md)。不是实验拟合或全族排除。后继按用户要求先合并已证条件，再攻跨部门接口；统一目标仍开放。'
next_steps='\n\n### 第547轮后：先合并已证条件，再攻共同接口\n\n547完成同一物质候选与味不变树级弱流条件的连接；人造MS目标、历史规范起点、阈值及匹配范围保持明示。\n\n按用户最新要求，先提取全局已证连接、条件性蕴涵、共同实现与剩余独立输入，不自动继续精密电弱支线。优先检查能同时约束组合、物质、实际参考或几何反作用的接口；已有严格矛盾及时处理，电弱数据提取仍保留为必要待办。整理不计新轮次；下一正式编号548应包含新构造、证明或反例。\n\n维数、Lorentz、σ来源、共同几何、场论连接、设备与资源以及完整观测相容性仍开放。保持统一目标，不另立任务或自动化。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第547轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 193.' not in body
        body+='\n\n## 193. 共同物质谱与味不变测量条件\n\n'+block+next_steps
    else:
        assert '\n## 98.' not in body
        body+='\n\n## 98. 弱流条件不由味重分配消失\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为546／2744','最新科学轮次与检查数为547／2749')
        body=body.replace('完成231—546轮。','完成231—547轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第547轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|547|[共同谱的中性弱流条件](research_note_547.md)|[代码](joint_neutral_current_flavour_invariant.py)、'
            '[结果](joint_neutral_current_flavour_invariant_results.json)、[核验](research_round_547_checks.json)|\n')
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body: continue
    temp=path.with_name(path.name+'.round547.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=547,cumulative_tests=2749,numbered_scientific_files=953,
    unique_protected_evidence_files=1246,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=548,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
