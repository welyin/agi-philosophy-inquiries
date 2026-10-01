"""Publish round 571; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_gauss_continuum_sampling as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round571_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_571_checks.json')
assert checked['round']==571 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round571_20261001'
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
        with (folder/key).open('xb') as stream: stream.write(raw)
        manifest[str(p.relative_to(ROOT))]=dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    with (folder/'manifest.json').open('x',encoding='utf8') as stream:
        json.dump(manifest,stream,ensure_ascii=False,indent=2)

summary='**第571轮完成：** [连续共同源与精确Gauss初值]({p}research_note_571.md)给固定光滑周期源到同一图上精确Gauss数据的映射：半边输运、全局总荷修正及原Higgs／圆通量共同补全，动能范数误差O(epsilon²)，保留原径向、横向和谐和来源。八组、十四式及独立终审通过，最新571／2907，1025份编号科学文件，1475份保护证据；[核验]({p}research_round_571_checks.json)、[条件账]({p}unified_physics_condition_ledger_571.md)。限声明的非零Higgs和背景分支；几何、时间演化、量子连续及统一仍开放。'
next_steps='\n\n### 第571轮后：接同一来源的演化与几何\n\n固定连续Gauss源与同一离散物质的精确初值已建立受控连接。原径向及圆无散来源须保留，不能从零重建后只凭Gauss成立签收。误差常数依赖固定光滑背景、Higgs下界与分支；不把特定源的可行映射升级为所有背景统一保证。\n\n572优先核同一来源的实际演化与几何接口：对接成熟约束保持定理，检查完整磁势、Higgs和尺度是否满足其前提；将Gauss密度动量与真实应力及Einstein初始约束分别交代。549无规范场构造及550—553真空／曲率限制直接继承，不能暗换绝对势常数或裸／有效参数。\n\n不继续优化本轮Poisson常数或增加相同网格；只为真实新映射、冲突或构造编号。维数、手征物质、量子连续、实际装置及引力来源开放，目标保持。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第571轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 217.' not in body
        body+='\n\n## 217. 共同连续源的精确Gauss取样\n\n'+block+next_steps
    else:
        assert '\n## 122.' not in body
        body+='\n\n## 122. 精确初值映射不选择空间维数\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[572共同来源与演化几何候选](archive_231_/round572_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为570／2899','最新科学轮次与检查数为571／2907')
        body=body.replace('完成231—570轮。','完成231—571轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第571轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|571|[连续共同源与精确Gauss初值](research_note_571.md)|[代码](joint_gauss_continuum_sampling.py)、'
            '[结果](joint_gauss_continuum_sampling_results.json)、[核验](research_round_571_checks.json)|\n')
        body+='\n**572候选：** [共同来源与演化几何](round572_drafts/STATUS.md)，尚未计轮次；先核成熟演化接口及同物质应力，保留绝对真空能和各类约束。\n'
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
    temp=path.with_name(path.name+'.round571.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=571,cumulative_tests=2907,numbered_scientific_files=1025,
    unique_protected_evidence_files=1475,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=572,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
