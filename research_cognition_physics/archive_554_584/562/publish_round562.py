"""Publish round 562; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_radial_link_alignment as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round562_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_562_checks.json')
assert checked['round']==562 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round562_20260930'
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

summary='**第562轮完成：** [协变边读数与径向参考的共同能源条件]({p}research_note_562.md)证明原F读口含内部取向贡献；在明确幅度窗口及固定H系数下，全源角向能源限制该贡献消失，未经修正的径向解释在细化时要求来源预算至少Ω(ε⁻²)。给正规Gauss对齐族，保留两端物质与电能；另读半径、扣角项及同步尺度变化仍允许。6组、12式及独立终审通过，最新562／2847，998份编号科学文件，1388份保护证据；[核验]({p}research_round_562_checks.json)、[条件账]({p}unified_physics_condition_ledger_562.md)。ε仅为接口标定，未建立连续极限、三维或统一。'
next_steps='\n\n### 第562轮后：检验同一全源的直接径向读口\n\n已经将未经修正的协变边读数、内部取向准备和完整源能源合并为必要条件；固定预算下不能把角项自动当零。该限制不覆盖直接径向读取。\n\n后继563复用556的径向域边界、557指针及561有限时工具，检查R=(r1-r2)²/2的真实不变仪器与同一来源。用Cartesian弱梯度处理轴点，核定非选择瞬时与有限时电能变化的区别；明确新增RP耦合、准备、开关和读者输入。不继续优化角向对齐常数。\n\n统一目标不改；径向记录不等于完整连续参考或几何反馈。冻结历史、已证反例及替代方案全部保留。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第562轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 208.' not in body
        body+='\n\n## 208. 协变读数、径向解释与来源能源\n\n'+block+next_steps
    else:
        assert '\n## 113.' not in body
        body+='\n\n## 113. 内部取向不等于径向坐标\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[563同一源直接径向读取候选](archive_231_/round563_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为561／2841','最新科学轮次与检查数为562／2847')
        body=body.replace('完成231—561轮。','完成231—562轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第562轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|562|[协变边读数与径向参考](research_note_562.md)|[代码](joint_radial_link_alignment.py)、'
            '[结果](joint_radial_link_alignment_results.json)、[核验](research_round_562_checks.json)|\n')
        body+='\n**563候选：** [同一完整源的直接径向读取](round563_drafts/STATUS.md)，尚未完成；复用556／557／561，先核域、同一来源、实际记录及反作用，不将原F的角向成本限制扩大为全部径向测量。\n'
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
    temp=path.with_name(path.name+'.round562.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=562,cumulative_tests=2847,numbered_scientific_files=998,
    unique_protected_evidence_files=1388,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=563,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))

