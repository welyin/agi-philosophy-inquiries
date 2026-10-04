"""Publish round 574; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_curved_quantum_source as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round574_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_574_checks.json')
assert checked['round']==574 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round574_20261001'
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

summary='**第574轮完成：** [曲目标中的共同物质量子初值]({p}research_note_574.md)在声明曲目标量子化及固定初片几何下，将572同一物质的精确Gauss采样接到正规有限能量半经典态；新曲图能量连接原连续初片物质能量。普通平直Gaussian直接移植有范数／能量失败，紧支撑不变态可行。六组、十四式及独立终审通过，最新574／2927，1034份编号科学文件，1504份保护证据；[核验]({p}research_round_574_checks.json)、[条件账]({p}unified_physics_condition_ledger_574.md)。量子化和离散化仍输入，非实际制备、量子连续或量子Einstein约束，统一仍开放。'
next_steps='\n\n### 第574轮后：复用共同量子初值，接实际记录与几何反作用\n\n曲目标、精确Gauss和同一572来源已有固定图半经典初值连接。量子化次序、测地边势与固定几何权是明示输入，两个极限不可交换；不把正规态存在当自主制备或有限物理ℏ的连续量子理论。\n\n已经完成的来源构造直接继承。575回查557—563的实际仪器与完整能源账，核同一曲目标上的记录怎样改变源能量、动量流及几何初始约束。一般Gaussian梯度平方注能公式只作成熟工具，不重报新轮；源、指针及控制应力须分开交代，不能凭Gauss合法宣布引力约束合法。\n\n若只能建立给定脉冲或源平均应力的接口，严格限定其范围，不宣称闭合装置、量子引力或Einstein作用生成。停止波包宽度、网格与单项成本扫描；低曲率、手征物质、维数、引力来源和全尺度统一目标保持。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第574轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 220.' not in body
        body+='\n\n## 220. 共同量子初值与几何来源的条件性连接\n\n'+block+next_steps
    else:
        assert '\n## 125.' not in body
        body+='\n\n## 125. 内部曲目标与量子态不选择时空维数\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[575同一记录与几何反作用候选](archive_231_/round575_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为573／2921','最新科学轮次与检查数为574／2927')
        body=body.replace('完成231—573轮。','完成231—574轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第574轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|574|[曲目标中的共同物质量子初值](research_note_574.md)|[代码](joint_curved_quantum_source.py)、'
            '[结果](joint_curved_quantum_source_results.json)、[核验](research_round_574_checks.json)|\n')
        body+='\n**575候选：** [同一记录与几何反作用](round575_drafts/STATUS.md)，尚未计轮次；复用共同量子初值，核实际记录对源、指针及几何约束的共同影响。\n'
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
    temp=path.with_name(path.name+'.round574.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=574,cumulative_tests=2927,numbered_scientific_files=1034,
    unique_protected_evidence_files=1504,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=575,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
