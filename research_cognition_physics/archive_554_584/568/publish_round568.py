"""Publish round 568; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_gauge_matter_propagation as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round568_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_568_checks.json')
assert checked['round']==568 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round568_20261001'
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

summary='**第568轮完成：** [同一Gauss物质的纵横传播与共同尺度]({p}research_note_568.md)保留原物质角动能及Gauss约束，得纵向传播；补磁面势并联立系数后，规范纵横与双标量共用谱斜率，完整非线性能量也有给定光滑场下的经典连接。八组、十四式及独立终审通过，最新568／2883，1016份编号科学文件，1446份保护证据；[核验]({p}research_round_568_checks.json)、[条件账]({p}unified_physics_condition_ledger_568.md)、[合并现状]({p}joint_condition_compression_update_568.md)。面胞、匹配与尺度输入仍在，未证明量子连续、三维、完整标准模型或统一。'
next_steps='\n\n### 第568轮后：接已有规范归一，保留未破缺约束\n\n同一V、真空和原Gauss物质已连接到共同玻色谱与经典能量尺度；纵向传播无需另设常数，横向需要明示磁面项及b beta/2=k_h/v。能量一致不等于流、量子相或全部Gauss解的收敛；三维、面胞、群及尺度仍有输入。\n\n后继569回查531—547，将当前SU(2)演示耦合与已有共同规范荷、动能归一和运行真正对应。共同Higgs接入超荷后有非平凡稳定子，须保留未破缺Gauss方向，不能直接套用SU(2)的全部角消去或把三个质量模当光子。\n\n先处理实际同对象约束和系数映射；教科书对角化、568斜率或网格重复计算不新增轮次。同物质装置、量子连续连接、手征费米子和几何反馈仍开放；冻结历史及统一目标保持。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第568轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 214.' not in body
        body+='\n\n## 214. 规范与物质的共同经典传播\n\n'+block+next_steps
    else:
        assert '\n## 119.' not in body
        body+='\n\n## 119. 共同玻色传播仍不选择空间维数\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[569共同超荷与规范归一候选](archive_231_/round569_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为567／2875','最新科学轮次与检查数为568／2883')
        body=body.replace('完成231—567轮。','完成231—568轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第568轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|568|[Gauss物质的纵横传播与尺度](research_note_568.md)|[代码](joint_gauge_matter_propagation.py)、'
            '[结果](joint_gauge_matter_propagation_results.json)、[核验](research_round_568_checks.json)|\n')
        body+='\n**569候选：** [共同超荷与规范归一](round569_drafts/STATUS.md)，尚未计算；回查已有群荷与运行，核Higgs稳定子及未破缺Gauss方向，不重复已完成谱和网格检查。\n'
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
    temp=path.with_name(path.name+'.round568.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=568,cumulative_tests=2883,numbered_scientific_files=1016,
    unique_protected_evidence_files=1446,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=569,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))

