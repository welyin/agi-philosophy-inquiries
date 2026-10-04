"""Publish round 537 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import seesaw_boundary_compatibility as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round537_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_537_checks.json')
assert checked['round']==537 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round537_20260930'
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
summary='**第537轮完成：** [轻中微子与共同Higgs边界]({p}research_note_537.md)给出三ν固定Dirac奇异值和Majorana上界下的最小树级Weinberg范数，并联立同一规范／Higgs矩。人工同尺度预算排除536的λ=0.2见证，但保留rank-one复相关边界；零预算时top＋三ν族锐界为λ≥g_w²。8组与独立终审通过，最新537／2679，923份编号科学文件，1172份保护证据；[核验]({p}research_round_537_checks.json)、[条件更新]({p}unified_physics_condition_ledger_537.md)。尚未做低能到高能运行、阈值或振荡数据拟合，整体目标开放。'
next_steps='\n\n### 第537轮后：存活复结构的对称性、运行与阈值\n\n537将536两个物质矩接到tree-level C₅=Y M⁻¹Yᵀ。三个右手中微子、正重质量≤M_max、给定z₁≥z₂≥z₃下，遍历所有具有该奇异值的复Y及允许M，精确min||C₅||=max(z₂,√(z₁z₃))/M_max；不仅是固定矩阵模长调相位。η是同尺度系数预算，非未经运行便直接搬来的实验界。\n\n联立top＋三ν的共同边界，有限预算给必要λ≥g_w²−4M_maxη，零预算锐界λ≥g_w²。536的λ=0.2见证不满足所给人工预算，但同一规范参数下λ=g_w²的rank-one复相关候选存活；零轻质量不是观测拟合。数值领先换算约355.63eV只用于旧见证的诊断，不是已经计算物理低能谱。\n\n下一正式编号538选择存活的近似轻子数分支，核同一作用与简并重质量下的RG、阈值和小破缺。对接Kersten–Smirnov 0705.3221及成熟矩阵RG，不重新证明已知rank-one分类；不能把一般精细抵消当自动辐射稳定。随后将人工η换为来源明确的低能到匹配尺度条件，同时检查Higgs和top，不新增无来源无菌场或另调规范权重。\n\n保留同尺度、规范化及三活跃ν_R条件；质量项、共同几何、实际记录、Lorentz和量子引力仍开放。继续及时编号报告和总账，不改应用目标、不设定时任务、不做图像检验。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第537轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 183.' not in body
        body+='\n\n## 183. 轻中微子条件对共同参数的收紧\n\n'+block+next_steps
    else:
        assert '\n## 88.' not in body
        body+='\n\n## 88. 联合物质条件不等于几何生成\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为536／2671','最新科学轮次与检查数为537／2679')
        body=body.replace('完成231—536轮。','完成231—537轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第537轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|537|[轻中微子与Higgs共同边界](research_note_537.md)|[代码](seesaw_boundary_compatibility.py)、'
            '[结果](seesaw_boundary_compatibility_results.json)、[核验](research_round_537_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round537.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=537,cumulative_tests=2679,numbered_scientific_files=923,
    unique_protected_evidence_files=1172,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=538,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
