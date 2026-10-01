"""Publish round 570; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_gauss_nonlinear_compatibility as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round570_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_570_checks.json')
assert checked['round']==570 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round570_20261001'
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

summary='**第570轮完成：** [线性传播的非线性总荷与内部补全]({p}research_note_570.md)证明当前闭图真空的一阶Gauss扰动还需二次总荷条件；单环切向有不可积反例，两环可用原Higgs角荷及内部圆通量精确补全，新增初值动能O(eta⁴)。小残差也不保证近稳定子时小动量修正。八组、十四式及独立终审通过，最新570／2899，1022份编号科学文件，1467份保护证据；[核验]({p}research_round_570_checks.json)、[条件账]({p}unified_physics_condition_ledger_570.md)。仅限声明的经典初值范围，未反驳特定连续采样或完成引力、量子连续及统一。'
next_steps='\n\n### 第570轮后：从兼容来源建立受控映射\n\n原Gauss已经约束线性扰动的二次总荷，不需要新增认知公理；反例只是当前闭图指定切向不可积。两环用同一物质和内部通量可精确补全，证明不是整个模型失败。近稳定子的小残差不能代替误差证明，固定图预算也不能外推规模一致性。\n\n571优先对接成熟约束兼容离散化，给精确连续Gauss源到离散物理初值的明确映射，并保留同一物质应力、边界、尺度与范数。不能只复述伪逆或把任意小残差当成功；先处理可继承的通量结构，再核真正非Abelian/Higgs交叉缺口。\n\n549的无规范场参考与Einstein约束仍直接继承。图源、维数、连续量子、手征物质、真实装置和引力动力学来源继续开放；目标保持。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第570轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 216.' not in body
        body+='\n\n## 216. 线性扰动的非线性总荷与内部补全\n\n'+block+next_steps
    else:
        assert '\n## 121.' not in body
        body+='\n\n## 121. 非线性Gauss兼容不选择空间维数\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[571连续源与离散约束候选](archive_231_/round571_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为569／2891','最新科学轮次与检查数为570／2899')
        body=body.replace('完成231—569轮。','完成231—570轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第570轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|570|[线性到非线性Gauss兼容](research_note_570.md)|[代码](joint_gauss_nonlinear_compatibility.py)、'
            '[结果](joint_gauss_nonlinear_compatibility_results.json)、[核验](research_round_570_checks.json)|\n')
        body+='\n**571候选：** [连续源与离散约束](round571_drafts/STATUS.md)，尚未计轮次；先核精确连续Gauss源的约束兼容映射，保留背景稳定子、边界及误差范围。\n'
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
    temp=path.with_name(path.name+'.round570.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=570,cumulative_tests=2899,numbered_scientific_files=1022,
    unique_protected_evidence_files=1467,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=571,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
