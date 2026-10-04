"""Publish round 528 with immutable navigation snapshots and guarded writes."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import verify_round528 as evidence

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
checked=evidence.verify()
assert core.read(HERE/'research_round_528_checks.json')==checked
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',
       HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round528_20260930'
resuming=folder.exists()
assert not (HERE/'round528_navigation_checks.json').exists(),'already published'
if resuming:
    saved_manifest=core.read(folder/'manifest.json')
    raws={p:(folder/saved_manifest[str(p.relative_to(ROOT))]['snapshot']).read_bytes() for p in paths}
    for p,raw in raws.items():
        assert hashlib.sha256(raw).hexdigest()==saved_manifest[str(p.relative_to(ROOT))]['sha256']
else:
    raws={p:p.read_bytes() for p in paths}
    folder.mkdir(exist_ok=False)
count=checked['cumulative_unique_protected_evidence_files']
summary=(
    '**第528轮完成：** [经典共同参考的量子反作用与关系门修订]({p}research_note_528.md)'
    '在527同一规则下证明：经典平坦共同方向在任意正规有限能量量子态中仍受严格同号的额外力，不能直接成为精确自由量子因子；'
    '单胞元归一化均值偏差可随b增大减小，不否定宏观近似。明示改用正关系门后，可同时得到精确自由共同部门与保留非Gaussian相对噪声的实际读口。'
    '非均匀空间来源、连续极限和引力反馈仍开放。'
    f'8组、11式及独立终审通过；最新528／2610，896份编号科学文件，{count}份保护证据；'
    '[核验]({p}research_round_528_checks.json)。本段为最新状态，旧接续保留为历史。')
next_steps='''

### 第528轮后：共同参考的量子稳定与实际空间来源

528已严格区分经典平坦分支与量子精确自由部门：原自门F的正力算符在全部正规有限能量态中严格非零，排除直接拼接；单胞元宏观均值误差又有O(b⁻²)界，所以不把它扩写为宏观空间失败。候选关系门G是明示的规则修订，给精确自由共同因子、全时相对能量界及实际χ读取；非Gaussian相对噪声和meter回冲均保留。有限共同位置二阶矩是读取置信的额外来源条件，不能由自由零模的有限能量替代。

**后继候选先去重：** 回查共同参考准备、热源、521来源及522—523的维数合同，再判断关系门模型能否在有来源的共同参考上保持空间差读的规模一致精度。优先接成熟有限图／无限体积统计力学定理，核清同一状态、实际χ读口及非线性相对部门的范围。重算自由共同模、旧热IR积分、一般二次型对角化或单纯扩大网格不新增529；须有新的来源、稳定性或接口排除结果。

空间来源仍需同一配置内可用的非均匀参考及与传播的适配，内部幅度、站点数或图度数不替代维数。经典自由参考—Einstein嵌入可条件性继承，但量子相对应力不能删除，连续量子场及几何反馈未闭合。研究目标不变，空间阶段未结项。

本轮完整research_note_528.md、代码、结果、终审与历史核验已及时保存；去重对照在报告第1节。保留所有旧稿和冻结文件，用户手动继续，不设置定时任务，不做图像检验。
'''
planned,manifest={},{}
for path,raw in raws.items():
    key=('root_' if path.parent==ROOT else 'research_' if path.parent==RESEARCH else 'archive_')+path.name
    if not resuming:
        with (folder/key).open('xb') as stream:stream.write(raw)
    manifest[str(path.relative_to(ROOT))]=dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第528轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 173.' not in body
        body+='\n\n## 173. 经典共同方向的量子反作用与关系门边界\n\n'+block+next_steps
    else:
        assert '\n## 78.' not in body
        body+='\n\n## 78. 第528轮：共同参考精确自由与宏观近似\n\n'+block
    if path==HERE/'README.md':
        body+=('\n\n## 第528轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|528|[经典共同参考的量子反作用与关系门修订](research_note_528.md)|[代码](quantum_common_reference.py)、'
            '[结果](quantum_common_reference_results.json)、[核验](research_round_528_checks.json)|\n')
    if path.name=='RESEARCH_STATE.md':body+=next_steps
    if path.name=='research_direction.md':
        assert '最新科学轮次与检查数为527／2602' in body
        body=body.replace('最新科学轮次与检查数为527／2602','最新科学轮次与检查数为528／2610')
        body=body.replace('完成231—527轮。','完成231—528轮。')
        body+='\n\n**528后接续（目标不变）：** 原自门的量子自由参考障碍及关系门候选修复已形成正式报告。先去重审查来源与热稳定，检验实际空间差读能否获得规模一致的精度；不能重算已证自由共同模或旧维数阈值。详见RESEARCH_STATE.md末尾。\n'
    planned[path]=body.replace('\n',newline).encode(encoding)
if not resuming:
    with (folder/'manifest.json').open('x',encoding='utf8') as stream:
        json.dump(manifest,stream,ensure_ascii=False,indent=2)
else:assert manifest==saved_manifest
links=0
for path,body in planned.items():
    for link in core.link_parser()(body.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent navigation edit',path)
for path,body in planned.items():
    if path.read_bytes()==body:continue
    temporary=path.with_name(path.name+'.round528.tmp')
    if temporary.exists():assert temporary.read_bytes()==body
    else:
        with temporary.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path],('concurrent navigation edit',path)
    os.replace(temporary,path)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[group].items():assert core.digest(HERE/name)==digest,name
report=dict(date='2026-09-30',latest_round=528,cumulative_tests=2610,
    numbered_scientific_files=896,unique_protected_evidence_files=count,
    navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=529,
    active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with (HERE/'round528_navigation_checks.json').open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
