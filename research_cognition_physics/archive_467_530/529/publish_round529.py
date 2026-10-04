"""Publish round 529, retaining snapshots and checking concurrent edits."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import verify_round529 as evidence

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
checked=evidence.verify()
assert core.read(HERE/'research_round_529_checks.json')==checked
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',
       HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round529_20260930'
resuming=folder.exists()
assert not (HERE/'round529_navigation_checks.json').exists(),'already published'
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
    '**第529轮完成：** [非线性关系门的量子热稳定性]({p}research_note_529.md)'
    '证明528同一关系门的全局凸性，并以精确量子路径表示和成熟BL工具，得到真实非Gaussian量子热态的全有限图统一点读／差读界及指数尾。'
    '实际χ仪器因此可继承522的热参考判据；有限能量密度的奇偶混合反例显示热来源条件不能省去。'
    '热化、仿射参考来源、几何与三维唯一性仍开放。'
    f'6组、11式及独立终审通过；最新529／2616，899份编号科学文件，{count}份保护证据；'
    '[核验]({p}research_round_529_checks.json)。本段为最新状态，旧接续保留为历史。')
next_steps='''

### 第529轮后：转向非均匀共同参考的来源

529关闭528明确遗留的相对热噪声规模问题：同一非线性G具有Hessian下界K_c，最低本征值1/4；用精确谐环路的BL估计及受支配虚时间极限，得到真正量子Gibbs态的等时协方差和subGaussian尾。实际χ差读保留非Gaussian相对项与meter噪声，因而522的d>2阈值可直接继承，不重复登记旧IR证明。量子热态是新增来源条件；371N/80能量上界的非热奇偶来源可使单站方差随N增长，故有限能量密度不替代该条件。

**下一项先回查来源：** 检索旧历史环路、仿射参考、twist／支路、参考对齐及287图册工具，判断非均匀且空间分离的共同参考须额外指定哪些边界或全局扇区，是否可由已有关系记录形成。必须在同一规则下给真正来源构造或准确来源障碍；已有环路和梯度恒等式直接引用，不能换成新轮。也不再重复谐模自由演化、BL、热IR阈值或更多噪声常数优化。只有实质来源／接口增量才正式开启530。

保留完整无质量零模的正规准备：它在同时差读中消去，并不具有归一化全Gibbs态。热化、仿射twist、读者控制及全局记录汇集仍有输入。非线性相对部门的有限图量子结论未自动给连续场、各向同性球面方向壳或量子应力—几何闭合；三维唯一性与空间阶段仍未完成。

完整research_note_529.md、代码、结果、独立终审及全历史核验已及时归档，原稿和冻结证据保留。目标不变，用户手动继续，不设置定时任务，不做图像检验。
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
    assert '**第529轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 174.' not in body
        body+='\n\n## 174. 非线性共同参考的量子热来源与统一噪声\n\n'+block+next_steps
    else:
        assert '\n## 79.' not in body
        body+='\n\n## 79. 第529轮：真实非Gaussian热噪声与旧维数阈值\n\n'+block
    if path==HERE/'README.md':
        body+=('\n\n## 第529轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|529|[非线性关系门的量子热稳定性](research_note_529.md)|[代码](nonlinear_thermal_reference.py)、'
            '[结果](nonlinear_thermal_reference_results.json)、[核验](research_round_529_checks.json)|\n')
    if path.name=='RESEARCH_STATE.md':body+=next_steps
    if path.name=='research_direction.md':
        assert '最新科学轮次与检查数为528／2610' in body
        body=body.replace('最新科学轮次与检查数为528／2610','最新科学轮次与检查数为529／2616')
        body=body.replace('完成231—528轮。','完成231—529轮。')
        body+='\n\n**529后接续（目标不变）：** 非线性相对部门的真实量子热态已有全规模差读与指数尾界；转向非均匀共同参考来源，先回查环路、twist和图册旧工具，避免重做噪声优化或IR阈值。详见RESEARCH_STATE.md末尾。\n'
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
    temporary=path.with_name(path.name+'.round529.tmp')
    if temporary.exists():assert temporary.read_bytes()==body
    else:
        with temporary.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path],('concurrent navigation edit',path)
    os.replace(temporary,path)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[group].items():assert core.digest(HERE/name)==digest,name
report=dict(date='2026-09-30',latest_round=529,cumulative_tests=2616,
    numbered_scientific_files=899,unique_protected_evidence_files=count,
    navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=530,
    active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with (HERE/'round529_navigation_checks.json').open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
