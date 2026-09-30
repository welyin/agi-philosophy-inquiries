"""Publish the same-reference limit report with snapshots and resumable writes."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import verify_round525 as evidence

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
checked=evidence.verify()
assert core.read(HERE/'research_round_525_checks.json')==checked
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',
       HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round525_20260930'
resuming=folder.exists()
assert not (HERE/'round525_navigation_checks.json').exists(),'already published'
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
    '**第525轮完成：** [同一参考自触发的宏观极限]({p}research_note_525.md)'
    '在同一有限三场正势量子模型中证明：增大参考幅度及标定，仍不能消除开关对参考的力；'
    '实际探针相对冻结目标的经典输出有严格四阶差。正确连接是完整非线性经典轨道及其Hessian量子涨落，'
    '已用真实Schrödinger传播复算。只排除该冻结近似，未反证认知框架；连续场、空间标定和几何反馈仍开放。'
    f'8组、11式及独立终审通过；最新525／2586，886份编号科学文件，{count}份保护证据；'
    '[核验]({p}research_round_525_checks.json)。本段为最新状态，旧接续保留为历史。')
next_steps='''

### 第525轮后：自洽参考的实际可辨识性

525没有以新增独立控制场替代“被读取的同一参考自触发”问题。有限胞元的正势H_b在每个b下固定且与态无关；跨b也改变F(R/b)标定，这一点保留为参数族输入。经典极限的额外参考力不消失，并传到真实探针；有限b的初始量子二阶偏差与极限四阶差分别记账。正确Hepp连接要求完整非线性轨道及Hessian，不能仅把F(X_cl)塞进原二次矩阵。

**下一轮526候选：** 对同一自洽参考—探针动力学，检验实际记录能否在一个明确参数区域内标定参考，求可辨识范围及量子涨落误差，再考虑多个区域的一致记录。先回查已有可观测性／参考反馈结论；参数可辨不等于物理空间已生成，三个内部场幅不等于三维空间。独立重控制场只作为另列假设的备选，它与被测R的共同标定不能省略。

第525轮正式报告、解析差、真实量子数值与正确宏观近似已保存。下一步不重算同一幅度曲线凑轮次；完整新结论继续及时形成research_note_编号.md。双向逻辑统一目标保持，空间阶段未结项；不设置定时任务、不做图像检验。
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
    assert '**第525轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 170.' not in body
        body+='\n\n## 170. 同源开关的反作用及正确宏观连接\n\n'+block+next_steps
    else:
        assert '\n## 75.' not in body
        body+='\n\n## 75. 第525轮：内部场幅与可用空间仍须区分\n\n'+block
    if path==HERE/'README.md':
        body+=('\n\n## 第525轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|525|[同源参考的自洽极限](research_note_525.md)|[代码](self_gated_reference_limit.py)、'
            '[结果](self_gated_reference_limit_results.json)、[核验](research_round_525_checks.json)|\n')
    if path.name=='RESEARCH_STATE.md':body+=next_steps
    if path.name=='research_direction.md':
        assert '最新科学轮次与检查数为524／2578' in body
        body=body.replace('最新科学轮次与检查数为524／2578','最新科学轮次与检查数为525／2586')
        body=body.replace('完成231—524轮。','完成231—525轮。')
        body+='\n\n**525后接续（目标不变）：** 同源自触发须采用自洽轨道和完整涨落矩阵；下一项检验其实际参考记录的可辨识与标定范围，继续区分内部幅度参数和物理空间。详见RESEARCH_STATE.md末尾。\n'
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
    temporary=path.with_name(path.name+'.round525.tmp')
    if temporary.exists():assert temporary.read_bytes()==body
    else:
        with temporary.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path],('concurrent navigation edit',path)
    os.replace(temporary,path)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[group].items():assert core.digest(HERE/name)==digest,name
report=dict(date='2026-09-30',latest_round=525,cumulative_tests=2586,
    numbered_scientific_files=886,unique_protected_evidence_files=count,
    navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=526,
    active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with (HERE/'round525_navigation_checks.json').open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
