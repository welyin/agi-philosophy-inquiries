"""Publish the finite self-gated decoder with snapshots and resumable writes."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import verify_round526 as evidence

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
checked=evidence.verify()
assert core.read(HERE/'research_round_526_checks.json')==checked
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',
       HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round526_20260930'
resuming=folder.exists()
assert not (HERE/'round526_navigation_checks.json').exists(),'already published'
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
    '**第526轮完成：** [自触发参考的有限标定与分支边界]({p}research_note_526.md)'
    '在525同一规则下，对未知参考幅度的明确区间证明统一可读：严格有限标定表、1201符号记录与有限资源下，'
    '单次恢复误差≤1/16的成功率至少99.5198%。完整量子输出误差≤19/b；全实轴存在同均值分支，'
    '其单次末读不能靠增大b取得全局一致恢复。只涉及内部参考幅度，未生成物理位置、三维或GR。'
    f'8组、13式及独立终审通过；最新526／2594，889份编号科学文件，{count}份保护证据；'
    '[核验]({p}research_round_526_checks.json)。本段为最新状态，旧接续保留为历史。')
next_steps='''

### 第526轮后：从本地可读参考到区域间一致记录

526补上525的实际可辨识性：固定同一H_b、同一T和同一有限解码表，对x∈[3/4,5/4]取得统一有限资源置信保证。未知x只改变准备；跨b仍改变门标定。区间算术给完整有限计算证书；b=2²¹的保证来自解析界，量子诊断仅模拟到4096。已有通用逆函数、Gaussian统计和维数阈值没有重报为新结果。全实轴折叠只排除该来源族／单χ读口的全局一致恢复。

**下一轮527候选：** 先对照局域参考、521—524及已有参考几何接口，选择同一参考结构，检验不同区域实际联合记录如何在同一交互下组成一致坐标图。必须保留来源关联、区域间耦合与反作用；不能把独立胞元复制后直接称为空间，也不重算521的线性读取。预给梯度或几何若仍需要，就明确为双向统一的接口输入并检查相容性。优先判断本地标定与共同参考是否能真正接起来，必要时给严格接口障碍；不继续优化同一个b或新增无来源维数。

本轮正式报告及时保存为research_note_526.md，代码、结果、区间证书、独立审查和全历史核验同目录归档。研究目标保持，不设置定时任务；无图像检验。三维及空间阶段尚未结项。
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
    assert '**第526轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 171.' not in body
        body+='\n\n## 171. 同一自触发参考的有限实际解码\n\n'+block+next_steps
    else:
        assert '\n## 76.' not in body
        body+='\n\n## 76. 第526轮：可读参考幅度仍不等于空间位置\n\n'+block
    if path==HERE/'README.md':
        body+=('\n\n## 第526轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|526|[自触发参考的有限解码](research_note_526.md)|[代码](self_gated_reference_calibration.py)、'
            '[结果](self_gated_reference_calibration_results.json)、[核验](research_round_526_checks.json)|\n')
    if path.name=='RESEARCH_STATE.md':body+=next_steps
    if path.name=='research_direction.md':
        assert '最新科学轮次与检查数为525／2586' in body
        body=body.replace('最新科学轮次与检查数为525／2586','最新科学轮次与检查数为526／2594')
        body=body.replace('完成231—525轮。','完成231—526轮。')
        body+='\n\n**526后接续（目标不变）：** 同源参考在明确区间已有实际有限解码；下一项检验多个区域的联合记录与共同标定能否在同一规则下相容，引用旧定位工具，区分内部幅度和空间位置。详见RESEARCH_STATE.md末尾。\n'
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
    temporary=path.with_name(path.name+'.round526.tmp')
    if temporary.exists():assert temporary.read_bytes()==body
    else:
        with temporary.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path],('concurrent navigation edit',path)
    os.replace(temporary,path)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[group].items():assert core.digest(HERE/name)==digest,name
report=dict(date='2026-09-30',latest_round=526,cumulative_tests=2594,
    numbered_scientific_files=889,unique_protected_evidence_files=count,
    navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=527,
    active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with (HERE/'round526_navigation_checks.json').open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
