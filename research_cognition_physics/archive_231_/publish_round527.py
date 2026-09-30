"""Publish the finite self-gated decoder with snapshots and resumable writes."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import verify_round527 as evidence

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
checked=evidence.verify()
assert core.read(HERE/'research_round_527_checks.json')==checked
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',
       HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round527_20260930'
resuming=folder.exists()
assert not (HERE/'round527_navigation_checks.json').exists(),'already published'
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
    '**第527轮完成：** [跨区域耦合与共同标定]({p}research_note_527.md)'
    '保留自触发参考并加入持续的三场图耦合，对任意有限有理图及加权度≤3证明联合逆增益、有限三步解码和相关量子噪声界。'
    '每站5601符号，在明确资源族下整组恢复成功率至少99.7868%；孤立表存在严格非零串扰偏差。'
    '同一模型也给配置映射满秩而各处参考值相同的对照：配置可读不推出空间位置，三维和GR仍开放。'
    f'8组、13式及独立终审通过；最新527／2602，893份编号科学文件，{count}份保护证据；'
    '[核验]({p}research_round_527_checks.json)。本段为最新状态，旧接续保留为历史。')
next_steps='''

### 第527轮后：从联合配置标定到空间分离的共同参考来源

527在同一真实跨区域耦合中，接上526的本地可读参考。任意有限有理图、加权度≤3、未知幅度立方体上的联合增益及有限解码已证；区域相关量子噪声完整保留。数值只传播经典中心和Hessian Gaussian，完整量子保证来自Duhamel解析界。自适应区间算法逐有限图终止，不把几个小图当全规模固定精度证明。常量来源说明配置逆不推出同一配置内位置分离；旧坐标定理及线性读取实验没有重复登记。

**下一轮528候选：** 检验空间分离的共同参考来源能否在本轮同一非线性正势规则中成立，并将它与已有关系几何、实际读取接成同一模型。先回查共同参考正势、521的来源差异和物理参考几何；可考察原M=0、cR−χ=0的共同方向，但不能仅作一次线性对角化就算新结果，须处理真实来源、涨落或非线性反作用中的新缺口。邻接与读数差的相容性、图册域和维数条件继续显式标注，不把增加参考分量当维数推导。

本轮正式报告及时保存为research_note_527.md，代码、结果、区间证书、独立审查和全历史核验同目录归档。研究目标保持，不设置定时任务；无图像检验。三维及空间阶段尚未结项。
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
    assert '**第527轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 172.' not in body
        body+='\n\n## 172. 跨区域联合配置可读与空间分离来源\n\n'+block+next_steps
    else:
        assert '\n## 77.' not in body
        body+='\n\n## 77. 第527轮：联合配置逆与同配置位置分离\n\n'+block
    if path==HERE/'README.md':
        body+=('\n\n## 第527轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|527|[跨区域耦合与共同标定](research_note_527.md)|[代码](coupled_reference_calibration.py)、'
            '[结果](coupled_reference_calibration_results.json)、[核验](research_round_527_checks.json)|\n')
    if path.name=='RESEARCH_STATE.md':body+=next_steps
    if path.name=='research_direction.md':
        assert '最新科学轮次与检查数为526／2594' in body
        body=body.replace('最新科学轮次与检查数为526／2594','最新科学轮次与检查数为527／2602')
        body=body.replace('完成231—526轮。','完成231—527轮。')
        body+='\n\n**527后接续（目标不变）：** 区域耦合中的联合有限标定已建立；下一项检验同一规则内空间分离的共同参考来源、实际读取与关系几何的相容性，保留配置逆和位置分离的量词区别。详见RESEARCH_STATE.md末尾。\n'
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
    temporary=path.with_name(path.name+'.round527.tmp')
    if temporary.exists():assert temporary.read_bytes()==body
    else:
        with temporary.open('xb') as stream:stream.write(body)
    assert path.read_bytes()==raws[path],('concurrent navigation edit',path)
    os.replace(temporary,path)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[group].items():assert core.digest(HERE/name)==digest,name
report=dict(date='2026-09-30',latest_round=527,cumulative_tests=2602,
    numbered_scientific_files=893,unique_protected_evidence_files=count,
    navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=528,
    active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with (HERE/'round527_navigation_checks.json').open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
