"""Publish the two completed reference reports without changing the active goal."""
from pathlib import Path
import hashlib
import json
import os
import verify_reference_acquisition_precision_interfaces as evidence
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
before = evidence.verify()
assert before['total_protected_including_this_review'] == 1047
paths = [ROOT/'README.md', RESEARCH/'README.md', RESEARCH/'research_direction.md',
         RESEARCH/'RESEARCH_STATE.md', HERE/'README.md',
         HERE/'spatial_premise_closure_audit.md', HERE/'three_dimensional_four_conditions_review.md']
raws = {p:p.read_bytes() for p in paths}
folder = HERE/'navigation_before_reference_acquisition_precision_20260930'
folder.mkdir(exist_ok=False)
summary = (
    '**正势接入与宏观精度报告已完成（未编号）：** '
    '[空白探针接入]({p}positive_reference_acquisition_review.md)在同一正势规则下修正双方源，'
    '证明任意有限盒中的实际坐标读取与有限记录；局域规则固定，读取资源随规模增长。'
    '相对模式仍振荡，未自动形成长期对齐。'
    '[热参考精度]({p}thermal_reference_dimension_review.md)把量子／热涨落接到真实双端仪器：'
    '固定正温、梯度和读取噪声下，逐点对一致的绝对方差上界当且仅当d>2；真空为d>1。'
    '局部或相对精度允许低维，四维以上也通过，尚未唯一选择三维。'
    '两项共12组诊断、21式及独立复核通过；'
    '[核验]({p}reference_acquisition_precision_checks.json)。编号仍520／2552，871份编号科学文件；保护证据1047份。')
next_steps = '''

### 当前接续：把维数上下界放到同一实际测量框架

两篇正式报告及代码、结果均已保存。空白探针报告关闭了正势稳定性与实际坐标读取的具体连接；热参考报告关闭了声明平衡来源下的实际读数噪声与下临界维数连接。两项使用不同准备与边界，不能把前者的保守演化当作后者的热化来源。

**下一实质检验：** 复用382已经证明的固定迹连续qubit反向端口上界d≤3，检验它能否与热参考的d>2条件在同一有效场、同一坐标差读和实际qubit仪器内同时实现。须交代连续方向边界、坐标与方向的对应、来源、控制及参考；不能把有限格点的少量方向当完整S^(d−1)，也不能只把两个独立定理并列。

热下界要求固定绝对精度，单纯方向识别未必需要它。若方向信号随距离增长可压过低维噪声，则应给出该区别的反例，而不能将“能够认同方向”偷换成绝对精度。允许把严格较强的要求作为候选原则研究，但要单列其动机和输入身份。采用更大方向载体也可能取消三维上界，须保留此边界。

本项不改变用户的研究目标。继续按双向统一路线推进，有新结果及时形成正式报告并更新索引；不积压为未编号草稿，不重复编号成熟定理。接合若只得到条件性三维，应明确哪些条件共同可实现、哪些尚未从底层生成。原动态图到连续场、参考准备与维护、共同时间及自洽几何反馈仍待接合。

最新冻结入口verify_reference_acquisition_precision_interfaces.py，1047份保护证据；编号520／2552保持。空间阶段及完整统一未结项。用户手动继续，不设置定时任务。
'''
manifest, planned = {}, {}
for p, raw in raws.items():
    key = ('root_' if p.parent==ROOT else 'research_' if p.parent==RESEARCH else 'archive_')+p.name
    with (folder/key).open('xb') as stream: stream.write(raw)
    manifest[str(p.relative_to(ROOT))] = dict(snapshot=key, sha256=hashlib.sha256(raw).hexdigest())
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl = '\r\n' if b'\r\n' in raw else '\n'
    body = raw.decode(enc).replace('\r\n','\n')
    assert '正势接入与宏观精度报告已完成（未编号）' not in body
    prefix = 'research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    block = summary.format(p=prefix)
    if p in paths[:5]:
        title, tail = body.split('\n\n',1)
        body = title+'\n\n'+block+'\n\n'+tail
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 165.' not in body and '\n## 166.' not in body
        body += '\n\n## 165. 正势下独立探针的实际坐标记录\n\n'+block
        body += '\n\n## 166. 热参考下界与实际方向上界的接合任务\n'+next_steps
    else:
        assert '\n## 70.' not in body and '\n## 71.' not in body
        body += '\n\n## 70. 读取稳定与相对对齐需要分开\n\n'+block
        body += '\n\n## 71. 已有上界与新精度下界不能直接视为共同模型\n\n382的上界可直接复用；热参考下界依赖另一组来源与精度条件。下一项必须给实际测量与连续方向的共同实现及条件账，不能自动宣布三维完成。\n'
    if p.name=='RESEARCH_STATE.md': body += next_steps
    if p.name=='research_direction.md':
        body += '\n\n**最新接续（目标不变）：** 正势空白探针与热参考精度已成正式报告。下一项将382方向上界与热参考下界接到共同测量框架；重点核验同一维数、实际仪器及绝对／方向精度区别。详见RESEARCH_STATE.md末尾。\n'
    planned[p] = body.replace('\n',nl).encode(enc)
with (folder/'manifest.json').open('x',encoding='utf8') as stream:
    json.dump(manifest,stream,ensure_ascii=False,indent=2)
count = 0
for p, body in planned.items():
    for link in core.link_parser()(body.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(), (p,link)
        count += 1
for p, raw in raws.items(): assert p.read_bytes()==raw, ('concurrent change',p)
for p, body in planned.items():
    temp = p.with_name(p.name+'.reference-acquisition.tmp')
    with temp.open('xb') as stream: stream.write(body)
    os.replace(temp,p)
assert evidence.verify()==before
report = dict(date='2026-09-30', scientific_base_through_round=520,
    unchanged_numbered_scientific_tests=2552, numbered_scientific_files=871,
    protected_evidence=1047, navigation_files=7, navigation_links=count, broken_links=0,
    previous_navigation_snapshots_preserved=True, independent_cross_review_completed=True,
    scientific_evidence_unchanged_after_navigation=True, active_goal_unchanged=True,
    manual_goal_start_preserved=True, stage_complete=False, all_checks_passed=True)
with (HERE/'reference_acquisition_precision_integration_checks.json').open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
