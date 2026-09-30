"""Publish reviewed round 524; preserve navigation snapshots and detect conflicts."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import verify_round524 as evidence

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
checked = evidence.verify()
assert core.read(HERE/'research_round_524_checks.json') == checked
paths = [ROOT/'README.md', RESEARCH/'README.md', RESEARCH/'research_direction.md',
         RESEARCH/'RESEARCH_STATE.md', HERE/'README.md',
         HERE/'spatial_premise_closure_audit.md', HERE/'three_dimensional_four_conditions_review.md']
folder = HERE/'navigation_before_round524_20260930'
resuming = folder.exists()
assert not (HERE/'round524_navigation_checks.json').exists(), 'already published'
if resuming:
    saved_manifest = core.read(folder/'manifest.json')
    raws = {p: (folder/saved_manifest[str(p.relative_to(ROOT))]['snapshot']).read_bytes()
            for p in paths}
    for p, raw in raws.items():
        assert hashlib.sha256(raw).hexdigest() == saved_manifest[str(p.relative_to(ROOT))]['sha256']
else:
    raws = {p: p.read_bytes() for p in paths}
    folder.mkdir(exist_ok=False)
count = checked['cumulative_unique_protected_evidence_files']
summary = (
    '**第524轮完成：** [紧支撑正势探针与实际局域读取]({p}research_note_524.md)'
    '在已给3+1背景上，以有限时空支撑的三场正势替代523的硬频率截断；'
    '有限非零耦合有严格正参考增益，同一散射给实际噪声、因果组合与有限末读回冲。'
    '热阈值和qubit上界直接引用旧轮次；终端读取、控制、来源及背景仍输入，未生成时空或闭合GR反馈。'
    f'8组复算、14式及独立终审通过；最新524／2578，883份编号科学文件，{count}份保护证据；'
    '[核验]({p}research_round_524_checks.json)。'
    '本段为最新状态，后文旧轮次与接续说明作为历史保留。')
next_steps = '''

### 第524轮后：动态参考开关的共同模型接口

完整报告已按research_note_524.md归档；8组新诊断没有重复登记热IR阈值、Borsuk–Ulam上界、一般Noether功账或旧Gaussian统计。历史索引、相关原报告及既有接口的去重对照见524第1节；图像检验未执行。

**下一轮525的候选问题：** 在同一动态参考模型中，以局域参考读数开关ρ=F(X)替代本轮外给ρ(t,x)，是否存在明确受控极限，能同时给实际读数、参考反作用及误差界？先对照参考钟仪器、物理参考几何接口及325，复用已有能量恒等式。只写另一个无关模型、换名或重新计算维数阈值不算进展；实际新推导及结果应及时进入编号报告。

本轮仍给定连续3+1背景，不声称原动态图已产生三维、Lorentz结构或Einstein方程。允许成熟物理输入的双向逻辑统一目标保持不变；末端读者暂作明示操作输入，不把重复向后追问读者作为唯一前进方式。空间阶段未结项，不改当前应用目标，不设置定时任务。
'''
planned, manifest = {}, {}
for path, raw in raws.items():
    key = ('root_' if path.parent == ROOT else 'research_' if path.parent == RESEARCH else 'archive_')+path.name
    if not resuming:
        with (folder/key).open('xb') as stream:
            stream.write(raw)
    manifest[str(path.relative_to(ROOT))] = dict(snapshot=key, sha256=hashlib.sha256(raw).hexdigest())
    encoding = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline = '\r\n' if b'\r\n' in raw else '\n'
    body = raw.decode(encoding).replace('\r\n', '\n')
    assert '**第524轮完成：**' not in body
    prefix = 'research_cognition_physics/archive_231_/' if path.parent == ROOT else 'archive_231_/' if path.parent == RESEARCH else ''
    block = summary.format(p=prefix)
    if path in paths[:5]:
        title, tail = body.split('\n\n', 1)
        body = title+'\n\n'+block+'\n\n'+tail
    elif path.name == 'spatial_premise_closure_audit.md':
        assert '\n## 169.' not in body
        body += '\n\n## 169. 紧支撑正势测量接续与重复研究边界\n\n'+block+next_steps
    else:
        assert '\n## 74.' not in body
        body += '\n\n## 74. 第524轮：去除硬截止后的条件性参考实现\n\n'+block
    if path == HERE/'README.md':
        body += ('\n\n## 第524轮研究索引\n\n'
            '| 轮次 | 完整报告 | 复算与核验 |\n|---|---|---|\n'
            '|524|[紧支撑正势探针](research_note_524.md)|'
            '[代码](compact_probe_reference_model.py)、[结果](compact_probe_reference_results.json)、'
            '[核验](research_round_524_checks.json)|\n')
    if path.name == 'RESEARCH_STATE.md':
        body += next_steps
    if path.name == 'research_direction.md':
        assert '最新科学轮次与检查数为523／2570' in body
        body = body.replace('最新科学轮次与检查数为523／2570', '最新科学轮次与检查数为524／2578')
        body = body.replace('完成231—523轮。', '完成231—524轮。')
        body += '\n\n**524后接续（目标不变）：** 紧支撑实际探针已记录；525候选检验动态参考读数开关的受控极限，先回查325及已有参考钟／参考几何接口，避免重证旧能量账。每个完整新研究及时形成编号报告；详见RESEARCH_STATE.md末尾。\n'
    planned[path] = body.replace('\n', newline).encode(encoding)
if not resuming:
    with (folder/'manifest.json').open('x', encoding='utf8') as stream:
        json.dump(manifest, stream, ensure_ascii=False, indent=2)
else:
    assert manifest == saved_manifest
links = 0
for path, body in planned.items():
    for link in core.link_parser()(body.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(), (path, link)
        links += 1
for path, raw in raws.items():
    assert path.read_bytes() in (raw, planned[path]), ('concurrent navigation edit', path)
for path, body in planned.items():
    if path.read_bytes() == body:
        continue
    temporary = path.with_name(path.name+'.round524.tmp')
    if temporary.exists():
        assert temporary.read_bytes() == body
    else:
        with temporary.open('xb') as stream:
            stream.write(body)
    assert path.read_bytes() == raws[path], ('concurrent navigation edit', path)
    os.replace(temporary, path)
for group in ('new_file_hashes', 'preserved_draft_hashes'):
    for name, digest in checked[group].items():
        assert core.digest(HERE/name) == digest, name
report = dict(date='2026-09-30', latest_round=524, cumulative_tests=2578,
    numbered_scientific_files=883, unique_protected_evidence_files=count,
    navigation_files=7, navigation_links=links, broken_links=0,
    navigation_snapshots_preserved=True, recovered_partial_navigation_publish=resuming,
    complete_report_published=True,
    existing_results_preserved=True, next_round=525, active_goal_unchanged=True,
    stage_complete=False, all_checks_passed=True)
with (HERE/'round524_navigation_checks.json').open('x', encoding='utf8', newline='\n') as stream:
    stream.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
print(json.dumps(report, ensure_ascii=False))
