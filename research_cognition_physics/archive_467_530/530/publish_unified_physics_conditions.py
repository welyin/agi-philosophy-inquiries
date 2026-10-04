"""Publish the user's joint-model direction, preserving history and round 530.

This registers an inventory and a familiar exact algebra example, not a new
scientific round. Run once; snapshots permit safe recovery of partial writes.
"""
from pathlib import Path
import hashlib
import json
import os

import unified_physics_condition_audit as inventory
import verify_interaction_rounds as core
import verify_round530 as evidence

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
TARGET = HERE / 'joint_conditions_navigation_checks.json'
SNAPSHOT = HERE / 'navigation_before_joint_conditions_20260930'


def main():
    assert not TARGET.exists(), 'already published; inspect the saved checks'
    checked = evidence.verify()
    assert core.read(HERE / 'research_round_530_checks.json') == checked
    audited = inventory.run()
    assert core.read(inventory.TARGET) == audited
    paths = [ROOT / 'README.md', RESEARCH / 'README.md',
             RESEARCH / 'research_direction.md', RESEARCH / 'RESEARCH_STATE.md',
             HERE / 'README.md', HERE / 'spatial_premise_closure_audit.md',
             HERE / 'three_dimensional_four_conditions_review.md']
    resuming = SNAPSHOT.exists()
    if resuming:
        manifest = core.read(SNAPSHOT / 'manifest.json')
        raws = {p: (SNAPSHOT / manifest[str(p.relative_to(ROOT))]['snapshot']).read_bytes()
                for p in paths}
        for p, raw in raws.items():
            assert hashlib.sha256(raw).hexdigest() == manifest[str(p.relative_to(ROOT))]['sha256']
    else:
        raws = {p: p.read_bytes() for p in paths}
        SNAPSHOT.mkdir(exist_ok=False)
        manifest = {}
        for p, raw in raws.items():
            key = ('root_' if p.parent == ROOT else 'research_' if p.parent == RESEARCH else 'archive_') + p.name
            with (SNAPSHOT / key).open('xb') as stream:
                stream.write(raw)
            manifest[str(p.relative_to(ROOT))] = dict(snapshot=key, sha256=hashlib.sha256(raw).hexdigest())
        with (SNAPSHOT / 'manifest.json').open('x', encoding='utf8') as stream:
            json.dump(manifest, stream, ensure_ascii=False, indent=2)

    marker = '**当前主线调整（2026-09-30，用户明确指示）：**'
    summary = (marker + ' [统一物理模型的条件总账]({p}unified_physics_condition_ledger.md)'
        '将3+1时空、引力、标准模型及其量子、记录和尺度条件纳入同一候选，'
        '先列条件、核交叉相容，再压缩并证明独立输入。三维不再是接入其它物理的前置关卡。'
        '本段与总账第7节优先于下方保留的旧“下一步”；后文530及更早成果作为科学历史继承。'
        '27项为审计入口，不是27条新认知公理；“只能联合证明”仍是待检验猜想。'
        '最新已完成科学轮次仍530／2622；[锚点候选草稿]({p}round531_drafts/STATUS.md)保留未结项，'
        '后继正式531围绕联合条件开展，不将整理和成熟结论重复计数。')
    next_steps = '''

### 2026-09-30当前接续：统一候选与联合条件

按照用户最新指示，以“3+1时空、GR和标准模型可能指向同一个答案”为研究假说，先建立共同模型及条件覆盖。猜想后文六条箭头是备选动机，不成为必须逐条证明的起点。空间维数仍是未知与验收对象，但不再要求先独立完成三维才进入物质和引力。

下一完整研究单元从统一条件总账出发，具体识别同一候选的可观测结构、态、动力学、区域组合与尺度映射，继承360—365的规范区域拼接，以及339、344的共同几何限制。优先联立组合／规范／手征物质条件与实际测量／物质／几何条件；可以先用GR＋规范物质有效理论作逆向比较目标，再核是否有更少共同结构承载这些要求。非交换几何的谱作用作为成熟对照，四维和有限谱结构等输入须显式标注；并未预先选它为认知本体论的唯一模型。

每个候选必须使用同一组几何、物质、态、参数和适用尺度，列出：原理输入、成熟定理接口、模型内推论、尚缺连接、可区分预测。不同重构路线可替代的条件按分支保留，不能机械合并。允许从既有物理反提候选认知原则，仍区分“定义相容”“存在一个实现”“唯一被迫如此”。

锚点条件来源候选保留代码、结果、审查和STATUS，因本次主线调整停止扩大该专项；不是研究目标暂停，也不签收为531。正式新研究仍须及时形成research_note_编号.md、可复算材料、结果与核验；本次条件总账及成熟反常代数复核不虚增轮次。最新已完成530／2622、902份编号科学文件、1114份既有保护证据不变。历史稿、旧路线和所有冻结结果保留；不改应用中的目标，不设置定时任务，不做图像检验。
'''
    planned = {}
    for path, raw in raws.items():
        encoding = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
        newline = '\r\n' if b'\r\n' in raw else '\n'
        body = raw.decode(encoding).replace('\r\n', '\n')
        assert marker not in body
        prefix = 'research_cognition_physics/archive_231_/' if path.parent == ROOT else 'archive_231_/' if path.parent == RESEARCH else ''
        block = summary.format(p=prefix)
        if path in paths[:5]:
            title, tail = body.split('\n\n', 1)
            body = title + '\n\n' + block + '\n\n' + tail
        elif path.name == 'spatial_premise_closure_audit.md':
            assert '\n## 176.' not in body
            body += '\n\n## 176. 用户调整主线：时空、物质与引力的联合条件\n\n' + block + next_steps
        else:
            assert '\n## 81.' not in body
            body += '\n\n## 81. 三维转为联合模型内的待解条件\n\n' + block
        if path.name in ('RESEARCH_STATE.md', 'research_direction.md'):
            body += next_steps
        if path == HERE / 'README.md':
            body += ('\n\n## 联合物理模型：条件与接续入口（非新增科学轮次）\n\n'
                '- [统一条件总账与跨部门约束](unified_physics_condition_ledger.md)\n'
                '- [条件与有理代数复核](unified_physics_condition_audit.py)、[保存结果](unified_physics_condition_audit_results.json)\n'
                '- [未结项锚点分支状态](round531_drafts/STATUS.md)\n')
        planned[path] = body.replace('\n', newline).encode(encoding)

    links = 0
    for path, body in planned.items():
        for link in core.link_parser()(body.decode('utf-8-sig')):
            assert (path.parent / link).resolve().exists(), (path, link)
            links += 1
    for path, raw in raws.items():
        assert path.read_bytes() in (raw, planned[path]), ('concurrent edit', path)
    for path, body in planned.items():
        if path.read_bytes() == body:
            continue
        temp = path.with_name(path.name + '.jointconditions.tmp')
        if temp.exists():
            assert temp.read_bytes() == body
        else:
            with temp.open('xb') as stream:
                stream.write(body)
        assert path.read_bytes() == raws[path], ('concurrent edit', path)
        os.replace(temp, path)
    for group in ('new_file_hashes', 'preserved_draft_hashes'):
        for name, sha in checked[group].items():
            assert core.digest(HERE / name) == sha, name
    assert inventory.run() == audited
    names = ('unified_physics_condition_ledger.md', 'unified_physics_condition_audit.py',
             'unified_physics_condition_audit_results.json', 'round531_drafts/STATUS.md',
             Path(__file__).name)
    report = dict(date='2026-09-30', kind='user_directed_joint_model_research_entry',
        latest_completed_round=530, cumulative_tests_unchanged=2622,
        numbered_scientific_files_unchanged=902, historical_protected_evidence_files=1114,
        complete_frozen_chain_verified=True, inventory_and_exact_algebra_verified=True,
        conditions_count=27, navigation_files=7, navigation_links=links, broken_links=0,
        navigation_snapshots_preserved=True, recovered_partial_publish=resuming,
        artifact_hashes={name: core.digest(HERE / name) for name in names},
        navigation_hashes={str(p.relative_to(ROOT)): core.digest(p) for p in paths},
        unfinished_anchor_candidate_preserved=True, next_scientific_round=531,
        application_goal_unchanged=True, stage_complete=False, all_checks_passed=True)
    with TARGET.open('x', encoding='utf8', newline='\n') as stream:
        stream.write(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ('latest_completed_round', 'conditions_count',
        'navigation_files', 'navigation_links', 'broken_links', 'all_checks_passed')}))


if __name__ == '__main__':
    main()
