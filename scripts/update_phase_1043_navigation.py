"""Update live navigation after closing phase 1009--1043; no scientific edits."""
from pathlib import Path
import argparse
import hashlib
import json

from organize_research_231_775 import links

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'research_cognition_physics'
OLD = 'archive_1009_/'
NEW = 'archive_1009_1043/'
PAPER = '条件生成链与物理选择的剩余自由_阶段论文.md'

TOPS = {
    'README.md': ('## 1042：', '''# 认知物理研究：阶段论文与研究档案

## 当前入口：1009—1043已结项，1044起的下一阶段已准备

[阶段论文](条件生成链与物理选择的剩余自由_阶段论文.md)汇总同一领先作用下的规范—物质—来源—引力条件生成链，并区分各层已证自由与未决输入。[1009—1043档案](archive_1009_1043/README.md)保留全部35份正式报告及配套材料；[1043综合报告](archive_1009_1043/research_note_1043.md)和[逐项验收](archive_1009_1043/1043/completion_audit.md)是阶段科学结项依据。

**下一阶段准备研究：认知操作要求对物理规律有没有独立选择力。** [1044起的准备计划](archive_1044_/README.md)规定首项共同量子物理过程检验、验收及停止条件。最新完成轮次仍为1043，累计3819；准备目录不代表1044已完成或新目标已启动。

阶段完成不等于完整物理由认知生成。34条依赖分别保留共同对象与前提，14类输入和16类现象的未决部分公开。33个原默认入口通过，1009的科学及固定旧资产另经[范围复核](archive_1009_1043/1043/legacy_1009_scope_results.json)通过；原全目录快照差异保留。

[原最终验收](archive_1009_1043/1043/research_round_1043_checks.json) · [归档维护与复算入口](_migration/closure_1009_1043_20261008/README.md) · [研究方向](research_direction.md) · [研究状态](RESEARCH_STATE.md)。以下保留各轮写入时的历史状态，旧active或旧“下一步”不覆盖本栏。

'''),
    'research_direction.md': ('## 本阶段原目标正文', '''# 研究方向：认知本体论与现代物理的共同模型

## 当前安排：关闭1009—1043阶段，准备检验独立选择力

2026-10-08按用户要求收尾。条件生成链、适用前提和剩余自由分类已经[逐项验收](archive_1009_1043/1043/completion_audit.md)；[阶段论文](条件生成链与物理选择的剩余自由_阶段论文.md)统一呈现成果，全部编号材料保留于[archive_1009_1043](archive_1009_1043/README.md)。本次整理不新增科学结果或认知公理，累计3819保持。

### 下一阶段建议主问题

**同一量子物理过程承担主体、记录、交互、资源、来源与跨尺度描述时，认知操作要求能否排除仅凭物理一致性仍然允许的自然动力学选择？**

先选择一项有真实物理含义的耦合／作用自由，在同一有限任务合同内检验其量子过程、实际仪器后态、物质来源和已声明反馈；优先复用已有共同接口。结果可以是新的必要关系，也可以是明确前提下的可辨自由；不能把未知写成独立，也不能把参数变化后预测不同本身当成认知选择。

准备区为[archive_1044_](archive_1044_/README.md)。其中的首项准入、判别标准及停止条件约束后续投入：先验证这个问题是否比929／1031旧证书真正更强，再补决定方向的桥；不将完整E_nat、全部器件或UV完成设为默认门槛。整体量子—时空—物质—引力目标保持，单项测试不代表全部物理已被选择或否定。

旧空间382—384、386或425、522—523以及1042区域恢复工具继续按原前提复用，按其能否解决共同选择问题决定接入顺序。允许双向推导和提出新原则，但须明列额外内容、避免循环并独立检验。

**当前仅完成下一阶段准备；没有创建新应用目标、任务或定时任务，1044没有正式研究报告。** 下文保存1009阶段原目标与历次方向记录，旧active及旧接续建议仅为当时状态。

'''),
    'RESEARCH_STATE.md': ('## 1042时点状态', '''# 研究状态

## 当前状态：1009—1043已结项，1044起的下一阶段已准备（2026-10-08）

- 最新完成轮次：**1043**；累计科学校准：**3819**。本次阶段收尾不增加轮次、科学组或认知公理。
- 本阶段科学结论：[1043综合报告](archive_1009_1043/research_note_1043.md)、[逐项验收](archive_1009_1043/1043/completion_audit.md)及[原最终验收](archive_1009_1043/1043/research_round_1043_checks.json)。[阶段论文](条件生成链与物理选择的剩余自由_阶段论文.md)是据此整理的阅读入口，没有扩大原结论。
- 全部35份编号报告与配套材料归入[archive_1009_1043](archive_1009_1043/README.md)。[归档维护记录](_migration/closure_1009_1043_20261008/README.md)交代路径变化及复算方式；科学代码、结果和历史验收不倒改。
- 应用中上一目标已完成：[结项回执](archive_1009_1043/1043/goal_completion.json)。本次核对没有进行中的目标，未自动新建目标。
- 下一阶段：[archive_1044_准备计划](archive_1044_/README.md)，检验认知操作要求对自然动力学的独立选择力。**尚未开始1044正式研究**，没有新科学结论；下次从准备计划的准入检验接续。

旧阶段的验收口径保持：33个原默认入口通过；1009原入口的出版范围比较失败原样保留，科学及7个固定旧资产经独立诊断通过，358份历史文件在当时回归前后相同。归档后的维护核验与这些历史结果分开记录。

同一领先作用下的顶点、动态来源与引力必要关系已经建立；929、1031、1042的自由证据分别保留其合同，不能合并为全部认知公理不蕴涵现实物理的证明。完整认知物理生成、全E_nat及现实全部共同预测仍未证明。

以下均为历史时点。旧active、待审、额度与“下一项”记录不覆盖本栏。

'''),
}

def update():
    assert (BASE / 'archive_1009_1043/1043/research_round_1043_checks.json').is_file()
    assert (BASE / PAPER).is_file()
    assert (BASE / 'archive_1044_/README.md').is_file()
    for name, (marker, top) in TOPS.items():
        path = BASE / name
        body = path.read_text('utf8')
        assert marker in body, name
        body = top + body[body.index(marker):]
        body = body.replace(OLD, NEW)
        path.write_text(body, encoding='utf8', newline='\n')
    paper = BASE / PAPER
    body = paper.read_text('utf8').replace(OLD, NEW)
    paper.write_text(body, encoding='utf8', newline='\n')
    print('Updated three live navigation files and stage-paper links.')


def check(record=False):
    migration = BASE / '_migration/closure_1009_1043_20261008'
    receipt = migration / 'handoff_checks.json'
    paths = [BASE/p for p in (*TOPS, PAPER)]
    paths.extend(sorted((BASE/'archive_1044_').glob('*.md')))
    missing, count = [], 0
    for path in paths:
        for _, _, target, local in links(path.read_text('utf-8-sig')):
            count += 1
            actual = (path.parent/local.replace('\\', '/')).resolve()
            if not actual.exists():
                missing.append(dict(source=str(path.relative_to(ROOT)), target=target))
    assert not missing, missing
    assert not (BASE/'archive_1009_').exists()
    notes = sorted((BASE/'archive_1009_1043').glob('research_note_*.md'))
    assert [int(p.stem.split('_')[-1]) for p in notes] == list(range(1009, 1044))
    assert not list((BASE/'archive_1044_').glob('research_note_*.md'))
    maintenance = migration/'checks.json'
    accepted = json.loads(maintenance.read_text('utf8'))
    assert accepted['all_checks_passed'] and accepted['cumulative_unchanged'] == 3819
    value = dict(schema='stage_closure_handoff_v1', date='2026-10-08',
                 all_checks_passed=True, checked_documents=len(paths),
                 local_links_checked=count, missing_local_links=missing,
                 closed_stage_formal_reports=35, latest_completed_round=1043,
                 next_stage_preparation_only=True, new_scientific_groups=0,
                 cumulative_unchanged=3819,
                 migration_checks_sha256=hashlib.sha256(maintenance.read_bytes()).hexdigest(),
                 document_sha256_at_handoff={p.relative_to(ROOT).as_posix():
                                            hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    if record:
        with receipt.open('x', encoding='utf8') as stream:
            json.dump(value, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    # Live navigation is allowed to evolve. Its recorded hashes are a dated
    # handoff record, never a new immutable scientific prerequisite.
    print(json.dumps({k:v for k,v in value.items() if k != 'document_sha256_at_handoff'}, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    if args.record and not args.check:
        parser.error('--record requires --check')
    check(args.record) if args.check else update()


if __name__ == '__main__':
    main()
