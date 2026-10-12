"""Maintain one reversible live-navigation layer for the active field stage."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
from organize_research_1044_1062 import text_edits
from research_layout import relocation_original_bytes

BASE = ROOT / 'research_cognition_physics'
LAYER = BASE / '_migration/field_active_1120/layer.json'
NAMES = ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')
NOTICE = '''## 1120—1156阶段：按现行公理严格不足的限定分支结项

**最新正式轮次1156；37轮、144组构造/代数复算，新采用公理0，科学校准新增0、历史累计3860。** [阶段结项](archive_1120_/_shared/stage_closure/stage_completion.md)与[独立范围审计](archive_1120_/_shared/stage_closure/closure_review.md)已完成，研究档案保留在[archive_1120_](archive_1120_/README.md)。

**已证明的限定结论：** 在现行采用的有限有效类型／任务语义下，F＋U＋C＋P、CO、SR及D-P不足以保证每项合法过程具有同事件、同期限的因果局域量子实现。[1156](archive_1120_/research_note_1156.md)给完整正常量子资源、原钟/端点、全部许可接续与完整四P的反模型；同一均匀四设置实验的联合结果误差下界为(2−√2)/4≈0.1464466。

这不是局域模型不存在、1155特定源任务失败或整个认知纲领被否定；也没有签收局域场、完整共同局部来源或引力的正向构造。数学证明、完整环境权限、过去集因果接续和四组复算均已审查；不把成熟PR/CHSH工具当成新的实验发现。

1120—1155的正常来源、实际反冲、原钟与共同运动、区域交换等条件成果按各自原合同保留。其中[1155](archive_1120_/research_note_1155.md)的分站四交换属于明确记录源，不自动等于自治全装置四P。[旧条件定理连接账](archive_1120_/_shared/stage_closure/conditional_links.md)继续复用1015—1019、1027—1029；1029已消去的一阶提升不重新列为缺口。

[下一阶段候选](archive_1120_/_shared/stage_closure/next_stage_candidates.md)：优先研究能排除该反模型、保留纠缠及完整资源的局域组合条件及其认知含义，再检查共同有效实现。尚未采用任何新公理，未自动创建新目标、任务或定时任务；1157不分配。

[轮次索引](archive_1120_/文件索引.md) · [正式依赖与验收账](archive_1120_/_shared/dependency_ledger.json) · [原目标](archive_1120_/_shared/goal.md) · [结项收据](archive_1120_/_shared/stage_closure/completion_record.json) · [历史只读复算入口](_migration/field_active_1120/README.md)。

以下“当前／active／下一步”保留历史研究时点，不覆盖本栏。前阶段1106—1119和五篇编号论文保持原结项范围；完整物理统一目标尚未完成。

---

'''


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    frozen_path = BASE / '_migration/merge_1106_1119/manifest.json'
    frozen = json.loads(frozen_path.read_text(encoding='utf-8'))
    prior = {e['destination']: e for e in frozen['entries']}
    existing_bytes = LAYER.read_bytes() if LAYER.exists() else None
    existing = {e['original']: e for e in json.loads(existing_bytes)['entries']} if existing_bytes else {}
    entries, changes = [], []
    for name in NAMES:
        path = BASE / name
        relative = path.relative_to(ROOT).as_posix()
        current = path.read_bytes()
        original = relocation_original_bytes(current, existing[relative]) if existing else current
        assert sha(original) == prior[relative]['current_sha256']
        assert len(original) == prior[relative]['current_bytes']
        encoding = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf-8'
        before = original.decode(encoding)
        newline = '\r\n' if '\r\n' in before else '\n'
        split = before.index('\n') + 1
        after = before[:split] + newline + NOTICE.replace('\n', newline) + before[split:]
        new = after.encode(encoding)
        entry = dict(original=relative, destination=relative,
                     original_sha256=sha(original), current_sha256=sha(new),
                     original_bytes=len(original), current_bytes=len(new),
                     edits=text_edits(before, after))
        assert relocation_original_bytes(new, entry) == original
        entries.append(entry)
        changes.append((path, current, new))
    support_names = ['scripts/run_research_merged.py', frozen_path.relative_to(ROOT).as_posix()]
    layer = dict(schema='field_active_navigation_v1', date='2026-10-12', entries=entries,
                 unchanged_support_sha256={name: sha((ROOT / name).read_bytes()) for name in support_names})
    if args.apply:
        assert (LAYER.read_bytes() if LAYER.exists() else None) == existing_bytes
        for path, before, after in changes:
            assert path.read_bytes() == before
        LAYER.parent.mkdir(parents=True, exist_ok=True)
        LAYER.write_text(json.dumps(layer, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        for path, before, after in changes:
            if before != after:
                path.write_bytes(after)
    print(json.dumps(dict(status='APPLIED' if args.apply else 'CHECKED', files=4,
                          frozen_merge_unchanged=True, original_navigation_recoverable=True)))


if __name__ == '__main__':
    main()
