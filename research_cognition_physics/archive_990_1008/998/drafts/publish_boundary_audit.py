"""Publish the working audit without changing formal round/test counts."""
from pathlib import Path
import hashlib
import json
import os

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1]
RESEARCH = STAGE.parent
ROOT = RESEARCH.parent
TITLE = '## 998工作审计：共同边界与候选投入'
MARK = '## 997：有效Λ与实际可达范围'
NAV = [RESEARCH/n for n in ('README.md', 'research_direction.md', 'RESEARCH_STATE.md')]
NAV += [STAGE/n for n in ('README.md', '文件索引.md', '阶段成果总览.md',
                          '跨阶段主题索引.md', '_shared/notes/unified_physics_condition_ledger_current.md')]


def link(doc, target):
    return Path(os.path.relpath(target, doc.parent)).as_posix()


def main():
    updates = []
    for doc in NAV:
        raw = doc.read_bytes()
        bom = raw.startswith(b'\xef\xbb\xbf')
        content = raw.decode('utf-8-sig').replace('\r\n', '\n')
        assert TITLE not in content and content.count(MARK) == 1, doc
        note = link(doc, STAGE/'research_note_998_working.md')
        results = link(doc, STAGE/'998/boundary_scope_audit_results.json')
        checks = link(doc, HERE/'boundary_adoption_audit_checks.json')
        nxt = link(doc, HERE/'NEXT.md')
        block = (TITLE + '\n\n'
            f'[998工作报告]({note})直接复用旧约束、来源和热资源结论，明确低关联不能删掉物理约束关联。'
            '997加速校准起点的辐射份额约2.06%，不能字面接成996辐射主导窗口；这不反证两轮各自结果。'
            f'[复算结果]({results}) · [审计核验]({checks})。\n\n')
        block += ('本次为采用审计，正式997／累计3781保持，不新增科学试验组。当前应用目标已含方向优先与有限有效域要求；'
            f'接[下一步]({nxt})先审稳定组织、资源和宇宙阶段的共同机制，停止CP、Λ及器件调参。'
            '完整共同模型尚未完成，原998入口和全部历史保留。\n\n')
        content = content.replace(MARK, block + MARK, 1)
        if doc.name == 'research_direction.md':
            anchor = '### 目标：以整体机制和决定性检验推进认知操作与物理的共同模型'
            assert content.count(anchor) == 1
            decision = ('**本次执行次序已落实（998工作期）：** 先判断候选是否值得承担共同模型，再决定哪些有限预测要签收；'
                '只影响另选UV外推的问题单独登记。当前先审共同边界及稳定组织／资源机制，不续修单一CP、Λ或器件候选。'
                '同一理论的不同准备与同一宇宙历史分开验收。详见[工作报告](archive_764_/research_note_998_working.md)。'
                '此处更新执行安排；应用目标已包含该方向，本次没有另改正文。\n\n')
            content = content.replace(anchor, decision + anchor, 1)
        if doc.name == 'unified_physics_condition_ledger_current.md':
            anchor = '### 当前取舍\n'
            assert content.count(anchor) == 1
            content = content.replace(anchor, anchor + '\n'
                '998工作审计优先：采用959所需约束关联及973共同状态／来源边界；允许合法准备输入，'
                '不要求先生成全部宇宙过去。996与997尚非同一热史窗口；先审稳定组织与资源链的共同机制，'
                '不继续修CP、Λ或热器件。C01—C27状态不因本次整理被标为全部关闭；正式997保持。\n', 1)
        newline = '\r\n' if b'\r\n' in raw else '\n'
        changed = content.replace('\n', newline).encode('utf-8')
        if bom:
            changed = b'\xef\xbb\xbf' + changed
        updates.append((doc, raw, changed))
    # Check all originals immediately before writing; preserve unrelated edits.
    for doc, raw, _ in updates:
        assert doc.read_bytes() == raw, f'Concurrent edit: {doc}'
    before = {str(p.relative_to(ROOT)): hashlib.sha256(raw).hexdigest()
              for p,raw,_ in updates}
    for doc, _, changed in updates:
        doc.write_bytes(changed)
    record = dict(date='2026-10-07', kind='working_audit_navigation_update',
        formal_reports=997, new_scientific_test_groups=0,
        application_goal_rewritten=False, navigation_before_hashes=before)
    with (HERE/'publication_record.json').open('x', encoding='utf-8') as dest:
        json.dump(record, dest, ensure_ascii=False, indent=2)
        dest.write('\n')
    print('Updated 8 live navigation files; formal reports and frozen files preserved.')


if __name__ == '__main__':
    main()
