"""Apply one minimal, reversible four-document navigation update."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
from organize_research_1044_1062 import text_edits
from research_layout import relocation_original_bytes

MAINTENANCE = Path(__file__).resolve().parent
LAYER = MAINTENANCE / 'layer.json'
NAMES = ['README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md']
SECTION = '''## 当前目标：实际狭义相对论的共同连接（1106起）

**新目标已按用户授权开启，工具状态active。** 工作目录为[archive_1106_](archive_1106_/README.md)，完整范围与验收约束见[本阶段目标](archive_1106_/_shared/goal.md)。1086—1105已封存，阶段论文及历史研究不改。

继承F＋U＋C＋P四原则、CO1—CO5认知原子操作五公理与已采用的SR-T／A／R／O／F；SR-S只在明确启用的可选定标分支使用，生成完备G仍未采用。继续区别原子公理、组合权限与宏观读出，不将新输入追认为旧原则的推论。

推进次序：先建立实际钟与世界线的共同连接，再检查同一参考族的实际实现、全部可辨影响及物质交互的一致性。每阶段记录依赖并收敛冗余，保持实际有限尺度验收；阶段运动学结论尚不等于完整物理SR或全部理论完成。

**下一轮1106处于启动准入，尚无新的正式完成报告或科学校准。** 最新完成轮次仍为1105，累计科学校准3860。历史复算使用新的最外层入口 `scripts/run_research_motion_active.py`，冻结代码、结果与清单不改。

以下结项栏及历史“当前／下一步”保留其原研究时点；本栏和本阶段目标给出新的执行方向。

---

'''


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    if LAYER.exists():
        raise ValueError('Navigation layer already exists; do not apply twice')
    entries, changes = [], []
    for name in NAMES:
        path = ROOT / 'research_cognition_physics' / name
        raw = path.read_bytes()
        encoding = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
        before = raw.decode(encoding)
        if '## 当前目标：实际狭义相对论的共同连接（1106起）' in before:
            raise ValueError('New navigation already exists')
        match = re.match(r'(#[^\r\n]*(?:\r\n|\n))', before)
        if not match:
            raise ValueError('Expected a root document title')
        newline = '\r\n' if '\r\n' in before else '\n'
        after = match.group(1) + newline + SECTION.replace('\n', newline) + before[match.end():].lstrip('\r\n')
        updated = after.encode(encoding)
        relative = path.relative_to(ROOT).as_posix()
        entry = dict(original=relative, destination=relative, original_sha256=sha(raw),
                     current_sha256=sha(updated), original_bytes=len(raw),
                     current_bytes=len(updated), edits=text_edits(before, after))
        if relocation_original_bytes(updated, entry) != raw:
            raise ValueError('Inverse edit does not recover original navigation')
        entries.append(entry)
        changes.append((path, raw, updated))
    support = ['scripts/run_research_sr_archive.py',
               'research_cognition_physics/_migration/sr_closure_1086_1105/manifest.json',
               'research_cognition_physics/archive_1086_1105/_shared/stage_paper_closure/publication_checks.json',
               'research_cognition_physics/认知时空操作公理与洛伦兹运动学_阶段论文.md']
    data = dict(schema='motion_active_navigation_v1', date='2026-10-10', entries=entries,
                unchanged_support_sha256={name: sha((ROOT / name).read_bytes()) for name in support},
                scope='Four root navigation insertions only; new stage starts at 1106',
                actual_goal_status='active', new_formal_reports=0, scientific_calibration_increment=0)
    for path, raw, _ in changes:
        if path.read_bytes() != raw:
            raise ValueError('Navigation changed concurrently')
    with LAYER.open('x', encoding='utf8', newline='\n') as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
    for path, _, updated in changes:
        path.write_bytes(updated)
    print(json.dumps(dict(files=4, reversible=True, old_scientific_assets_modified=False), ensure_ascii=False))


if __name__ == '__main__':
    main()
