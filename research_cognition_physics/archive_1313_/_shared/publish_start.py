"""One-time publication of the reviewed goal; preserve old navigation bytes."""
from pathlib import Path
import hashlib
import json
import os
import sys

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE.parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from research_layout import relocation_original_bytes
import run_research_conditional_selection_closed as closed


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def atomic_write(path, raw):
    temp = path.with_name(path.name + '.interaction-1313.tmp')
    with temp.open('xb') as f:
        f.write(raw)
    try:
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()


def main():
    layer_path = ROOT / 'research_cognition_physics/_migration/interaction_generation_active_1313/layer.json'
    receipt_path = STAGE / '_shared/navigation_start_receipt.json'
    assert not layer_path.exists(), 'Already published; do not repeat'
    assert not receipt_path.exists(), 'Publication receipt already exists'
    receipt = json.loads((STAGE / '_shared/goal_start_receipt.json').read_text(encoding='utf8'))
    assert receipt['status'] == 'active'
    prompt = (STAGE / '_shared/goal_prompt.txt').read_text(encoding='utf8').strip()
    assert receipt['objective'] == prompt
    review = (STAGE / '_shared/goal_scope_review.md').read_text(encoding='utf8')
    assert 'PASS_WITH_SCOPE' in review
    for name in ('goal.md', 'start_plan.md'):
        assert sha((STAGE / '_shared' / name).read_bytes()) in review, name
    store = closed.ClosedConditionalSelectionLayout()
    closure_check = store.verify_closure()
    store.verify_conditional_navigation()
    prefix = ('\n## 当前目标：从认知交互机制生成物质分工、代际与质量结构\n\n'
              '2026-10-18，依用户授权启动[1313起新阶段](archive_1313_/README.md)。'
              '[正式目标与联合验收](archive_1313_/_shared/goal.md) · '
              '[启动计划](archive_1313_/_shared/start_plan.md) · '
              '[目标范围独审](archive_1313_/_shared/goal_scope_review.md)。应用目标已开启。\n\n'
              '长期目标为完整标准模型及现实参数。本阶段先检验同一内部交互如何约束物质角色、荷与代际，'
              '并生成质量、混合和耦合关系。首项比较单参考／右单态、多参考和带活动荷右伙伴的实际交接；'
              '保持共同逻辑身份、未知伙伴关联、全部结果与完整资源。'
              '不重复一般SWAP、单块质量秩或1311已经退化的预算路线。\n\n'
              '[计0事前筛查](archive_1313_/_admission/finite_handoff_role_screen_after1312.md)已回查1306等成果：'
              '下一项先明确旧参考真实更新／退出后的完整接收通道和合格竞争者，只有新命题才发布1313报告。'
              '最新验收科学成果仍为1312；目标设定与归档不计新轮，新公理0。\n\n'
              '上一阶段已归档为[archive_1308_1312](archive_1308_1312/README.md)，'
              '[论文14](14_单参考质量与有限交换任务下的规范结构选择_阶段论文.md)和全部证据保留。'
              '[迁移与复算](_migration/closure_1308_1312/README.md)保持历史核验。'
              '以下各时点记录保留，当前安排以本节和新阶段目标为准。\n')
    entries, updates = [], []
    for name in ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md'):
        rel = 'research_cognition_physics/' + name
        path = ROOT / rel
        raw = path.read_bytes()
        boundary = store.conditional_closure_entries[rel]
        assert sha(raw) == boundary['current_sha256']
        enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
        text = raw.decode(enc)
        at = text.index('\n') + 1
        new = (text[:at] + prefix + text[at:]).encode(enc)
        entry = dict(original=rel, destination=rel, original_sha256=sha(raw),
                     current_sha256=sha(new), original_bytes=len(raw), current_bytes=len(new),
                     edits=[dict(start=at, end=at, before='', after=prefix)])
        assert relocation_original_bytes(new, entry) == raw
        entries.append(entry)
        updates.append((path, raw, new))
    support = dict(store.conditional_closure_manifest['unchanged_support_sha256'])
    for path in (closed.MANIFEST, closed.THIS):
        support[path.relative_to(ROOT).as_posix()] = sha(path.read_bytes())
    layer = dict(schema='interaction_generation_active_navigation_v1', date='2026-10-18',
                 stage_path=STAGE.relative_to(ROOT).as_posix(), stage_first_round=1313,
                 entry_points=[STAGE.relative_to(ROOT).as_posix() + '/' + p
                               for p in ('README.md', '_shared/goal.md')],
                 latest_accepted_round=1312, in_progress_round=None, next_round=1313,
                 application_goal_status_last_observed='active', stage_complete=False,
                 startup_state='count-zero admission screening; no new scientific acceptance',
                 entries=entries, unchanged_support_sha256=support)
    layer_bytes = (json.dumps(layer, ensure_ascii=False, indent=2) + '\n').encode('utf8')
    result = dict(status='PUBLISHED', scientific_rounds_added=0, new_axioms=0,
                  stage_first_round=1313, latest_accepted_round=1312,
                  closure_before_navigation=closure_check, navigation_files=4,
                  layer_sha256=sha(layer_bytes),
                  source_sha256={p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in
                                 [STAGE / '_shared/goal.md', STAGE / '_shared/start_plan.md',
                                  STAGE / '_shared/goal_scope_review.md',
                                  STAGE / '_shared/goal_prompt.txt',
                                  STAGE / '_shared/goal_start_receipt.json']})
    receipt_bytes = (json.dumps(result, ensure_ascii=False, indent=2) + '\n').encode('utf8')
    for path, raw, new in updates:
        assert path.read_bytes() == raw, 'Concurrent navigation edit: ' + str(path)
    assert sha(closed.MANIFEST.read_bytes()) == store.conditional_closure_sha256
    layer_path.parent.mkdir(parents=True, exist_ok=True)
    written = []
    try:
        for path, raw, new in updates:
            assert path.read_bytes() == raw, 'Concurrent navigation edit: ' + str(path)
            atomic_write(path, new)
            written.append((path, raw, new))
        assert not layer_path.exists() and not receipt_path.exists()
        atomic_write(layer_path, layer_bytes)
        atomic_write(receipt_path, receipt_bytes)
    except BaseException:
        for path, expected in ((receipt_path, receipt_bytes), (layer_path, layer_bytes)):
            if path.exists() and path.read_bytes() == expected:
                path.unlink()
        for path, raw, new in reversed(written):
            if path.read_bytes() == new:
                atomic_write(path, raw)
        raise
    print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
