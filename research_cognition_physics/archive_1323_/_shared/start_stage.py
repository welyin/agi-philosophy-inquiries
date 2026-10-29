"""Publish the reviewed round-1323 start once; default mode only checks readiness."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import sys
import uuid

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE.parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True
from research_layout import relocation_original_bytes
import run_research_interaction_generation_closed as closed

LAYER_PATH = ROOT / 'research_cognition_physics/_migration/interaction_rules_active_1323/layer.json'
RECEIPT_PATH = STAGE / '_shared/navigation_start_receipt.json'
LOCK_PATH = STAGE / '_shared/navigation_start.lock'
RUNNER = ROOT / 'scripts/run_research_interaction_rules_active.py'
SCREEN = STAGE / '_admission/reference_access_screen_after1322.md'
ROOT_NAMES = ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def atomic_replace(path, raw):
    """Replace an existing navigation file; never overwrite another temporary file."""
    temp = path.with_name(path.name + '.interaction-rules-1323.tmp')
    created = False
    try:
        with temp.open('xb') as stream:
            created = True
            stream.write(raw)
        os.replace(temp, path)
    finally:
        if created and temp.exists() and raw.startswith(temp.read_bytes()):
            temp.unlink()


@contextmanager
def publication_lock():
    """Serialize preparation, writes and rollback; never remove another owner's lock."""
    token = (str(uuid.uuid4()) + '\n' + str(os.getpid()) + '\n').encode('ascii')
    acquired = False
    try:
        with LOCK_PATH.open('xb') as stream:
            acquired = True
            stream.write(token)
            stream.flush()
            os.fsync(stream.fileno())
        yield
    finally:
        if acquired and LOCK_PATH.exists() and LOCK_PATH.read_bytes() == token:
            LOCK_PATH.unlink()


def prepare():
    require(not LAYER_PATH.exists(), 'Already published; do not repeat')
    require(not RECEIPT_PATH.exists(), 'Publication receipt already exists')
    stage_name = STAGE.relative_to(ROOT).as_posix()
    require(stage_name == 'research_cognition_physics/archive_1323_', 'Unexpected stage')
    goal_path = STAGE / '_shared/goal.md'
    review_path = STAGE / '_shared/goal_scope_review.md'
    goal_receipt_path = STAGE / '_shared/goal_start_receipt.json'
    goal_receipt = json.loads(goal_receipt_path.read_text(encoding='utf8'))
    require(goal_receipt['status'] == 'active', 'Application goal is not recorded active')
    require(goal_receipt['stage_path'] == stage_name, 'Goal receipt stage differs')
    require(goal_receipt['goal_document'] == '_shared/goal.md', 'Unexpected goal document')
    require(goal_receipt['title'] == '基本相互作用规则的认知来源', 'Goal title differs')
    require(goal_receipt['first_available_scientific_round'] == 1323
            and goal_receipt['latest_scientific_acceptance_at_start'] == 1322,
            'Goal receipt scientific boundary differs')
    require(goal_receipt['new_scientific_rounds'] == 0 and goal_receipt['new_axioms'] == 0,
            'Startup must add no scientific result or axiom')
    require(sha(goal_path.read_bytes()) == goal_receipt['goal_document_sha256'],
            'Goal differs from successful goal-creation receipt')
    require(sha(review_path.read_bytes()) == goal_receipt['scope_review_sha256'],
            'Scope review differs from goal-creation receipt')
    review = review_path.read_text(encoding='utf8')
    require('结论：**PASS**' in review, 'Goal scope review has not passed')
    for name in ('goal.md', 'start_plan.md', 'selection_ledger.md'):
        require(sha((STAGE / '_shared' / name).read_bytes()) in review,
                'Scope review does not bind ' + name)
    require(SCREEN.is_file(), 'Count-zero screening entry is missing')
    require((STAGE / 'README.md').is_file() and RUNNER.is_file(), 'Missing startup entry')
    require('计0' in SCREEN.read_text(encoding='utf8'), 'Screening status is not explicit')

    # Construction is intentionally delayed until the separate closure exists.
    store = closed.ClosedInteractionGenerationLayout()
    closure_check = store.verify_closure()
    prior_navigation_check = store.verify_interaction_navigation()
    require(closure_check['status'] == 'PASS' and prior_navigation_check['status'] == 'PASS',
            'Previous closure or recovered navigation failed')
    prefix = (
        '\n## 当前目标：基本相互作用规则的认知来源\n\n'
        '2026-10-19，依用户授权启动[1323起新阶段](archive_1323_/README.md)。'
        '[正式目标与联合验收](archive_1323_/_shared/goal.md) · '
        '[启动计划](archive_1323_/_shared/start_plan.md) · '
        '[目标范围独审](archive_1323_/_shared/goal_scope_review.md)。应用目标已开启。\n\n'
        '长期目标仍为完整标准模型及现实参数。本阶段沿“实际参考访问合同—'
        '基本相互作用与组合闭包—物质角色及代际／参数影响”推进，'
        '以同一父过程、完整来源与有限误差作联合核验。'
        '四原则与五公理共同保留；一次原子事件、一份参考和一次场插入不能直接等同。'
        '不先规定线性Higgs、右弱单态或恰三代，不因首个弱机制反例就结束整个目标。\n\n'
        '[计0事前筛查](archive_1323_/_admission/reference_access_screen_after1322.md)'
        '已回查原子操作、组合层及参考访问的旧成果。下一项先明确独立操作机制与选择力，'
        '只有通过准入的新命题才发布1323报告。'
        '最新验收科学成果仍为1322；目标设定、归档和筛查不计新轮，新公理0。\n\n'
        '上一阶段已归档为[archive_1313_1322](archive_1313_1322/README.md)，'
        '[论文15](15_认知交接机制的物质角色选择与参数剩余自由_阶段论文.md)和全部证据保留。'
        '[迁移与复算](_migration/closure_1313_1322/README.md)保持历史核验。'
        '以下各时点记录保留，当前安排以本节和新阶段目标为准。\n')
    entries, updates = [], []
    for name in ROOT_NAMES:
        rel = 'research_cognition_physics/' + name
        path = ROOT / rel
        raw = path.read_bytes()
        boundary = store.interaction_closure_entries[rel]
        require((sha(raw), len(raw)) == (boundary['current_sha256'], boundary['current_bytes']),
                'Navigation does not match closure boundary: ' + rel)
        enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
        text = raw.decode(enc)
        at = text.index('\n') + 1
        new = (text[:at] + prefix + text[at:]).encode(enc)
        entry = dict(original=rel, destination=rel, original_sha256=sha(raw),
                     current_sha256=sha(new), original_bytes=len(raw), current_bytes=len(new),
                     edits=[dict(start=at, end=at, before='', after=prefix)])
        require(relocation_original_bytes(new, entry) == raw,
                'Navigation is not exactly reversible: ' + rel)
        entries.append(entry)
        updates.append((path, raw, new))
    support = dict(store.interaction_closure_manifest['unchanged_support_sha256'])
    for path in (closed.MANIFEST, closed.THIS, RUNNER):
        support[path.relative_to(ROOT).as_posix()] = sha(path.read_bytes())
    layer = dict(
        schema='interaction_rules_active_navigation_v1', date='2026-10-19',
        stage_path=stage_name, stage_first_round=1323,
        entry_points=[stage_name + '/' + p for p in ('README.md', '_shared/goal.md')],
        latest_accepted_round=1322, in_progress_round=None, next_round=1323,
        application_goal_status_last_observed='active', stage_complete=False,
        startup_state='count-zero admission screening; no new scientific acceptance',
        entries=entries, unchanged_support_sha256=support)
    layer_bytes = (json.dumps(layer, ensure_ascii=False, indent=2) + '\n').encode('utf8')
    source_paths = [
        goal_path, STAGE / '_shared/start_plan.md', STAGE / '_shared/selection_ledger.md',
        review_path, goal_receipt_path, STAGE / 'README.md', SCREEN, RUNNER,
        Path(__file__).resolve()]
    result = dict(
        status='PUBLISHED', scientific_rounds_added=0, new_axioms=0,
        stage_first_round=1323, latest_accepted_round=1322,
        application_goal_status_last_observed='active',
        application_goal_created_at=goal_receipt['created_at'],
        closure_before_navigation=closure_check,
        recovered_previous_navigation=prior_navigation_check,
        navigation_files=4, layer_sha256=sha(layer_bytes),
        original_navigation_recoverable=True, scientific_python_json_modified=False,
        source_sha256={p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in source_paths})
    receipt_bytes = (json.dumps(result, ensure_ascii=False, indent=2) + '\n').encode('utf8')
    return updates, layer_bytes, receipt_bytes, result, store.interaction_closure_sha256


def publish_locked():
    # A second publisher cannot prepare against an in-flight transaction.
    updates, layer_bytes, receipt_bytes, result, closure_sha = prepare()
    for path, raw, new in updates:
        require(path.read_bytes() == raw, 'Concurrent navigation edit: ' + str(path))
    require(sha(closed.MANIFEST.read_bytes()) == closure_sha, 'Closure changed before publication')
    for rel, expected in result['source_sha256'].items():
        require(sha((ROOT / rel).read_bytes()) == expected, 'Startup source changed: ' + rel)
    LAYER_PATH.parent.mkdir(parents=True, exist_ok=True)
    written, created = [], []
    try:
        for path, raw, new in updates:
            require(path.read_bytes() == raw, 'Concurrent navigation edit: ' + str(path))
            atomic_replace(path, new)
            written.append((path, raw, new))
        # Exclusive creation prevents replacement of another publisher's layer/receipt.
        for path, content in ((LAYER_PATH, layer_bytes), (RECEIPT_PATH, receipt_bytes)):
            with path.open('xb') as stream:
                created.append((path, content))
                stream.write(content)
    except BaseException:
        for path, expected in reversed(created):
            # Own xb-created files can be partial after a failed write. The
            # transaction lock prevents a concurrent publisher from owning them.
            if path.exists() and expected.startswith(path.read_bytes()):
                path.unlink()
        for path, raw, new in reversed(written):
            if path.read_bytes() == new:
                atomic_replace(path, raw)
        raise
    print(json.dumps(result, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--check', action='store_true', help='Check readiness without writes (default)')
    modes.add_argument('--publish', action='store_true', help='Publish the four prefixes exactly once')
    args = parser.parse_args()
    if args.publish:
        with publication_lock():
            publish_locked()
    else:
        updates, layer_bytes, receipt_bytes, result, closure_sha = prepare()
        print(json.dumps(dict(result, status='READY', files_written=0), ensure_ascii=False))


if __name__ == '__main__':
    main()
