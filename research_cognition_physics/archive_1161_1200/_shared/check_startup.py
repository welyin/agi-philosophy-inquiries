"""Audit stage launch only; this is not a scientific research round."""
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
RESEARCH = STAGE.parent
ROOT = RESEARCH.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    receipt = json.loads((HERE / 'goal_start_receipt.json').read_text(encoding='utf-8-sig'))
    ledger = json.loads((HERE / 'dependency_ledger.json').read_text(encoding='utf-8-sig'))
    goal = (HERE / 'goal.md').read_text(encoding='utf-8-sig')
    review = (STAGE / '_admission/scope_review.md').read_text(encoding='utf-8-sig')
    assert receipt['goal']['status'] == ledger['status'] == 'active'
    assert receipt['goal']['objective'] in goal
    assert receipt['first_round'] == ledger['first_round'] == 1161
    assert ledger['accepted_stage_rounds'] == []
    assert not list(STAGE.glob('research_note_*.md'))
    assert not (RESEARCH / 'archive_1159_').exists()
    assert len(list(p for p in (RESEARCH / 'archive_1159_1160').rglob('*') if p.is_file())) == 41
    bound = [HERE / 'goal.md', STAGE / '_admission/reuse_and_first_question.md']
    for path in bound:
        assert sha(path) in review, path
    migration = RESEARCH / '_migration/closure_1159_1160/verification.json'
    assert migration.is_file()
    json.loads(migration.read_text(encoding='utf-8-sig'))
    links = []
    output = HERE / 'startup_checks.json'
    for path in STAGE.rglob('*.md'):
        base = RESEARCH if path.name == 'navigation_notice.md' else path.parent
        body = re.sub(r'```.*?```', '', path.read_text(encoding='utf-8-sig'), flags=re.S)
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)', body):
            target = unquote(target.strip().strip('<>').split('#', 1)[0])
            if not target or re.match(r'^[a-zA-Z]+://', target):
                continue
            resolved = (base / target).resolve()
            assert resolved.exists() or resolved == output, (path, target)
            links.append({'source': path.relative_to(RESEARCH).as_posix(), 'target': target})
    result = {
        'schema': 'cognitive_gr_stage_1161_startup_checks_v1',
        'status': 'PASS_STARTUP_ONLY',
        'date': '2026-10-12',
        'goal_status_at_launch': 'active',
        'objective_matches_saved_goal': True,
        'scope_review_hash_bindings_match': True,
        'moved_files': 41,
        'old_directory_absent': True,
        'local_links_checked': len(links),
        'migration_verification': migration.relative_to(ROOT).as_posix(),
        'migration_verification_sha256': sha(migration),
        'goal_receipt_sha256': sha(HERE / 'goal_start_receipt.json'),
        'scope_review_sha256': sha(STAGE / '_admission/scope_review.md'),
        'scientific_results_added': 0,
        'first_research_round': 1161,
        'note': 'Administrative launch audit only. Goal status is bound to the actual launch receipt, not a permanent claim about future status.',
        'links': links,
    }
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k != 'links'}, ensure_ascii=False))


if __name__ == '__main__':
    main()
