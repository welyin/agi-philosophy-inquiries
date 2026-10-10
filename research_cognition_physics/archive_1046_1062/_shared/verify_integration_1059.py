"""Verify round 1059 acceptance and navigation without rewriting evidence."""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import runpy
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PHASE = HERE.parent
ROOT = PHASE.parent
ROUND = PHASE / '1059'
RECEIPT = HERE / 'integration_1059_checks.json'
BASE = runpy.run_path(str(HERE / 'verify_adoptions_after1058.py'))
NAV = BASE['BASE']['NAV']
EXTRA = [
    'integration_1059.md', 'integration_1059_review.md',
    'verify_integration_1059.py',
    '../_admission/after1058_qed_prediction/selection.md',
    '../_admission/after1058_qed_prediction/statistical_pre_review.md',
    '../_admission/after1058_qed_prediction/physical_pre_review.md',
    '../_admission/after1058_rescaling_selection/selection.md',
    '../_admission/after1058_dispersion_selection/selection.md',
]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def run(freeze=False):
    author_path = ROUND / 'research_round_1059_checks.json'
    author = load(author_path)
    acceptance_path = ROUND / 'mainline_acceptance.json'
    acceptance = load(acceptance_path)
    assert acceptance['accepted'] and acceptance['round'] == 1059
    assert acceptance['author_receipt_sha256'] == digest(author_path)
    assert acceptance['cumulative_scientific_calibrations'] == 3835
    assert acceptance['stage_groups'] == 14
    assert acceptance['new_project_empirical_groups'] == 1
    assert acceptance['new_cognitive_axioms'] == 0
    assert not acceptance['full_GR_or_roadmap_completed']
    for rel, expected in author['owned_sha256'].items():
        assert digest(ROUND / rel) == expected, rel
    for rel, expected in author['historical_sha256'].items():
        assert digest(ROOT / rel) == expected, rel
    for rel, expected in acceptance['review_assets_sha256'].items():
        assert digest(ROUND / rel) == expected, rel
    for label in ('physical', 'statistical'):
        checks = load(ROUND / f'{label}_review_checks.json')
        field = 'all_checks_passed' if label == 'physical' else 'all_computational_checks_passed'
        assert checks[field]
        text = (ROUND / f'{label}_review.md').read_text(encoding='utf-8')
        assert digest(author_path) in text and '通过' in text
    integration_review = (HERE / 'integration_1059_review.md').read_text(encoding='utf-8')
    assert digest(HERE / 'integration_1059.md') in integration_review
    assert '通过' in integration_review

    scripts = ['verify_round1059.py', 'physical_review_check.py', 'statistical_review_check.py']
    for name in scripts:
        result = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(ROUND / name)],
                                check=True, capture_output=True, text=True, encoding='utf-8')
        output = json.loads(result.stdout)
        field = {'verify_round1059.py': 'author_checks_passed',
                 'physical_review_check.py': 'all_checks_passed',
                 'statistical_review_check.py': 'pass'}[name]
        assert output[field]

    for name, (marker, expected) in BASE['BASE']['TAILS'].items():
        raw = (ROOT / name).read_bytes()
        assert sha256(raw[raw.index(marker.encode('utf-8')):]).hexdigest() == expected
    raw = (ROOT / 'ROADMAP.md').read_bytes()
    marker = '## M1 检验认知要求的独立选择力'.encode('utf-8')
    assert sha256(raw[raw.index(marker):]).hexdigest() == BASE['ROADMAP_SHA']

    assets = {str((ROUND / rel).relative_to(HERE.parent)): expected
              for rel, expected in author['owned_sha256'].items() if not rel.startswith('..')}
    assets['research_note_1059.md'] = author['owned_sha256']['../research_note_1059.md']
    for path in [author_path, acceptance_path] + [ROUND / rel for rel in acceptance['review_assets_sha256']]:
        assets[str(path.relative_to(PHASE))] = digest(path)
    for rel in EXTRA:
        path = (HERE / rel).resolve()
        assets[str(path.relative_to(PHASE.resolve()))] = digest(path)
    links = sum(BASE['local_links'](PHASE / rel) for rel in assets if rel.endswith('.md'))
    navlinks = sum(BASE['local_links'](p) for p in NAV)
    stable = {
        'round': 1059, 'date': '2026-10-09', 'all_checks_passed': True,
        'assets_sha256': assets, 'historical_sha256': author['historical_sha256'],
        'prior_batch_sha256': digest(HERE / 'adoptions_after1058_checks.json'),
        'asset_local_links_checked': links, 'three_readonly_reproductions_passed': True,
        'two_independent_scientific_reviews_passed': True,
        'integration_scope_review_passed': True,
        'original_roadmap_requirements_sha256': BASE['ROADMAP_SHA'],
        'latest_accepted_round_at_batch': 1059,
        'cumulative_scientific_calibrations': 3835, 'stage_groups': 14,
        'new_project_empirical_groups': 1, 'new_cognitive_axioms': 0,
        'goal_tool_status_observed': 'active', 'goal_objective_or_status_mutated': False,
        'whole_roadmap_complete': False,
    }
    if freeze:
        record = dict(stable)
        record['navigation_snapshot'] = {
            'sha256': {str(p.relative_to(ROOT)): digest(p) for p in NAV},
            'local_links_checked': navlinks,
        }
        with RECEIPT.open('x', encoding='utf-8') as stream:
            json.dump(record, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        record = load(RECEIPT)
        record.pop('navigation_snapshot')
        assert record == stable, 'Frozen integration changed; preserve it and add a new batch.'
    return stable, navlinks


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    args = parser.parse_args()
    record, navlinks = run(args.freeze)
    print(json.dumps({'all_checks_passed': True, 'round': 1059,
                      'assets': len(record['assets_sha256']),
                      'history': len(record['historical_sha256']),
                      'current_navigation_links': navlinks,
                      'new_project_empirical_groups': 1,
                      'whole_roadmap_complete': False}))
