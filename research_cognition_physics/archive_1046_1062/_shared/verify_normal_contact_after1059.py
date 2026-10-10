"""Read-only verification of the count-zero normal-contact adoption batch.

--freeze exclusively creates the first receipt. Current navigation may evolve;
the original roadmap requirements and frozen evidence may not.
"""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import re
import runpy
import subprocess
import sys
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
PHASE = HERE.parent
ROOT = PHASE.parent
RECEIPT = HERE / 'normal_contact_after1059_checks.json'
BASE = runpy.run_path(str(HERE / 'verify_adoptions_after1058.py'))
NAV = BASE['BASE']['NAV']
OWN = [
    'normal_contact_transport_adoption.md',
    'normal_contact_transport/sources.md',
    'normal_contact_transport/check.py',
    'normal_contact_transport/results.json',
    'normal_contact_transport_adoption_review.md',
    'normal_contact_transport/review_checks.json',
    '../_admission/after1059_macroscopic_bridge/selection.md',
    '../_admission/after1059_macroscopic_bridge/scope_review.md',
    '../_admission/after1059_principle_selection/selection.md',
    '../_admission/after1059_dynamic_source/selection.md',
    'integration_normal_contact_after1059.md',
    'integration_normal_contact_after1059_review.md',
    'verify_normal_contact_after1059.py',
]
HISTORY = [
    'archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md',
    'archive_956_989/981/drafts/common_parent_contract_v1.md',
    'archive_990_1008/research_note_1000.md',
    'archive_990_1008/research_note_1006.md',
    'archive_990_1008/1008/overall_operation_hypothesis_v2_2.md',
    'archive_1009_1043/1043/completion_audit.md',
    'archive_1046_/research_note_1047.md',
    'archive_1046_/research_note_1048.md',
    'archive_1046_/research_note_1059.md',
    'archive_1046_/1059/mainline_acceptance.json',
    'archive_1046_/_shared/integration_1059_checks.json',
    'archive_1046_/_admission/after1059_selection_checks.json',
    'archive_1046_/_shared/goal_start.json',
]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def local_links(path, allow_receipt=False):
    body = path.read_text(encoding='utf-8')
    if path.parent == ROOT and path.name in BASE['BASE']['TAILS']:
        body = body.split(BASE['BASE']['TAILS'][path.name][0], 1)[0]
    body = re.sub(r'\\\[[\s\S]*?\\\]', '', body)
    body = re.sub(r'\$\$[\s\S]*?\$\$', '', body)
    count = 0
    for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)', body):
        target = target.strip().strip('<>')
        if re.match(r'[a-zA-Z][\w+.-]*://', target) or target.startswith('#'):
            continue
        dest = (path.parent / unquote(target.split('#', 1)[0])).resolve()
        assert dest.is_file() or (allow_receipt and dest == RECEIPT.resolve()), (path, target)
        count += 1
    return count


def run(freeze=False):
    author = HERE / 'normal_contact_transport_adoption.md'
    checks = load(HERE / 'normal_contact_transport/review_checks.json')
    assert checks['review_status'] == 'passed' and checks['science_count'] == 0
    assert checks['numerical_scalar_comparisons'] == 32
    assert checks['max_absolute_difference'] < checks['comparison_tolerance']
    for rel, expected in checks['author_sha256'].items():
        assert digest(HERE / rel) == expected, rel
    review = (HERE / 'normal_contact_transport_adoption_review.md').read_text(encoding='utf-8')
    scope = (HERE / '../_admission/after1059_macroscopic_bridge/scope_review.md').read_text(encoding='utf-8')
    assert '通过' in review and '通过' in scope
    for expected in checks['author_sha256'].values():
        assert expected in review and expected in scope
    assert digest(HERE / '../_admission/after1059_macroscopic_bridge/selection.md') in scope
    assert digest(HERE / 'normal_contact_transport/review_checks.json') in review
    integration_review = (HERE / 'integration_normal_contact_after1059_review.md').read_text(encoding='utf-8')
    assert '通过' in integration_review
    assert digest(HERE / 'integration_normal_contact_after1059.md') in integration_review

    output = subprocess.run(
        [sys.executable, '-B', '-X', 'utf8', str(HERE / 'normal_contact_transport/check.py')],
        check=True, capture_output=True, text=True, encoding='utf-8')
    diagnostic = json.loads(output.stdout)
    assert diagnostic['passed'] and diagnostic['new_science'] == 0

    # Preserve accepted evidence without repeating unchanged empirical calculations.
    baseline = load(HERE / 'integration_1059_checks.json')
    assert baseline['all_checks_passed'] and baseline['latest_accepted_round_at_batch'] == 1059
    for rel, expected in baseline['assets_sha256'].items():
        assert digest(PHASE / rel) == expected, rel
    for rel, expected in baseline['historical_sha256'].items():
        assert digest(ROOT / rel) == expected, rel
    for name, (marker, expected) in BASE['BASE']['TAILS'].items():
        raw = (ROOT / name).read_bytes()
        assert sha256(raw[raw.index(marker.encode('utf-8')):]).hexdigest() == expected
    raw = (ROOT / 'ROADMAP.md').read_bytes()
    marker = '## M1 检验认知要求的独立选择力'.encode('utf-8')
    assert sha256(raw[raw.index(marker):]).hexdigest() == BASE['ROADMAP_SHA']

    assets = {rel: digest(HERE / rel) for rel in OWN}
    links = sum(local_links(HERE / rel, freeze) for rel in OWN if rel.endswith('.md'))
    navlinks = sum(local_links(path, freeze) for path in NAV)
    stable = {
        'date': '2026-10-09', 'kind': 'normal_contact_mature_adoption_and_scope_audits',
        'latest_accepted_round_at_batch': 1059,
        'cumulative_scientific_calibrations': 3835, 'stage_groups': 14,
        'new_scientific_groups': 0, 'new_empirical_groups': 0, 'new_cognitive_axioms': 0,
        'assets_sha256': assets,
        'historical_sha256': {rel: digest(ROOT / rel) for rel in HISTORY},
        'asset_local_links_checked': links,
        'original_roadmap_requirements_sha256': BASE['ROADMAP_SHA'],
        'physical_and_scope_reviews_passed': True,
        'integration_scope_review_passed': True,
        'author_readonly_reproduction_passed': True,
        'accepted_1059_frozen_evidence_unchanged': True,
        'new_joint_electrical_thermal_noise_experiment': False,
        'full_finite_reservoir_or_apparatus_history_certified': False,
        'goal_tool_status_observed': 'active', 'goal_objective_or_status_mutated': False,
        'whole_roadmap_complete': False, 'all_checks_passed': True,
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
        assert record == stable, 'Frozen batch changed; preserve it and create a new batch.'
    # The only temporarily pending target above must exist after the first freeze.
    assert RECEIPT.is_file()
    return stable, navlinks


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    args = parser.parse_args()
    record, navlinks = run(args.freeze)
    print(json.dumps({'all_checks_passed': True,
                      'assets': len(record['assets_sha256']),
                      'history': len(record['historical_sha256']),
                      'current_navigation_links': navlinks,
                      'latest_round_at_batch': 1059, 'new_science_groups': 0,
                      'whole_roadmap_complete': False}))
