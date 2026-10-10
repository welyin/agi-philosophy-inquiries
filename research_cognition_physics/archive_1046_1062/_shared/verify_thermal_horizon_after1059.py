"""Verify mature thermal-history and horizon interfaces without rewriting them."""
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
RECEIPT = HERE / 'thermal_horizon_after1059_checks.json'
BASE = runpy.run_path(str(HERE / 'verify_adoptions_after1058.py'))
NAV = BASE['BASE']['NAV']
OWN = [
    'electroweak_thermal_history_adoption.md',
    'electroweak_thermal_history/sources.md',
    'electroweak_thermal_history/check.py',
    'electroweak_thermal_history/results.json',
    'electroweak_thermal_history/review_check.py',
    'electroweak_thermal_history/review_checks.json',
    'electroweak_thermal_history_adoption_review.md',
    'maxwell_horizon_scattering_adoption.md',
    'maxwell_horizon_scattering_adoption_review.md',
    '../_admission/after1059_asymmetry_history/selection.md',
    '../_admission/after1059_horizon_scattering/selection.md',
    '../_admission/after1059_joint_parameters/selection.md',
    'integration_thermal_horizon_after1059.md',
    'integration_thermal_horizon_after1059_review.md',
    'verify_thermal_horizon_after1059.py',
]
HISTORY = [
    'archive_702_741/research_note_709.md',
    'archive_956_989/981/drafts/common_parent_contract_v1.md',
    'archive_990_1008/research_note_994.md',
    'archive_990_1008/research_note_995.md',
    'archive_990_1008/research_note_996.md',
    'archive_990_1008/research_note_1004.md',
    'archive_990_1008/1008/overall_operation_hypothesis_v2_2.md',
    'archive_1046_/1059/mainline_acceptance.json',
    'archive_1046_/_shared/normal_contact_after1059_checks.json',
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
    old = load(HERE / 'normal_contact_after1059_checks.json')
    assert old['all_checks_passed'] and old['latest_accepted_round_at_batch'] == 1059
    for rel, expected in old['assets_sha256'].items():
        assert digest(HERE / rel) == expected, rel
    for rel, expected in old['historical_sha256'].items():
        assert digest(ROOT / rel) == expected, rel

    ew_review = (HERE / 'electroweak_thermal_history_adoption_review.md').read_text(encoding='utf-8')
    independent = load(HERE / 'electroweak_thermal_history/review_checks.json')
    assert independent['all_checks_passed'] and independent['new_science_groups'] == 0
    for rel, expected in independent['author_sha256'].items():
        assert digest(HERE / rel) == expected, rel
        assert expected in ew_review, rel
    assert '通过' in ew_review
    for name in ('maxwell_horizon_scattering_adoption', 'integration_thermal_horizon_after1059'):
        review = (HERE / f'{name}_review.md').read_text(encoding='utf-8')
        assert '通过' in review and digest(HERE / f'{name}.md') in review, name
    for name, field in [('check.py', 'passed'), ('review_check.py', 'all_checks_passed')]:
        output = subprocess.run(
            [sys.executable, '-B', '-X', 'utf8', str(HERE / 'electroweak_thermal_history' / name)],
            check=True, capture_output=True, text=True, encoding='utf-8')
        assert json.loads(output.stdout)[field]

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
        'date': '2026-10-09', 'kind': 'mature_electroweak_thermal_and_maxwell_horizon_interfaces',
        'previous_goal_turn_classification': 'progress',
        'latest_accepted_round_at_batch': 1059,
        'cumulative_scientific_calibrations': 3835, 'stage_groups': 14,
        'new_scientific_groups': 0, 'new_empirical_groups': 0, 'new_cognitive_axioms': 0,
        'assets_sha256': assets,
        'historical_sha256': {rel: digest(ROOT / rel) for rel in HISTORY},
        'asset_local_links_checked': links,
        'original_roadmap_requirements_sha256': BASE['ROADMAP_SHA'],
        'two_adoption_reviews_passed': True,
        'integration_review_passed': True,
        'central_equation_independent_reproduction_passed': True,
        'lattice_or_greybody_calculation_repeated': False,
        'real_baryon_yield_generated': False,
        'full_semiclassical_backreaction_certified': False,
        'new_observational_conflict_from_lithium': False,
        'goal_objective_or_status_mutated': False,
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
        assert record == stable, 'Preserve frozen evidence and create a new batch for changes.'
    assert RECEIPT.is_file()
    return stable, navlinks


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    args = parser.parse_args()
    record, navlinks = run(args.freeze)
    print(json.dumps({'all_checks_passed': True, 'assets': len(record['assets_sha256']),
                      'history': len(record['historical_sha256']),
                      'current_navigation_links': navlinks, 'new_science_groups': 0,
                      'whole_roadmap_complete': False}))
