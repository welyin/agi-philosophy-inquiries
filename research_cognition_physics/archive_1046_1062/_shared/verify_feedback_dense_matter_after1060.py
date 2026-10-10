"""Freeze or read-only verify two mature physical interfaces after round 1060."""
from pathlib import Path
import argparse
import hashlib
import json
import runpy
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PHASE = HERE.parent
ROOT = PHASE.parent
RECEIPT = HERE / 'feedback_dense_matter_after1060_checks.json'
BASE = runpy.run_path(str(PHASE / '1060/verify_round1060.py'))
NAV = BASE['NAV']
GROUPS = {
    'autonomous_feedback': [
        'autonomous_feedback_adoption.md',
        'autonomous_feedback/sources.md',
        'autonomous_feedback/check.py',
        'autonomous_feedback/results.json',
        '../_admission/after1060_autonomous_feedback/selection.md',
    ],
    'neutron_star_tidal': [
        'neutron_star_tidal_adoption.md',
        'neutron_star_tidal/sources.md',
        '../_admission/after1060_dense_matter/selection.md',
    ],
}
INTEGRATION = 'integration_feedback_dense_matter_after1060.md'
HISTORY = [
    'archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md',
    'archive_956_989/981/drafts/common_parent_contract_v1.md',
    'archive_956_989/research_note_975.md',
    'archive_956_989/research_note_980.md',
    'archive_990_1008/research_note_1002.md',
    'archive_1046_/_shared/normal_contact_transport_adoption.md',
    'archive_1046_/_shared/nucleon_source_adoption.md',
    'archive_1046_/1060/mainline_acceptance.json',
    'archive_1046_/1060/research_round_1060_checks.json',
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def links(path, freeze):
    # Reuse the already reviewed local-link rule, allowing only this pending receipt.
    if not freeze or RECEIPT.exists():
        return BASE['links'](path)
    import re
    from urllib.parse import unquote
    text = path.read_text(encoding='utf-8')
    tails = BASE['BASE']['BASE']['TAILS']
    if path.parent == ROOT and path.name in tails:
        text = text.split(tails[path.name][0], 1)[0]
    text = re.sub(r'\\\[[\s\S]*?\\\]|\$\$[\s\S]*?\$\$', '', text)
    count = 0
    for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)', text):
        target = target.strip().strip('<>')
        if re.match(r'[a-zA-Z][\w+.-]*://', target) or target.startswith('#'):
            continue
        dest = (path.parent / unquote(target.split('#', 1)[0])).resolve()
        assert dest.is_file() or dest == RECEIPT.resolve(), (path, target)
        count += 1
    return count


def run(freeze):
    old = load(PHASE / '1060/research_round_1060_checks.json')
    assert old['all_checks_passed'] and old['round'] == 1060
    for field in ('author_assets_sha256', 'review_assets_sha256'):
        for rel, expected in old[field].items():
            assert sha(PHASE / '1060' / rel) == expected, rel
    for rel, expected in old['historical_sha256'].items():
        assert sha(ROOT / rel) == expected, rel
    assert sha(PHASE / '1060/verify_round1060.py') == old['verifier_sha256']
    for name, (marker, expected) in BASE['BASE']['BASE']['TAILS'].items():
        raw = (ROOT / name).read_bytes()
        assert hashlib.sha256(raw[raw.index(marker.encode('utf-8')):]).hexdigest() == expected
    raw = (ROOT / 'ROADMAP.md').read_bytes()
    marker = '## M1 检验认知要求的独立选择力'.encode('utf-8')
    assert hashlib.sha256(raw[raw.index(marker):]).hexdigest() == old['original_roadmap_requirements_sha256']

    groups = {k: list(v) for k, v in GROUPS.items()}
    # The dense-matter author may include a finite algebra diagnostic; never a fitted EOS.
    for name in ('check.py', 'results.json'):
        path = HERE / 'neutron_star_tidal' / name
        if path.exists():
            groups['neutron_star_tidal'].append('neutron_star_tidal/' + name)
    assert (HERE / 'neutron_star_tidal/check.py').exists() == (HERE / 'neutron_star_tidal/results.json').exists()
    own = [INTEGRATION, Path(__file__).name]
    reviews = []
    for key, assets in groups.items():
        review_name = key + '_adoption_review.md'
        review = (HERE / review_name).read_text(encoding='utf-8')
        assert '通过' in review, review_name
        for name in assets:
            assert sha(HERE / name) in review, name
        own.extend(assets + [review_name])
        reviews.append(review)
        code = HERE / key / 'check.py'
        if code.exists():
            cp = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(code)],
                                capture_output=True, text=True, encoding='utf-8', check=True)
            if key == 'neutron_star_tidal':
                assert cp.stdout.startswith('PASS: exact units, tidal weights and phase conventions;'), cp.stdout
                assert load(HERE / key / 'results.json')['all_assertions_passed']
            else:
                assert cp.stdout.startswith('Read-only four-state, energy, information and coarse-channel checks passed.'), cp.stdout
                assert load(HERE / key / 'results.json')['all_checks_passed']
    assert any(sha(HERE / INTEGRATION) in review for review in reviews)
    own_links = sum(links(HERE / name, freeze) for name in own if name.endswith('.md'))
    nav_links = sum(links(path, freeze) for path in NAV)
    stable = {
        'date': '2026-10-09', 'kind': 'feedback_and_dense_matter_mature_adoptions',
        'latest_accepted_round_at_batch': 1060,
        'cumulative_scientific_calibrations': 3836, 'stage_groups': 15,
        'new_scientific_groups': 0, 'new_empirical_groups': 0, 'new_cognitive_axioms': 0,
        'assets_sha256': {name: sha(HERE / name) for name in own},
        'historical_sha256': {name: sha(ROOT / name) for name in HISTORY},
        'original_roadmap_requirements_sha256': old['original_roadmap_requirements_sha256'],
        'asset_local_links_checked': own_links,
        'independent_physics_and_scope_reviews_passed': True,
        'complete_finite_reservoir_history_proved': False,
        'high_density_EOS_derived_from_cognition_or_SM': False,
        'two_interfaces_combined_into_one_history': False,
        'whole_roadmap_completed': False, 'goal_objective_or_status_mutated': False,
        'all_checks_passed': True,
    }
    if freeze:
        record = dict(stable)
        record['navigation_snapshot'] = {
            'sha256': {str(path.relative_to(ROOT)): sha(path) for path in NAV},
            'local_links_checked': nav_links,
        }
        with RECEIPT.open('x', encoding='utf-8') as stream:
            json.dump(record, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        record = load(RECEIPT)
        record.pop('navigation_snapshot')
        assert record == stable
    print(json.dumps({'all_checks_passed': True, 'assets': len(own),
                      'history': len(HISTORY), 'navigation_links': nav_links,
                      'new_scientific_groups': 0, 'whole_roadmap_completed': False}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    run(parser.parse_args().freeze)
