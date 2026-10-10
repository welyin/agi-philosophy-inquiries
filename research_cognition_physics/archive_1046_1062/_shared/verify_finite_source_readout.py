"""Check the count-zero finite source-readout adoption; --freeze is exclusive."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import runpy
import subprocess
import sys
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RECEIPT = HERE/'finite_source_readout_checks.json'
PREVIOUS = runpy.run_path(str(HERE/'verify_thermal_horizon_after1059.py'))
BASE = PREVIOUS['BASE']
NAV = PREVIOUS['NAV']
OWN = ['finite_source_readout_adoption.md', 'finite_source_readout/check.py',
       'finite_source_readout/results.json', 'finite_source_readout_review.md',
       'verify_finite_source_readout.py']
HISTORY = ['archive_078_099/research_note_93.md',
           'archive_956_989/research_note_962.md',
           'archive_956_989/research_note_974.md',
           'archive_956_989/research_note_977.md',
           'archive_1009_1043/research_note_1033.md',
           'archive_1046_/research_note_1053.md',
           'archive_1046_/1053/proof.md',
           'archive_1046_/_shared/thermal_horizon_after1059_checks.json']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def links(p, freeze):
    text = p.read_text(encoding='utf-8')
    if p.parent == ROOT and p.name in BASE['BASE']['TAILS']:
        text = text.split(BASE['BASE']['TAILS'][p.name][0], 1)[0]
    text = re.sub(r'\\\[[\s\S]*?\\\]|\$\$[\s\S]*?\$\$', '', text)
    n = 0
    for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)', text):
        target = target.strip().strip('<>')
        if re.match(r'[a-zA-Z][\w+.-]*://', target) or target.startswith('#'):
            continue
        dest = (p.parent/unquote(target.split('#', 1)[0])).resolve()
        assert dest.is_file() or (freeze and dest == RECEIPT.resolve()), (p, target)
        n += 1
    return n


def run(freeze):
    # The predecessor checks immutable scientific assets and original roadmap tails.
    # During first creation its navigation scanner needs the new receipt, so perform
    # the same immutable predecessor checks directly and check current navigation here.
    previous = json.loads((HERE/'thermal_horizon_after1059_checks.json').read_text(encoding='utf-8'))
    assert previous['all_checks_passed']
    for rel, expected in previous['assets_sha256'].items():
        assert sha(HERE/rel) == expected, rel
    for rel, expected in previous['historical_sha256'].items():
        assert sha(ROOT/rel) == expected, rel
    for name, (marker, expected) in BASE['BASE']['TAILS'].items():
        raw = (ROOT/name).read_bytes()
        assert hashlib.sha256(raw[raw.index(marker.encode('utf-8')):]).hexdigest() == expected
    raw = (ROOT/'ROADMAP.md').read_bytes()
    marker = '## M1 检验认知要求的独立选择力'.encode('utf-8')
    assert hashlib.sha256(raw[raw.index(marker):]).hexdigest() == BASE['ROADMAP_SHA']

    review = (HERE/'finite_source_readout_review.md').read_text(encoding='utf-8')
    assert '通过' in review
    for name in OWN[:3]:
        assert sha(HERE/name) in review, name
    check = subprocess.run([sys.executable, '-B', '-X', 'utf8',
                            str(HERE/'finite_source_readout/check.py')],
                           capture_output=True, text=True, encoding='utf-8', check=True)
    assert json.loads(check.stdout)['passed']
    own_links = sum(links(HERE/name, freeze) for name in OWN if name.endswith('.md'))
    nav_links = sum(links(p, freeze) for p in NAV)
    stable = {'date': '2026-10-09', 'kind': 'finite_source_readout_mature_adoption',
              'latest_accepted_round_at_batch': 1059,
              'cumulative_scientific_calibrations': 3835, 'stage_groups': 14,
              'new_scientific_groups': 0, 'new_empirical_groups': 0,
              'new_cognitive_axioms': 0,
              'assets_sha256': {name: sha(HERE/name) for name in OWN},
              'historical_sha256': {name: sha(ROOT/name) for name in HISTORY},
              'original_roadmap_requirements_sha256': BASE['ROADMAP_SHA'],
              'asset_local_links_checked': own_links,
              'mathematical_and_scope_review_passed': True,
              'actual_gravitational_instrument_certified': False,
              'new_cognitive_selection_proved': False,
              'goal_objective_or_status_mutated': False,
              'whole_roadmap_complete': False, 'all_checks_passed': True}
    if freeze:
        record = dict(stable)
        record['navigation_snapshot'] = {'sha256': {str(p.relative_to(ROOT)): sha(p) for p in NAV},
                                         'local_links_checked': nav_links}
        with RECEIPT.open('x', encoding='utf-8') as f:
            json.dump(record, f, ensure_ascii=False, indent=2)
            f.write('\n')
    else:
        record = json.loads(RECEIPT.read_text(encoding='utf-8'))
        record.pop('navigation_snapshot')
        assert record == stable
    print(json.dumps({'all_checks_passed': True, 'assets': len(OWN),
                      'history': len(HISTORY), 'current_navigation_links': nav_links,
                      'new_science_groups': 0, 'whole_roadmap_complete': False}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    run(parser.parse_args().freeze)
