"""Read-only verification of round 1087 assets and independent review.

Numerics check the declared effective model. They do not certify autonomous
generation, an empirical error budget, or absence of a Lorentz completion.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess
import sys
from urllib.parse import unquote

BASE = Path(__file__).resolve().parent
STAGE = BASE.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--recompute', action='store_true')
args = parser.parse_args()
receipt = json.loads((BASE / 'acceptance.json').read_text('utf8'))
assert receipt['status'] == 'PASS'
assert receipt['science_increment'] == 0
assert receipt['goal_status'] == 'active'
assert receipt['FUCP_plus_CO_counterexample_qualified'] is False
assert receipt['current_goal_conclusion_proved'] is False
for name, digest in receipt['assets_sha256'].items():
    assert hashlib.sha256((STAGE / name).read_bytes()).hexdigest() == digest, name
review = json.loads((BASE / 'independent_review_checks.json').read_text('utf8'))
assert review['review_status'] == 'PASS'
if args.recompute:
    process = subprocess.run(
        [sys.executable, '-B', '-X', 'utf8', str(BASE / 'check.py'), '--check'],
        check=True, capture_output=True, text=True, encoding='utf8',
    )
    actual = json.loads(process.stdout)
    assert actual['passed']
    assert actual['exact']['histories'][0]['reciprocity_gap'] == '9/10'
    assert actual['exact']['histories'][1]['reciprocity_gap'] == '189/85'
    rank_process = subprocess.run(
        [sys.executable, '-B', '-X', 'utf8', str(BASE / 'fucp_interface_rank_check.py'), '--check'],
        check=True, capture_output=True, text=True, encoding='utf8',
    )
    assert json.loads(rank_process.stdout)['passed']
links = 0
for file in [STAGE / 'research_note_1087.md', *BASE.glob('*.md')]:
    text = file.read_text('utf8')
    for target in re.findall(r'\]\(([^)]+)\)', text):
        target = unquote(target.strip().strip('<>').split('#')[0])
        if not target or re.match(r'^(?:https?://|mailto:|codex:)', target):
            continue
        path = Path(target)
        if not path.is_absolute():
            path = file.parent / path
        assert path.exists(), (str(file), target)
        links += 1
    assert text.count('\n$$\n') % 2 == 0, str(file)
print(json.dumps({
    'round': 1087, 'passed': True, 'frozen_round_assets': len(receipt['assets_sha256']),
    'independent_review': 'PASS', 'local_links': links,
    'recomputed': args.recompute, 'scope': receipt['scope'],
    'generative_completeness_disproved': False,
    'FUCP_plus_CO_counterexample_qualified': False,
    'current_goal_conclusion_proved': False,
    'fundamental_lorentz_theory_excluded': False, 'goal_status': 'active',
}, ensure_ascii=False))
