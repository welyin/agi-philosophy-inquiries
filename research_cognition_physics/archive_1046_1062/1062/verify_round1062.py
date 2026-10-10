"""Read-only verification of the accepted, explicitly limited round 1062."""
from pathlib import Path
from hashlib import sha256
from urllib.parse import unquote
import json
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
SHARED = HERE.parent / '_shared'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    acceptance = json.loads((HERE / 'mainline_acceptance.json').read_text(encoding='utf-8'))
    assert acceptance['accepted'] and acceptance['round'] == 1062
    assert not acceptance['whole_roadmap_completed']
    assert not acceptance['full_cognitive_contract_counterexample']
    assert not acceptance['standard_thermal_history_refuted']
    assert acceptance['new_empirical_groups'] == acceptance['new_cognitive_axioms'] == 0
    for group in ('author_assets_sha256', 'review_assets_sha256', 'mainline_assets_sha256'):
        for name, expected in acceptance[group].items():
            assert digest(HERE / name) == expected, name
    for name, expected in acceptance['adoption_assets_sha256'].items():
        assert digest(SHARED / name) == expected, name
    for name, expected in acceptance['baseline_receipts_sha256'].items():
        assert digest(SHARED / name) == expected, name
    results = json.loads((HERE / 'results.json').read_text(encoding='utf-8'))
    counts = results['count']
    assert results['passed'] and results['round'] == 1062
    assert counts['candidate_scientific_calibration_groups'] == 1
    assert acceptance['cumulative_scientific_calibrations'] == counts['baseline_cumulative'] + 1 == 3838
    assert acceptance['stage_groups'] == counts['baseline_stage_groups'] + 1 == 17
    for name in ('mathematical_review.md', 'physical_review.md'):
        text = (HERE / name).read_text(encoding='utf-8')
        assert '通过' in text
        assert all(h in text for h in acceptance['author_assets_sha256'].values())
    local_links = 0
    paths = [HERE / name for group in ('author_assets_sha256', 'review_assets_sha256', 'mainline_assets_sha256')
             for name in acceptance[group] if name.endswith('.md')]
    paths += [SHARED / name for name in acceptance['adoption_assets_sha256'] if name.endswith('.md')]
    for path in paths:
        for ref in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)', path.read_text(encoding='utf-8')):
            ref = ref.strip().strip('<>')
            if re.match(r'[a-zA-Z][\w+.-]*://', ref) or ref.startswith('#'):
                continue
            assert (path.parent / unquote(ref.split('#', 1)[0])).is_file(), (path, ref)
            local_links += 1
    for program in (HERE / 'check.py', SHARED / 'verify_integration_1061.py'):
        run = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(program)],
                             capture_output=True, text=True, encoding='utf-8', check=True)
        assert run.stdout.strip()
        print(run.stdout.strip())
    print(json.dumps({'round': 1062, 'accepted_limited_contract': True,
                      'asset_links_checked': local_links,
                      'cumulative': 3838, 'stage_groups': 17,
                      'whole_roadmap_completed': False}))


if __name__ == '__main__':
    main()
