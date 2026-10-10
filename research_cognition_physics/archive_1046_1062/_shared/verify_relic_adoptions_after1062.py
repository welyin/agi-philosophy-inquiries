"""Verify two mature relic interfaces without creating a scientific round."""
from pathlib import Path
from hashlib import sha256
from urllib.parse import unquote
import argparse
import json
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PHASE = HERE.parent
RECEIPT = HERE / 'relic_adoptions_after1062_checks.json'
GROUPS = ('relic_neutrino_capture', 'scalar_relic_fraction')


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main(freeze=False):
    if freeze:
        assert not RECEIPT.exists(), 'Existing receipts must not be overwritten.'
    baseline = PHASE / '1062/mainline_acceptance.json'
    base = json.loads(baseline.read_text(encoding='utf-8'))
    assert base['accepted'] and base['round'] == 1062
    assert base['cumulative_scientific_calibrations'] == 3838
    assert base['stage_groups'] == 17 and not base['whole_roadmap_completed']
    owned = ['integration_relic_interfaces_after1062.md', Path(__file__).name]
    for group in GROUPS:
        assets = [group + '_adoption.md', group + '/sources.md']
        code = HERE / group / 'check.py'
        result = HERE / group / 'results.json'
        assert code.exists() == result.exists(), group
        if code.exists():
            assets += [group + '/check.py', group + '/results.json']
            run = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(code)],
                                 check=True, capture_output=True, text=True, encoding='utf-8')
            assert run.stdout.strip(), group
            print(run.stdout.strip())
        review_name = group + '_adoption_review.md'
        review = (HERE / review_name).read_text(encoding='utf-8')
        assert '通过' in review and all(digest(HERE / name) in review for name in assets)
        owned.extend(assets + [review_name])
    local_links = 0
    for name in owned:
        path = HERE / name
        if path.suffix != '.md':
            continue
        for ref in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)', path.read_text(encoding='utf-8')):
            ref = ref.strip().strip('<>')
            if re.match(r'[a-zA-Z][\w+.-]*://', ref) or ref.startswith('#'):
                continue
            dest = (path.parent / unquote(ref.split('#', 1)[0])).resolve()
            assert dest.is_file() or (freeze and dest == RECEIPT.resolve()), (name, ref)
            local_links += 1
    stable = {
        'accepted': True, 'new_scientific_groups': 0, 'new_empirical_groups': 0,
        'new_cognitive_axioms': 0, 'latest_accepted_round': 1062,
        'cumulative_scientific_calibrations': 3838, 'stage_groups': 17,
        'whole_roadmap_completed': False,
        'assets_sha256': {name: digest(HERE / name) for name in owned},
        'history_sha256': {
            '../1062/mainline_acceptance.json': digest(baseline),
            'joint_scope_after1062.md': digest(HERE / 'joint_scope_after1062.md'),
            'nucleon_source_adoption.md': digest(HERE / 'nucleon_source_adoption.md'),
            'neutrino_decoupling_adoption.md': digest(HERE / 'neutrino_decoupling_adoption.md'),
        },
    }
    run = subprocess.run([sys.executable, '-B', '-X', 'utf8',
                          str(PHASE / '1062/verify_round1062.py')],
                         check=True, capture_output=True, text=True, encoding='utf-8')
    print(run.stdout.strip())
    if freeze:
        with RECEIPT.open('x', encoding='utf-8', newline='\n') as out:
            json.dump(stable, out, ensure_ascii=False, indent=2)
            out.write('\n')
    else:
        assert json.loads(RECEIPT.read_text(encoding='utf-8')) == stable
    print(json.dumps({'relic_interfaces_accepted': True,
                      'assets': len(owned), 'local_links': local_links,
                      'new_scientific_groups': 0, 'whole_roadmap_completed': False}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze-exclusive', action='store_true')
    main(parser.parse_args().freeze_exclusive)
