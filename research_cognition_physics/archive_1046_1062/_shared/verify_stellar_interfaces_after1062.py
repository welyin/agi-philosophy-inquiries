"""Read-only verification of two mature stellar interfaces and a source audit."""
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
ROOT = PHASE.parent
RECEIPT = HERE / 'stellar_interfaces_after1062_checks.json'
AUDIT = PHASE / '_admission/dynamical_tidal_source_after1062'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main(freeze=False):
    if freeze:
        assert not RECEIPT.exists(), 'A frozen receipt cannot be overwritten.'
    base_path = PHASE / '1062/mainline_acceptance.json'
    base = json.loads(base_path.read_text(encoding='utf-8'))
    assert base['accepted'] and base['round'] == 1062
    assert base['cumulative_scientific_calibrations'] == 3838
    assert base['stage_groups'] == 17 and not base['whole_roadmap_completed']

    owned = [HERE / 'integration_stellar_interfaces_after1062.md', Path(__file__)]
    for group in ('neutron_star_dynamics', 'solar_energy_balance'):
        assets = [HERE / (group + '_adoption.md')]
        assets += [HERE / group / name for name in ('sources.md', 'check.py', 'results.json')]
        review = HERE / (group + '_adoption_review.md')
        review_text = review.read_text(encoding='utf-8')
        assert '通过' in review_text
        assert all(digest(asset) in review_text for asset in assets)
        owned.extend(assets + [review])
        run = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(HERE / group / 'check.py')],
                             capture_output=True, text=True, encoding='utf-8', check=True)
        print(run.stdout.strip())

    audit_assets = [AUDIT / name for name in ('source_audit.md', 'check.py', 'results.json')]
    audit_review = AUDIT / 'integration_review.md'
    audit_text = audit_review.read_text(encoding='utf-8')
    assert '通过' in audit_text and all(digest(asset) in audit_text for asset in audit_assets)
    owned.extend(audit_assets + [audit_review])
    run = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(AUDIT / 'check.py')],
                         capture_output=True, text=True, encoding='utf-8', check=True)
    print(run.stdout.strip())

    links = 0
    for path in owned:
        if path.suffix != '.md':
            continue
        for ref in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)', path.read_text(encoding='utf-8')):
            ref = ref.strip().strip('<>')
            if re.match(r'[a-zA-Z][\w+.-]*://', ref) or ref.startswith('#'):
                continue
            target = (path.parent / unquote(ref.split('#', 1)[0])).resolve()
            assert target.is_file() or (freeze and target == RECEIPT.resolve()), (path, ref)
            links += 1

    markers = {
        'README.md': '## 当前入口：1044—1045阶段已结项',
        'research_direction.md': '以下保存上一阶段结项及历史方向，不覆盖当前目标。',
        'RESEARCH_STATE.md': '以下为1044—1045结项快照和更早历史，不覆盖本栏。',
        'ROADMAP.md': '## M1 检验认知要求的独立选择力',
    }
    expected_tails = {
        'README.md': '69ac60ee0679e7587adec3e75f884de91b79700d8a24b13ecdef4ab44363ce7f',
        'research_direction.md': 'b06d7acbb9e2308e491ac9cbdea1f48400083abb454ff359031b2737d17231c2',
        'RESEARCH_STATE.md': '15b1d1d926f6c4de40e891a7712308c904b4a4a4fd4a7c99578d0d7c3db68838',
        'ROADMAP.md': '0572fd4469c75ffa86d4e08162ca11784dfdbaa983762c520076a5e152f10fa4',
    }
    assert all(len(value) == 64 for value in expected_tails.values())
    actual_tails = {}
    for name, marker in markers.items():
        data = (ROOT / name).read_bytes()
        actual_tails[name] = sha256(data[data.index(marker.encode('utf-8')):]).hexdigest()
    assert actual_tails == expected_tails, 'Protected history or original milestone changed.'
    stable = {
        'accepted': True, 'date': '2026-10-09',
        'new_scientific_groups': 0, 'new_empirical_groups': 0, 'new_cognitive_axioms': 0,
        'latest_accepted_round': 1062, 'cumulative_scientific_calibrations': 3838,
        'stage_groups': 17, 'whole_roadmap_completed': False,
        'assets_sha256': {str(p.relative_to(PHASE)).replace('\\', '/'): digest(p) for p in owned},
        'history_sha256': {
            '1062/mainline_acceptance.json': digest(base_path),
            '_shared/neutron_star_tidal_adoption.md': digest(HERE / 'neutron_star_tidal_adoption.md'),
            '_shared/neutrino_propagation_adoption.md': digest(HERE / 'neutrino_propagation_adoption.md'),
        },
        'protected_tails_sha256': actual_tails,
        'local_links_verified': links,
    }
    if freeze:
        with RECEIPT.open('x', encoding='utf-8', newline='\n') as handle:
            json.dump(stable, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
    else:
        assert json.loads(RECEIPT.read_text(encoding='utf-8')) == stable
    print(json.dumps({'stellar_interfaces_accepted': True, 'assets': len(owned),
                      'local_links': links, 'new_scientific_groups': 0,
                      'whole_roadmap_completed': False}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze-exclusive', action='store_true')
    main(parser.parse_args().freeze_exclusive)
