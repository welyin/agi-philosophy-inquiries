"""Read-only checks for two adopted response relations and an M3C scope correction."""
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
RECEIPT = HERE / 'response_resources_after1062_checks.json'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main(freeze=False):
    if freeze:
        assert not RECEIPT.exists(), 'Do not overwrite a frozen receipt.'
    base_path = PHASE / '1062/mainline_acceptance.json'
    base = json.loads(base_path.read_text(encoding='utf-8'))
    assert base['accepted'] and base['round'] == 1062
    assert base['cumulative_scientific_calibrations'] == 3838
    assert base['stage_groups'] == 17 and not base['whole_roadmap_completed']
    owned = [HERE / 'integration_response_resources_after1062.md', Path(__file__)]
    for group in ('qcd_axion_common_response', 'classical_gravity_interferometer'):
        assets = [HERE / (group + '_adoption.md')]
        assets += [HERE / group / name for name in ('sources.md', 'check.py', 'results.json')]
        review = HERE / (group + '_adoption_review.md')
        review_text = review.read_text(encoding='utf-8')
        assert '通过' in review_text and all(digest(p) in review_text for p in assets)
        owned.extend(assets + [review])
        run = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(HERE / group / 'check.py')],
                             capture_output=True, text=True, encoding='utf-8', check=True)
        print(run.stdout.strip())

    scope = HERE / 'm3c_scope_after1062.md'
    scope_review = HERE / 'm3c_scope_after1062_review.md'
    text = scope_review.read_text(encoding='utf-8')
    assert '通过' in text and digest(scope) in text
    owned += [scope, scope_review]
    # Verify the precise historical receipts quoted by the scope audit;
    # do not reframe those old scientific checks as new experiments.
    old = {
        'archive_956_989/980/research_round_980_checks.json':
            'cf9aaba7348cacdd4cd906c7cebdf1959e8be13d90431ed8121aa04e89c2fefe',
        'archive_990_1008/1002/research_round_1002_checks.json':
            'ceebdc9af4ae627743ca01a1e2d3053db6f2848964802b5cbc4979381dad4475',
        'archive_1044_/1045/research_round_1045_checks.json':
            '86db9340df6869855dc18e4b945d71dc16497e5a2862fadd354495e3876d6236',
        'archive_1044_/1045/mainline_acceptance.json':
            '661b42b37b45f74b6a0085fe899d1386ed6cf5c05f716064757909820428f1c6',
        'archive_1044_/1045/independent_review.md':
            '3c6baf4566c9e76bd7b7308e58fa12bfbb77f1aa92f470ff13b27679f51239f7',
    }
    scope_text = scope.read_text(encoding='utf-8')
    for name, expected in old.items():
        assert len(expected) == 64 and digest(ROOT / name) == expected
        assert expected in scope_text

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

    prior = HERE / 'stellar_interfaces_after1062_checks.json'
    expected_tails = json.loads(prior.read_text(encoding='utf-8'))['protected_tails_sha256']
    markers = {
        'README.md': '## 当前入口：1044—1045阶段已结项',
        'research_direction.md': '以下保存上一阶段结项及历史方向，不覆盖当前目标。',
        'RESEARCH_STATE.md': '以下为1044—1045结项快照和更早历史，不覆盖本栏。',
        'ROADMAP.md': '## M1 检验认知要求的独立选择力',
    }
    tails = {}
    for name, marker in markers.items():
        data = (ROOT / name).read_bytes()
        tails[name] = sha256(data[data.index(marker.encode('utf-8')):]).hexdigest()
    assert tails == expected_tails, 'Original standards or protected history changed.'
    navs = [ROOT / name for name in markers] + [PHASE / 'README.md', PHASE / '文件索引.md']
    for path in navs:
        assert 'integration_response_resources_after1062.md' in path.read_text(encoding='utf-8')
    stable = {
        'accepted': True, 'date': '2026-10-09',
        'new_scientific_groups': 0, 'new_empirical_groups': 0, 'new_cognitive_axioms': 0,
        'latest_accepted_round': 1062, 'cumulative_scientific_calibrations': 3838,
        'stage_groups': 17, 'whole_roadmap_completed': False,
        'assets_sha256': {str(p.relative_to(PHASE)).replace('\\', '/'): digest(p) for p in owned},
        'historical_scope_receipts_sha256': old,
        'baseline_sha256': {'1062/mainline_acceptance.json': digest(base_path),
                            '_shared/stellar_interfaces_after1062_checks.json': digest(prior)},
        'protected_tails_sha256': tails,
        'local_links_verified': links, 'navigation_entries_verified': len(navs),
    }
    if freeze:
        with RECEIPT.open('x', encoding='utf-8', newline='\n') as handle:
            json.dump(stable, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
    else:
        assert json.loads(RECEIPT.read_text(encoding='utf-8')) == stable
    print(json.dumps({'response_and_resource_interfaces_accepted': True,
                      'assets': len(owned), 'local_links': links,
                      'new_scientific_groups': 0, 'whole_roadmap_completed': False}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze-exclusive', action='store_true')
    main(parser.parse_args().freeze_exclusive)
