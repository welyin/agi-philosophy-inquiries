"""Check two count0 physical adoptions without changing accepted historical assets."""
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
RECEIPT = HERE / 'antimatter_lnv_after1062_checks.json'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main(freeze=False):
    if freeze:
        assert not RECEIPT.exists(), 'Never replace a frozen receipt.'
    baseline = PHASE / '1062/mainline_acceptance.json'
    base = json.loads(baseline.read_text(encoding='utf-8'))
    assert base['accepted'] and base['round'] == 1062
    assert base['cumulative_scientific_calibrations'] == 3838
    assert base['stage_groups'] == 17 and not base['whole_roadmap_completed']
    prior = HERE / 'response_resources_after1062_checks.json'
    prior_data = json.loads(prior.read_text(encoding='utf-8'))
    assert digest(baseline) == prior_data['baseline_sha256']['1062/mainline_acceptance.json']
    owned = [HERE / 'integration_antimatter_lnv_after1062.md', Path(__file__)]
    for group in ('antihydrogen_gravity', 'neutrino_lnv_matching'):
        assets = [HERE / (group + '_adoption.md')]
        assets += [HERE / group / name for name in ('sources.md', 'check.py', 'results.json')]
        review = HERE / (group + '_adoption_review.md')
        review_text = review.read_text(encoding='utf-8')
        assert '通过' in review_text and all(digest(p) in review_text for p in assets)
        owned.extend(assets + [review])
        run = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(HERE / group / 'check.py')],
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
    tails = {}
    for name, marker in markers.items():
        data = (ROOT / name).read_bytes()
        tails[name] = sha256(data[data.index(marker.encode('utf-8')):]).hexdigest()
    assert tails == prior_data['protected_tails_sha256'], 'Protected history or original standards changed.'
    navs = [ROOT / name for name in markers] + [PHASE / 'README.md', PHASE / '文件索引.md']
    for path in navs:
        assert 'integration_antimatter_lnv_after1062.md' in path.read_text(encoding='utf-8')
    stable = {
        'accepted': True, 'date': '2026-10-09',
        'new_scientific_groups': 0, 'new_empirical_groups': 0, 'new_cognitive_axioms': 0,
        'latest_accepted_round': 1062, 'cumulative_scientific_calibrations': 3838,
        'stage_groups': 17, 'whole_roadmap_completed': False,
        'assets_sha256': {str(p.relative_to(PHASE)).replace('\\', '/'): digest(p) for p in owned},
        'baseline_sha256': {'1062/mainline_acceptance.json': digest(baseline),
                            '_shared/response_resources_after1062_checks.json': digest(prior)},
        'protected_tails_sha256': tails,
        'local_links_verified': links, 'navigation_entries_verified': len(navs),
    }
    if freeze:
        with RECEIPT.open('x', encoding='utf-8', newline='\n') as handle:
            json.dump(stable, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
    else:
        assert json.loads(RECEIPT.read_text(encoding='utf-8')) == stable
    print(json.dumps({'antimatter_and_lnv_interfaces_accepted': True,
                      'assets': len(owned), 'local_links': links,
                      'new_scientific_groups': 0, 'whole_roadmap_completed': False}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze-exclusive', action='store_true')
    main(parser.parse_args().freeze_exclusive)
