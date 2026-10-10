"""Read-only round validation. Hash checks are provenance, not a proof audit."""
import hashlib
import importlib.util
import json
import re
from pathlib import Path
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent


def main():
    receipt = json.loads((HERE / 'mainline_acceptance.json').read_text(encoding='utf8'))
    assert receipt['accepted'] is True
    assert receipt['whole_roadmap_completed'] is False
    assert receipt['full_cognitive_contract_counterexample'] is False
    count = 0
    for group in ('author_assets_sha256', 'review_assets_sha256', 'baseline_assets_sha256', 'mainline_assets_sha256'):
        for name, expected in receipt[group].items():
            path = HERE / name
            actual = hashlib.sha256(path.read_bytes()).hexdigest()
            assert actual == expected, name
            count += 1
    spec = importlib.util.spec_from_file_location('round1063_check', HERE / 'check.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result = module.run()
    module.compare(result, json.loads((HERE / 'results.json').read_text(encoding='utf8')))
    other_spec = importlib.util.spec_from_file_location('round1063_independent', HERE / 'independent_check.py')
    other = importlib.util.module_from_spec(other_spec)
    other_spec.loader.exec_module(other)
    independent = other.run()
    module.compare(independent, json.loads((HERE / 'independent_results.json').read_text(encoding='utf8')))
    checked_links = 0
    names = list(receipt['author_assets_sha256']) + list(receipt['review_assets_sha256']) + list(receipt['mainline_assets_sha256'])
    for name in names:
        path = HERE / name
        if path.suffix != '.md':
            continue
        for target in re.findall(r'\]\(([^\n)]+)\)', path.read_text(encoding='utf8')):
            if '://' in target or target.startswith('#'):
                continue
            target = unquote(target.split('#', 1)[0].strip('<>'))
            assert (path.parent / target).exists(), (name, target)
            checked_links += 1
    print(json.dumps({'round': 1063, 'passed': True, 'independent_algorithm_passed': independent['passed'], 'asset_hashes': count,
                      'existing_local_links': checked_links, 'maximum_residual': result['maximum_residual'],
                      'scope': receipt['scope']}, ensure_ascii=False))


if __name__ == '__main__':
    main()
