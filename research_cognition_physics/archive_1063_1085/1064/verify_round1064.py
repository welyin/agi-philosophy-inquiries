"""Read-only frozen provenance, current links and two algorithm checks."""
from pathlib import Path
import hashlib
import importlib.util
import json
import re
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    receipt = json.loads((HERE / 'mainline_acceptance.json').read_text(encoding='utf8'))
    assert receipt['accepted'] is True
    assert receipt['whole_roadmap_completed'] is False
    assert receipt['full_cognitive_contract_counterexample'] is False
    total, links = 0, 0
    for group in ('author_assets_sha256', 'review_assets_sha256', 'baseline_assets_sha256'):
        for name, digest in receipt[group].items():
            path = HERE / name
            assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, name
            total += 1
            if group == 'baseline_assets_sha256' or path.suffix != '.md':
                continue
            for target in re.findall(r'\]\(([^\n)]+)\)', path.read_text(encoding='utf8')):
                if '://' in target or target.startswith('#'):
                    continue
                local = unquote(target.split('#', 1)[0].strip('<>'))
                assert (path.parent / local).exists(), (name, target)
                links += 1
    author = module('round1064_author', HERE / 'check.py')
    result = author.run()
    author.compare(result, json.loads((HERE / 'results.json').read_text(encoding='utf8')))
    independent = module('round1064_independent', HERE / 'independent_check.py').run()
    author.compare(independent, json.loads((HERE / 'independent_results.json').read_text(encoding='utf8')))
    assert independent['passed']
    print(json.dumps({'round':1064, 'passed':True, 'frozen_assets':total,
                      'local_links':links, 'author_and_independent_recomputed':True,
                      'maximum_author_residual':result['maximum_residual'],
                      'cumulative_scientific_calibrations':3840,
                      'whole_roadmap_completed':False}, ensure_ascii=False))


if __name__ == '__main__':
    main()
