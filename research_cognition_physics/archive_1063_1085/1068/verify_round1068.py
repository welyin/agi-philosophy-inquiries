"""Read-only verification of round 1068 evidence, formulas, and local links."""
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    item = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(item)
    return item

def main():
    receipt = json.loads((HERE / 'mainline_acceptance.json').read_text(encoding='utf8'))
    assert receipt['accepted'] and not receipt['spatial_generation_goal_completed']
    assets = 0
    for category in ('author_assets_sha256', 'review_assets_sha256', 'baseline_assets_sha256'):
        for relative, expected in receipt[category].items():
            assert hashlib.sha256((HERE / relative).read_bytes()).hexdigest() == expected, relative
            assets += 1
    author = module('author1068', HERE / 'check.py')
    author.compare(author.run(), json.loads((HERE / 'results.json').read_text(encoding='utf8')))
    independent = subprocess.run(
        [sys.executable, '-B', '-X', 'utf8', str(HERE / 'independent_check.py')],
        capture_output=True, text=True, encoding='utf8', check=True)
    body, script_hash = independent.stdout.rsplit('SCRIPT_SHA256 ', 1)
    assert script_hash.strip() == receipt['review_assets_sha256']['independent_check.py']
    assert json.loads(body) == json.loads((HERE / 'independent_results.json').read_text(encoding='utf8'))
    links = 0
    for relative in ('../research_note_1068.md', 'plan.md', 'proof.md', 'source_scope.md',
                     'dependency_update.md', 'independent_review.md', 'mathematical_review.md'):
        path = HERE / relative
        for target in re.findall(r'\]\(([^)]+)\)', path.read_text(encoding='utf8')):
            if re.match(r'^[a-zA-Z]+://', target) or target.startswith('#'):
                continue
            target = unquote(target.split('#')[0].strip('<>'))
            assert (path.parent / target).exists(), (relative, target)
            links += 1
    output = {'round': 1068, 'passed': True, 'frozen_assets': assets, 'local_links': links,
              'author_recomputed': True, 'independent_recomputed': True,
              'spatial_generation_goal_completed': False, 'whole_roadmap_completed': False}
    print(json.dumps(output, ensure_ascii=False))
    return output

if __name__ == '__main__':
    main()
