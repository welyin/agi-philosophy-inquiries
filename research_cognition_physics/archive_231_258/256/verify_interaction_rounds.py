"""Read-only science/text audit for rounds 255 and 256; no image rendering."""
import argparse
import ast
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import re
import unittest

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
CONFIG = {255: ('mutual_proposal_interaction_audit', 10),
          256: ('interaction_record_bell_audit', 10)}


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def text_checks(path):
    content = path.read_text(encoding='utf-8')
    assert not any(ord(c) < 32 and c not in '\r\n\t' for c in content), path
    delimiters = re.findall(r'^\$\$\s*$', content, re.MULTILINE)
    assert len(delimiters) % 2 == 0, path
    tags = [int(x) for x in re.findall(r'\\tag\{(\d+)\}', content)]
    assert tags == list(range(1, len(tags)+1)), (path, tags)
    assert len(delimiters)//2 == len(tags), path
    return {'path': path.name, 'display_formulas': len(tags),
            'numbering_sequential': True, 'delimiters_paired': True,
            'sha256': digest(path)}


def link_parser():
    spec = importlib.util.spec_from_file_location('old_links', BASE/'archive_001_222/verify_stage1.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.local_links


def verify(number, pending=False):
    stem, count = CONFIG[number]
    target = HERE/f'research_round_{number}_checks.json'
    if target.exists():
        for name, sha in read(target)['new_file_hashes'].items():
            assert digest(HERE/name) == sha, name
    old_count = 0
    for previous in range(231, number):
        for name, sha in read(HERE/f'research_round_{previous}_checks.json')['new_file_hashes'].items():
            assert digest(HERE/name) == sha, name
            old_count += 1
    source = HERE/(stem+'.py')
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location(stem, source)
    science = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(science)
    saved = read(HERE/(stem+'_results.json'))
    assert {k: v for k, v in saved.items() if k not in ('checks', 'runtime')} == json.loads(json.dumps(science.report()))
    output = io.StringIO()
    tests = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromModule(science))
    assert tests.wasSuccessful() and tests.testsRun == count, output.getvalue()
    assert saved['checks'] == {'run': count, 'failures': 0, 'errors': 0}
    note = HERE/f'research_note_{number}.md'
    checked = text_checks(note)
    links = 0
    for link in link_parser()(note.read_text(encoding='utf-8')):
        dest = (note.parent/link).resolve()
        if pending and dest == target and not dest.exists():
            continue
        assert dest.exists(), (note.name, link)
        links += 1
    names = [source.name, stem+'_results.json', note.name]
    return {'date': '2026-09-22', 'round': number,
            'fresh_tests': {'run': count, 'failures': 0, 'errors': 0},
            'saved_results_reproduced': True, 'scientific_results_rewritten': False,
            'previous_scientific_file_hashes_verified': old_count,
            'text_checks': checked, 'local_links_checked': links, 'broken_links': 0,
            'visual_rendering_performed': False,
            'visual_check_policy': 'User requested no further article image checks; text checks only.',
            'new_file_hashes': {name: digest(HERE/name) for name in names},
            'all_reported_checks_passed': True,
            'scope': saved['scope']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round', type=int, choices=CONFIG)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.round, args.write_checks)
    if args.write_checks:
        path = HERE/f'research_round_{args.round}_checks.json'
        if path.exists():
            raise RuntimeError('Check report already exists; use read-only mode.')
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
