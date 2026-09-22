"""Read-only recomputation and evidence audit for round 251 or 252."""
import argparse
import ast
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
CONFIG = {251: ('quantum_controlled_birth_audit', 10), 252: ('record_preserving_feedback_audit', 12)}


def read(p):
    return json.loads(p.read_text(encoding='utf-8'))


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify(round_number, pending_report=False):
    stem, count = CONFIG[round_number]
    target = HERE/f'research_round_{round_number}_checks.json'
    names = [stem+'.py', stem+'_results.json', f'research_note_{round_number}.md']
    if target.exists():
        for name, sha in read(target)['new_file_hashes'].items():
            assert digest(HERE/name) == sha, name
    old_frozen = 0
    for old in range(231, round_number):
        for name, sha in read(HERE/f'research_round_{old}_checks.json')['new_file_hashes'].items():
            assert digest(HERE/name) == sha, name
            old_frozen += 1
    spec = importlib.util.spec_from_file_location(stem, HERE/(stem+'.py'))
    science = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(science)
    ast.parse((HERE/(stem+'.py')).read_text(encoding='utf-8'))
    saved = read(HERE/(stem+'_results.json'))
    computed = json.loads(json.dumps(science.report()))
    assert {k:v for k,v in saved.items() if k not in ('checks', 'runtime')} == {
        k:v for k,v in computed.items() if k != 'runtime'}
    output = io.StringIO()
    checks = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromModule(science))
    assert checks.wasSuccessful() and checks.testsRun == count, output.getvalue()
    assert saved['checks'] == {'run': count, 'failures': 0, 'errors': 0}
    math = read(HERE/f'round{round_number}_math_checks.json')
    assert math['all_passed']
    for document in math['documents']:
        assert digest(HERE/document['path']) == document['sha256']
    spec = importlib.util.spec_from_file_location('old_link_parser', BASE/'archive_001_222/verify_stage1.py')
    links_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(links_module)
    links = 0
    for name in [f'research_note_{round_number}.md', 'causal_set_bridge_review.md']:
        p = HERE/name
        text = p.read_text(encoding='utf-8')
        assert not any(ord(c)<32 and c not in '\r\n\t' for c in text)
        for link in links_module.local_links(text):
            dest = (p.parent/link).resolve()
            if pending_report and dest == target and not dest.exists():
                continue
            assert dest.exists(), (name, link)
            links += 1
    return {'date': '2026-09-22', 'round': round_number,
            'fresh_tests': {'run': count, 'failures': 0, 'errors': 0},
            'saved_results_reproduced': True, 'scientific_results_rewritten': False,
            'previous_scientific_file_hashes_verified': old_frozen,
            'math_documents_checked': len(math['documents']), 'local_links_checked': links,
            'broken_links': 0, 'new_file_hashes': {name:digest(HERE/name) for name in names},
            'all_reported_checks_passed': True,
            'scope': 'Finite state-dependent growth and common-future feedback protocols. Classical-mixture limitation, protected-record algebra, disturbance and legal-order compatibility are checked. No autonomous local growth law, continuum or gravity is derived.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round', type=int, choices=CONFIG)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.round, args.write_checks)
    if args.write_checks:
        (HERE/f'research_round_{args.round}_checks.json').write_text(
            json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
