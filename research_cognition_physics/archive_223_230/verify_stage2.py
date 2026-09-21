"""Read-only integrity audit for the closed research stage, rounds 223--228.

--run-tests explicitly reruns the six suites without invoking result writers.
--write-checks updates only the phase closure report, never scientific results.
"""
import argparse
import ast
import hashlib
import importlib
import importlib.util
import io
import json
from pathlib import Path
import platform
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
MODULES = {
    223: 'reversible_dynamics_bridge',
    224: 'tensor_process_bridge',
    225: 'global_orientation_bridge',
    226: 'steering_permission_bridge',
    227: 'reversible_control_bridge',
    228: 'instrument_completion_bridge',
}


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(run_tests=False):
    manifest = read_json(HERE/'STAGE2_MANIFEST.json')
    assert manifest['status'] == 'closed'
    assert manifest['rounds'] == list(MODULES)
    expected = {name for n, module in MODULES.items() for name in
                (f'research_note_{n}.md', f'{module}.py', f'{module}_results.json',
                 f'research_round_{n}_checks.json')}
    assert {x['path'] for x in manifest['research_files']} == expected
    for entry in manifest['research_files']+manifest['support_files']+[manifest['paper']]:
        path = (HERE/entry['path']).resolve()
        assert path.is_relative_to(BASE), entry['path']
        assert path.stat().st_size == entry['bytes'], entry['path']
        assert digest(path) == entry['sha256'], entry['path']
    for module in MODULES.values():
        ast.parse((HERE/f'{module}.py').read_text(encoding='utf-8'))

    spec = importlib.util.spec_from_file_location('stage1_archive_verifier', BASE/'archive_001_222/verify_stage1.py')
    stage1 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(stage1)
    previous = stage1.verify()
    assert digest(BASE/'可组合认知结构与复量子状态空间_阶段论文.md') == manifest['stage1_paper_sha256']

    # Check each historical round's science hashes; old README hashes are snapshots.
    historical_science_hashes = set()
    for n, module in MODULES.items():
        checks = read_json(HERE/f'research_round_{n}_checks.json')
        assert checks['all_reported_checks_passed']
        for name, sha in checks.get('new_file_hashes', {}).items():
            if Path(name).name == 'README.md':
                continue
            path = HERE/name
            if not path.exists():
                path = BASE/name
            assert digest(path) == sha, (n, name)
            historical_science_hashes.add(path.name)

    paper = (HERE/manifest['paper']['path']).resolve()
    markdown = [paper, HERE/'README.md', BASE/'README.md', BASE/'research_direction.md',
                BASE/'RESEARCH_STATE.md']+[HERE/f'research_note_{n}.md' for n in MODULES]
    link_count = 0
    for path in markdown:
        source = path.read_text(encoding='utf-8')
        assert not any(ord(c)<32 and c not in '\r\n\t' for c in source), str(path)
        for target in stage1.local_links(source):
            link_count += 1
            assert (path.parent/target).resolve().exists(), (path.name,target)

    math = read_json(HERE/'stage2_math_review_checks.json')
    assert math['paper_sha256'] == digest(paper)
    assert math['all_passed']
    stored = {n: read_json(HERE/f'{module}_results.json')['checks'] for n,module in MODULES.items()}
    assert all(c['failures'] == c['errors'] == 0 for c in stored.values())
    expected_tests = sum(c['run'] for c in stored.values())
    assert expected_tests == 55
    current = {'performed': False, 'scope': 'Saved round reports verified; use --run-tests for a fresh six-module run.'}
    if run_tests:
        sys.path.insert(0, str(HERE))
        suite = unittest.TestSuite()
        for module in MODULES.values():
            suite.addTests(unittest.defaultTestLoader.loadTestsFromModule(importlib.import_module(module)))
        output = io.StringIO()
        result = unittest.TextTestRunner(stream=output, verbosity=1).run(suite)
        if not result.wasSuccessful():
            raise AssertionError(output.getvalue())
        assert result.testsRun == expected_tests
        current = {'performed': True, 'modules': len(MODULES), 'run': result.testsRun,
                   'failures': len(result.failures), 'errors': len(result.errors),
                   'skipped': len(result.skipped), 'scientific_results_rewritten': False}
        # A fresh test run must preserve every frozen artifact.
        for entry in manifest['research_files']:
            assert digest(HERE/entry['path']) == entry['sha256'], entry['path']
    return {'date': '2026-09-20', 'stage': 2, 'status': 'closed', 'rounds': list(MODULES),
            'scope': 'Phase synthesis, archive integrity, links, formula review and optional integrated regression; not a new research round or independent peer review.',
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
            'paper': manifest['paper'], 'frozen_research_files_verified': len(expected),
            'historical_scientific_hashes_verified': len(historical_science_hashes),
            'stage1_original_versions_verified': previous['original_research_versions_preserved_and_verified'],
            'stage1_c0_artifacts_verified': previous['original_c0_artifact_hashes_verified'],
            'stage1_paper_unchanged': True, 'stage1_numerical_suite_rerun': False,
            'saved_round_test_counts': stored, 'saved_tests_total': expected_tests,
            'fresh_regression': current, 'formula_review': math,
            'markdown_files_checked': len(markdown), 'local_links_checked': link_count,
            'new_broken_links': 0, 'top_level_files': previous['top_level_files'],
            'open_problem_status': 'Minimal regularity is retained as an open problem, not an active round229 task.',
            'all_reported_checks_passed': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-tests', action='store_true')
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    report = verify(args.run_tests)
    if args.write_checks:
        (HERE/'stage2_closure_checks.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in {'formula_review','saved_round_test_counts'}}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
