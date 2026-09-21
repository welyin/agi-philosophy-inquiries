"""Read-only verification of stage 223--230 and its two frozen evidence batches.

--run-tests reruns all eight suites, without invoking scientific result writers.
--write-checks may write the completion audit report, but never freezes new hashes
or replaces any research result, paper, original manifest or round report.
"""
import argparse
import ast
import hashlib
import importlib
import io
import json
from pathlib import Path
import platform
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
NEW = {229:'continuous_seed_bridge', 230:'finite_protocol_closure'}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def verify(run_tests=False):
    manifest = read_json(HERE/'STAGE2_CLOSURE_ADDENDUM.json')
    assert manifest['status'] == 'closed'
    assert manifest['rounds'] == list(range(223,231))
    assert not (BASE/'archive_229_').exists()
    assert digest(HERE/'STAGE2_MANIFEST.json') == manifest['original_manifest_sha256']
    for entry in manifest['new_research_files']+manifest['support_files']:
        target = (HERE/entry['path']).resolve()
        assert target.is_relative_to(HERE)
        assert target.stat().st_size == entry['bytes'], entry['path']
        assert digest(target) == entry['sha256'], entry['path']
    for name in (*NEW.values(), 'verify_stage2_completion'):
        ast.parse((HERE/f'{name}.py').read_text(encoding='utf-8'))

    # This old verifier is frozen; its historical policy text is not a new task.
    import verify_stage2
    previous = verify_stage2.verify(run_tests=run_tests)
    assert previous['all_reported_checks_passed']
    saved = {}
    for number, module in NEW.items():
        data = read_json(HERE/f'{module}_results.json')
        check = read_json(HERE/f'research_round_{number}_checks.json')
        assert check['all_reported_checks_passed']
        assert data['checks']['failures'] == data['checks']['errors'] == 0
        assert check['new_tests'] == data['checks']
        for name, sha in check['new_file_hashes'].items():
            assert digest(HERE/name) == sha, name
        saved[number] = data['checks']
    assert sum(c['run'] for c in saved.values()) == 15
    current = {'performed':False, 'saved_tests':70}
    if run_tests:
        suite = unittest.TestSuite()
        for name in NEW.values():
            suite.addTests(unittest.defaultTestLoader.loadTestsFromModule(importlib.import_module(name)))
        output = io.StringIO()
        result = unittest.TextTestRunner(stream=output,verbosity=1).run(suite)
        assert result.wasSuccessful(), output.getvalue()
        assert result.testsRun == 15
        current = {'performed':True,'modules':8,'run':previous['fresh_regression']['run']+result.testsRun,
                   'failures':0,'errors':0,'skipped':len(result.skipped),
                   'scientific_results_rewritten':False}
    math = read_json(HERE/'completion_math_checks.json')
    assert math['all_passed']
    for document in math['documents']:
        assert digest(HERE/document['path']) == document['sha256']

    sys.path.insert(0,str(BASE/'archive_001_222'))
    from verify_stage1 import local_links
    markdown = [HERE/'README.md',HERE/'research_note_229.md',HERE/'research_note_230.md',
                HERE/'STAGE2_ADDENDUM.md',HERE/'round229_initial_index.md',
                BASE/'README.md',BASE/'research_direction.md',BASE/'RESEARCH_STATE.md',BASE.parent/'README.md']
    links = 0
    for path in markdown:
        source = path.read_text(encoding='utf-8')
        assert not any(ord(c)<32 and c not in '\r\n\t' for c in source), str(path)
        for target in local_links(source):
            links += 1
            assert (path.parent/target).resolve().exists(), (path.name,target)
    for entry in manifest['new_research_files']+manifest['support_files']:
        assert digest(HERE/entry['path']) == entry['sha256'], entry['path']
    return {'date':'2026-09-21','status':'closed','rounds':list(range(223,231)),
            'user_instruction':'All closing work is in archive_223_; report completion and await user direction conjectures. No next research round started.',
            'runtime':{'python':platform.python_version(),'numpy':np.__version__},
            'current_tests':current,'new_round_saved_tests':saved,
            'original_stage2_files_verified':previous['frozen_research_files_verified'],
            'new_frozen_research_files_verified':len(manifest['new_research_files']),
            'stage1_original_versions_verified':previous['stage1_original_versions_verified'],
            'stage1_c0_artifacts_verified':previous['stage1_c0_artifacts_verified'],
            'original_stage1_and_stage2_papers_unchanged':True,
            'old_numerical_results_rewritten':False,'stage1_numerical_suite_rerun':False,
            'formula_documents':[{k:v for k,v in x.items() if k!='screenshot'} for x in math['documents']],
            'markdown_files_checked':len(markdown),'local_links_checked':links,'broken_links':0,
            'top_level_files':previous['top_level_files'],
            'open_problem':'Bare F+U+C+P may or may not force a nonconstant reversible path; no complete countermodel or proof is claimed.',
            'all_reported_checks_passed':True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-tests',action='store_true')
    parser.add_argument('--write-checks',action='store_true')
    args = parser.parse_args()
    report = verify(args.run_tests)
    if args.write_checks:
        (HERE/'stage2_completion_checks.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
