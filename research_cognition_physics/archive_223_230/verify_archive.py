"""Audit the renamed archive and integrated v1.1 without rewriting old evidence.

Old manifests still identify the preserved v1.0 paper. Their science/support
entries are checked unchanged; only that paper and the old README use explicit
backup mappings. --run-tests reruns the original eight suites.
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

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
PAPER=BASE/'可组合认知结构与有限维量子理论_阶段论文.md'
REPORT=HERE/'integration_checks.json'
MODULES={223:'reversible_dynamics_bridge',224:'tensor_process_bridge',
         225:'global_orientation_bridge',226:'steering_permission_bridge',
         227:'reversible_control_bridge',228:'instrument_completion_bridge',
         229:'continuous_seed_bridge',230:'finite_protocol_closure'}
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path): return json.loads(path.read_text(encoding='utf-8'))
def check(entry,path):
    assert path.exists(),str(path)
    assert path.stat().st_size==entry['bytes'],str(path)
    assert digest(path)==entry['sha256'],str(path)


def verify(run_tests=False,pending_report=False):
    assert HERE.name=='archive_223_230'
    assert not (BASE/'archive_223_').exists()
    assert not (BASE/'archive_229_').exists()
    original=read(HERE/'STAGE2_MANIFEST.json')
    addition=read(HERE/'STAGE2_CLOSURE_ADDENDUM.json')
    assert original['status']==addition['status']=='closed'
    assert addition['rounds']==list(MODULES)
    assert digest(HERE/'STAGE2_MANIFEST.json')==addition['original_manifest_sha256']
    science=original['research_files']+addition['new_research_files']
    assert len(science)==32
    frozen=science+original['support_files']+addition['support_files']
    for entry in frozen:
        p=(HERE/entry['path']).resolve()
        assert p.is_relative_to(HERE)
        check(entry,p)
    old_paper=HERE/'paper_versions/finite_quantum_v1.0.md'
    check(original['paper'],old_paper)
    snapshot=read(HERE/'PRE_RENAME_FILES.json')
    assert snapshot['new_directory']==HERE.name and len(snapshot['files'])==47
    for entry in snapshot['files']:
        p=HERE/('paper_versions/archive_readme_before_rename.md'
                if entry['path']=='README.md' else entry['path'])
        check(entry,p)
    for p in HERE.glob('*.py'): ast.parse(p.read_text(encoding='utf-8'))
    import integrate_stage2_paper
    assert PAPER.read_text(encoding='utf-8')==integrate_stage2_paper.build()
    assert digest(BASE/'可组合认知结构与复量子状态空间_阶段论文.md')==original['stage1_paper_sha256']
    math=read(HERE/'integrated_math_checks.json')
    assert math['all_passed'] and len(math['documents'])==1
    assert math['documents'][0]['sha256']==digest(PAPER)
    old_math=read(HERE/'stage2_math_review_checks.json')
    assert old_math['all_passed'] and old_math['paper_sha256']==digest(old_paper)
    closing_math=read(HERE/'completion_math_checks.json')
    assert closing_math['all_passed']
    for entry in closing_math['documents']: assert digest(HERE/entry['path'])==entry['sha256']
    spec=importlib.util.spec_from_file_location('stage1_archive_verifier',BASE/'archive_001_222/verify_stage1.py')
    stage1=importlib.util.module_from_spec(spec);spec.loader.exec_module(stage1)
    previous=stage1.verify()
    counts={}
    for number,module in MODULES.items():
        saved=read(HERE/f'{module}_results.json')['checks']
        assert saved['failures']==saved['errors']==0
        counts[number]=saved['run']
        checks=read(HERE/f'research_round_{number}_checks.json')
        assert checks['all_reported_checks_passed']
        for name,sha in checks.get('new_file_hashes',{}).items():
            if Path(name).name=='README.md': continue
            p=HERE/name
            if not p.exists(): p=BASE/name
            assert digest(p)==sha,(number,name)
    assert sum(counts.values())==70
    current={'performed':False,'saved_tests_verified':70}
    if run_tests:
        suite=unittest.TestSuite()
        for module in MODULES.values():
            suite.addTests(unittest.defaultTestLoader.loadTestsFromModule(importlib.import_module(module)))
        output=io.StringIO()
        result=unittest.TextTestRunner(stream=output).run(suite)
        assert result.wasSuccessful() and result.testsRun==70,output.getvalue()
        current={'performed':True,'run':result.testsRun,'failures':len(result.failures),
                 'errors':len(result.errors),'skipped':len(result.skipped)}
    documents=[PAPER,HERE/'README.md',HERE/'STAGE2_ADDENDUM.md',
               HERE/'round229_initial_index.md',HERE/'paper_versions/README.md',
               BASE/'README.md',BASE/'research_direction.md',BASE/'RESEARCH_STATE.md',
               BASE.parent/'README.md']+[HERE/f'research_note_{n}.md' for n in MODULES]
    documents+=list((BASE/'archive_231_').glob('*.md'))
    links=0
    historical_link_mappings=[]
    for p in documents:
        source=p.read_text(encoding='utf-8')
        assert not any(ord(c)<32 and c not in '\r\n\t' for c in source),str(p)
        for target in stage1.local_links(source):
            dest=(p.parent/target).resolve()
            if p==HERE/'round229_initial_index.md' and target=='../archive_223_/README.md':
                dest=HERE/'README.md'
                historical_link_mappings.append({'source':p.name,'old_target':target,
                                                 'current_target':'README.md'})
            if pending_report and dest==REPORT and not dest.exists(): continue
            assert dest.exists(),(str(p),target)
            links+=1
    # Tests must not modify the evidence that was verified before execution.
    for entry in frozen: check(entry,HERE/entry['path'])
    check(original['paper'],old_paper)
    return {'date':'2026-09-22','archive':'archive_223_230','rounds':list(MODULES),
            'paper_version':'v1.1','paper_sha256':digest(PAPER),
            'integration_reproducible_from_v1_0':True,'v1_0_sha256':digest(old_paper),
            'old_manifest_paper_mapping':'paper_versions/finite_quantum_v1.0.md',
            'old_readme_mapping':'paper_versions/archive_readme_before_rename.md',
            'pre_rename_versions_verified':len(snapshot['files']),
            'frozen_scientific_files_verified':len(science),
            'old_manifests_and_scientific_results_rewritten':False,
            'stage1_original_versions_verified':previous['original_research_versions_preserved_and_verified'],
            'stage1_c0_artifacts_verified':previous['original_c0_artifact_hashes_verified'],
            'stage1_paper_unchanged':True,'stage1_science_tests_rerun':False,
            'runtime':{'python':platform.python_version(),'numpy':np.__version__},
            'saved_test_counts':counts,'fresh_regression':current,
            'formula_review':math,'markdown_files_checked':len(documents),
            'local_links_checked':links,'broken_current_links':0,
            'historical_link_mappings':historical_link_mappings,
            'historical_backup_links_not_rewritten':True,
            'top_level_files':previous['top_level_files'],
            'current_stage':'Round231 starts physical-generation research; old stopping instructions are historical.',
            'all_reported_checks_passed':True}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-tests',action='store_true')
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args()
    report=verify(args.run_tests,args.write_checks)
    if args.write_checks:
        REPORT.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='formula_review'},ensure_ascii=False,indent=2))
