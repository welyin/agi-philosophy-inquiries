"""Run original early-stage tests in their original, verified directory layout."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import zipfile

HERE=Path(__file__).resolve().parent
MIG=HERE.parent
BASE=MIG.parent.parent
ROOT=BASE.parent


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('1','2','all'),default='all')
    parser.add_argument('--module',action='append',help='Selected stage-1 module, repeatable; otherwise all 218 test modules')
    parser.add_argument('--report',type=Path,help='New receipt path; never overwrites an existing report')
    args=parser.parse_args()
    if args.report:
        assert not args.report.exists(), 'Choose a new report path'
    manifest=json.loads((HERE/'manifest.json').read_text('utf8'))
    code_entries=[e for e in manifest['entries'] if Path(e['original']).suffix in ('.py','.mjs','.json')]
    for e in code_entries:
        assert sha((ROOT/e['destination']).read_bytes())==e['original_sha256'],e['destination']
    snapshot=json.loads((MIG/'snapshot_checks.json').read_text('utf8'))
    frozen=MIG/snapshot['file']
    assert sha(frozen.read_bytes())==snapshot['sha256']
    inventory={e['path']:e for e in snapshot['entries']}
    for e in code_entries:
        assert inventory[e['original']]['sha256']==e['original_sha256'], 'Replay source differs from the migrated original: '+e['original']
    jobs=[]
    if args.stage in ('1','all'):
        jobs.append(('stage1',['research_cognition_physics/archive_001_222/run_stage1.py',*(args.module or ['--all'])],1915 if not args.module else None))
    if args.stage in ('2','all'):
        modules=['reversible_dynamics_bridge','tensor_process_bridge','global_orientation_bridge',
                 'steering_permission_bridge','reversible_control_bridge','instrument_completion_bridge',
                 'continuous_seed_bridge','finite_protocol_closure']
        jobs.append(('stage2',['-m','unittest',*modules],70))
    runtime=ROOT/'.research_runtime'
    runtime.mkdir(exist_ok=True)
    env=dict(os.environ)
    # Reuse the already-installed optional SciPy runtime without copying/installing it.
    env['PYTHONPATH']=str(runtime)+(os.pathsep+env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
    reports=[]
    with tempfile.TemporaryDirectory(prefix='phase_replay_001_230_',dir=runtime) as temporary:
        work=Path(temporary).resolve()
        assert work.is_relative_to(runtime.resolve()) and work!=runtime.resolve()
        with zipfile.ZipFile(frozen) as z:
            for name,e in inventory.items():
                target=(work/name).resolve()
                assert target.is_relative_to(work)
                raw=z.read(name)
                assert sha(raw)==e['sha256'],name
                target.parent.mkdir(parents=True,exist_ok=True)
                target.write_bytes(raw)
        for label,command,expected in jobs:
            cwd=work if label=='stage1' else work/'research_cognition_physics/archive_223_230'
            start=time.monotonic()
            print(f'Running {label} original unittest suites...',flush=True)
            result=subprocess.run([sys.executable,'-B','-X','utf8',*command],cwd=cwd,env=env,
                                  text=True,encoding='utf8',stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=600)
            report=dict(stage=label,command=command,expected_tests=expected,returncode=result.returncode,
                        seconds=round(time.monotonic()-start,3),stdout=result.stdout,stderr=result.stderr)
            if label=='stage1' and result.returncode==0:
                result_data=json.loads(result.stdout)
                report['tests_run']=result_data['tests_run']
            elif label=='stage2':
                import re
                match=re.search(r'Ran (\d+) tests?',result.stderr)
                report['tests_run']=int(match[1]) if match else 0
            report['passed']=result.returncode==0 and (expected is None or report.get('tests_run')==expected)
            reports.append(report)
            print(f'{label}: {"PASS" if report["passed"] else "FAIL"}, tests={report.get("tests_run")}, seconds={report["seconds"]}',flush=True)
            if not report['passed']:
                print(result.stdout[-4000:]);print(result.stderr[-4000:])
    for e in code_entries:
        assert sha((ROOT/e['destination']).read_bytes())==e['original_sha256'],e['destination']
    report=dict(date='2026-10-04',scope='Original unittest suites; no result writers or historical publishing checks',
                python=sys.executable,existing_optional_runtime=str(runtime),frozen_layout_verified=True,
                current_code_results_and_receipts_verified=len(code_entries),jobs=reports,
                all_tests_passed=all(r['passed'] for r in reports),scientific_evidence_unchanged=True)
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True)
        with args.report.open('x',encoding='utf8') as f:
            json.dump(report,f,ensure_ascii=False,indent=2);f.write('\n')
    assert report['all_tests_passed'],'See the report for the original test failure'


if __name__=='__main__':
    main()
