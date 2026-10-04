"""Replay unchanged historical scripts in a verified, isolated original layout."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import zipfile

HERE=Path(__file__).resolve().parent
ARCHIVE=HERE.parent
LAYOUT=ARCHIVE.parent
ROOT=LAYOUT.parent


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(script=None, args=(), suite=False, science_followup=False):
    manifest=json.loads((HERE/'manifest.json').read_text('utf8'))
    snapshot=json.loads((HERE/'snapshot_checks.json').read_text('utf8'))
    zip_path=HERE/snapshot['file']
    assert digest(zip_path)==snapshot['sha256'], 'Snapshot changed'
    entries=manifest['entries']
    by_old={e['original']:e for e in entries}
    by_new={e['destination']:e for e in entries}
    code_count=0
    for e in entries:
        if Path(e['original']).suffix in ('.py','.mjs','.json'):
            assert digest(LAYOUT/e['destination'])==e['original_sha256'], e['destination']
            code_count+=1
    if science_followup:
        jobs=[('local_halving_coordinate_audit.py',['--check']),
              ('thermal_direction_bridge_model.py',['--check'])]
    elif suite:
        jobs=[('verify_growth_rounds.py',['240']),
              ('verify_local_halving_round.py',[]),
              ('verify_thermal_direction_bridge.py',[]),
              ('joint_offdiagonal_haar_certificate.py',[]),
              ('verify_round741.py',[]),
              ('verify_round775.py',[]),
              ('postcheck_round775.py',[])]
    else:
        logical=by_new[script]['original'] if script in by_new else script
        assert logical in by_old and logical.endswith('.py'), logical
        jobs=[(logical,list(args))]
    runtime=ROOT/'.research_runtime'
    runtime.mkdir(exist_ok=True)
    reports=[]
    with tempfile.TemporaryDirectory(prefix='phase_replay_231_775_',dir=runtime) as temporary:
        workspace=Path(temporary).resolve()
        # The cleanup target is a newly-created, verified child of our workspace runtime.
        assert workspace.is_relative_to(runtime.resolve()) and workspace!=runtime.resolve()
        started=time.time()
        with zipfile.ZipFile(zip_path) as z:
            for item in snapshot['entries']:
                target=(workspace/item['path']).resolve()
                assert target.is_relative_to(workspace)
                raw=z.read(item['path'])
                assert hashlib.sha256(raw).hexdigest()==item['sha256'],item['path']
                target.parent.mkdir(parents=True,exist_ok=True)
                target.write_bytes(raw)
        original=workspace/'research_cognition_physics/archive_231_'
        print(f'Restored and verified {snapshot["files"]} original files in {time.time()-started:.1f}s',flush=True)
        for name,arguments in jobs:
            selected=original/name
            assert selected.resolve().is_relative_to(original.resolve())
            assert digest(selected)==by_old[name]['original_sha256']
            start=time.time()
            result=subprocess.run([sys.executable,'-B','-X','utf8',str(selected),*arguments],
                                  cwd=selected.parent,text=True,encoding='utf8',
                                  stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=240)
            reports.append(dict(script=name,arguments=arguments,returncode=result.returncode,
                                seconds=round(time.time()-start,3),stdout=result.stdout,stderr=result.stderr,
                                original_code_sha256=digest(selected),passed=result.returncode==0))
            print(f'{name}: {"PASS" if result.returncode==0 else "FAIL"} ({reports[-1]["seconds"]}s)',flush=True)
            if not suite and not science_followup:
                print(result.stdout)
                if result.stderr:
                    print(result.stderr,file=sys.stderr)
        result=dict(date='2026-10-04',original_snapshot_verified=True,
                    restored_files=snapshot['files'],current_code_and_json_files_unchanged=code_count,
                    run_in_isolated_original_layout=True,original_runtime_reused=True,
                    original_evidence_not_overwritten=True,full_545_rounds_rerun=False,
                    jobs=reports,all_checks_passed=all(r['passed'] for r in reports))
    if suite or science_followup:
        target=HERE/('replay_science_followup_checks.json' if science_followup else 'replay_checks.json')
        with target.open('x',encoding='utf8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
    assert result['all_checks_passed'], 'See replay_checks.json for exact failed checks'
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--script')
    p.add_argument('--suite',action='store_true')
    p.add_argument('--science-followup',action='store_true')
    p.add_argument('arguments',nargs=argparse.REMAINDER)
    a=p.parse_args()
    assert sum((bool(a.script),a.suite,a.science_followup))==1
    arguments=a.arguments[1:] if a.arguments[:1]==['--'] else a.arguments
    run(a.script,arguments,a.suite,a.science_followup)
