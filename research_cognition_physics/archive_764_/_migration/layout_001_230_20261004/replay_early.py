"""Run early-stage tests directly from the current physical stage directories."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[4]
RUNNER = ROOT / 'scripts/run_research_current.py'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('1', '2', 'all'), default='all')
    parser.add_argument('--module', action='append')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    if args.report and args.report.exists():
        parser.error('Choose a new report path')
    jobs = []
    for stage in ('1', '2') if args.stage == 'all' else (args.stage,):
        command = [sys.executable, '-B', '-X', 'utf8', str(RUNNER), '--stage', stage]
        if stage == '1':
            for name in args.module or []:
                command += ['--module', name]
        result = subprocess.run(command, cwd=ROOT, text=True, encoding='utf8',
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=600)
        jobs.append(dict(stage=stage, returncode=result.returncode, stdout=result.stdout,
                         stderr=result.stderr, passed=result.returncode == 0))
        print(result.stdout, flush=True)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
    report = dict(source='current_stage_directories', zip_dependency=False,
        restored_old_directory=False, jobs=jobs, all_tests_passed=all(r['passed'] for r in jobs))
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with args.report.open('x', encoding='utf8') as stream:
            json.dump(report, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    if not report['all_tests_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
