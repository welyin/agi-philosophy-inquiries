"""Replay archived scripts directly from the current stage directories."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
RUNNER = ROOT / 'scripts/run_research_current.py'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--script')
    group.add_argument('--suite', action='store_true')
    group.add_argument('--science-followup', action='store_true')
    parser.add_argument('--report', type=Path)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.report and args.report.exists():
        parser.error('Choose a new report path')
    extra = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
    if args.science_followup:
        jobs = [('local_halving_coordinate_audit.py', ['--check']),
                ('thermal_direction_bridge_model.py', ['--check'])]
    elif args.suite:
        # Publication-time navigation hashes are historical receipts, not current science.
        jobs = [('verify_growth_rounds.py', ['240']),
                ('local_halving_coordinate_audit.py', ['--check']),
                ('thermal_direction_bridge_model.py', ['--check']),
                ('joint_offdiagonal_haar_certificate.py', []),
                ('joint_wick_source_anomaly.py', [])]
    else:
        jobs = [(args.script, extra)]
    reports = []
    for script, arguments in jobs:
        result = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(RUNNER),
            '--script', script, '--', *arguments], cwd=ROOT, text=True, encoding='utf8',
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=600)
        reports.append(dict(script=script, arguments=arguments, returncode=result.returncode,
                            stdout=result.stdout, stderr=result.stderr, passed=result.returncode == 0))
        print(script + ': ' + ('PASS' if result.returncode == 0 else 'FAIL'), flush=True)
        if not result.returncode and len(jobs) == 1:
            print(result.stdout)
        if result.returncode:
            print(result.stdout); print(result.stderr, file=sys.stderr)
    report = dict(source='current_stage_directories', zip_dependency=False,
        restored_old_directory=False, publication_snapshots_required=False,
        jobs=reports, all_checks_passed=all(r['passed'] for r in reports))
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with args.report.open('x', encoding='utf8') as stream:
            json.dump(report, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    if not report['all_checks_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
