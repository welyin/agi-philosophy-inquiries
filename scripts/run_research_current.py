"""Read-only research replay using files in the current stage directories."""
import argparse
import importlib
import json
import os
from pathlib import Path
import sys
import time
import unittest

from research_layout import Layout, ResearchRuntime, ROOT, no_archive_or_research_writes, test_modules


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--script')
    choice.add_argument('--stage', choices=('1', '2'))
    parser.add_argument('--module', action='append')
    parser.add_argument('--report', type=Path)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.report and args.report.exists():
        parser.error('Choose a new report path')
    store = Layout()
    runtime = ResearchRuntime(store, 'early' if args.stage else 'late')
    sys.dont_write_bytecode = True
    sys.path.append(str(ROOT / '.research_runtime'))
    sys.addaudithook(no_archive_or_research_writes)
    start = time.monotonic()
    details = {}
    with runtime.installed():
        if args.script:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)
        else:
            names = args.module or test_modules(store, args.stage)
            suite = unittest.TestSuite()
            for name in names:
                suite.addTests(unittest.defaultTestLoader.loadTestsFromModule(importlib.import_module(name)))
            result = unittest.TextTestRunner(verbosity=1).run(suite)
            details = dict(module_count=len(names), tests_run=result.testsRun,
                           failures=len(result.failures), errors=len(result.errors),
                           skipped=len(result.skipped))
            if not result.wasSuccessful():
                raise SystemExit(1)
    report = dict(source='current_stage_directories', zip_reads_forbidden=True,
                  restored_old_directory=False, research_files_read_only=True,
                  physical_files_read=sorted(store.read_paths),
                  historical_markdown_link_views=sorted(store.restored_markdown),
                  seconds=round(time.monotonic()-start, 3), all_checks_passed=True, **details)
    if args.report:
        # The report is a new maintenance receipt, not a write by the old script.
        if args.report.absolute().is_relative_to(ROOT / 'research_cognition_physics'):
            parser.error('Write this runtime report outside the research evidence tree')
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with args.report.open('x', encoding='utf8') as stream:
            json.dump(report, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    print(json.dumps({k:v for k,v in report.items()
                      if k not in ('physical_files_read', 'historical_markdown_link_views')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
