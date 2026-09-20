"""Run archived round 1..222 unittest cases without rewriting result files."""
import argparse
import ast
import importlib
import json
from pathlib import Path
import platform
import re
import sys
import time
import unittest

ARCHIVE = Path(__file__).resolve().parent
ROOT = ARCHIVE.parent.parent
LIBRARY = ARCHIVE / 'research_process'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('modules', nargs='*', help='Module names, optionally ending in .py')
    parser.add_argument('--all', action='store_true', help='All archived modules declaring test cases')
    parser.add_argument('--report', type=Path, help='Write a new report outside research_process; the archive root is allowed')
    args = parser.parse_args()
    if args.all and args.modules:
        parser.error('Choose explicit modules or --all')
    if not args.all and not args.modules:
        parser.error('Supply one or more modules, or --all')
    if args.report:
        destination = args.report.resolve()
        if (destination.is_relative_to(LIBRARY.resolve()) or destination.exists()
                or destination.parent == ARCHIVE.parent):
            parser.error('Report must be new, outside research_process and outside the papers-only parent directory')
    sys.dont_write_bytecode = True
    sys.path[:0] = [str(LIBRARY), str(ROOT / '.research_runtime')]
    if args.all:
        names = []
        for path in sorted(LIBRARY.glob('*.py')):
            tree = ast.parse(path.read_text(encoding='utf-8-sig'))
            if any(isinstance(n, ast.ClassDef) and any(
                    (isinstance(b, ast.Attribute) and b.attr == 'TestCase') or
                    (isinstance(b, ast.Name) and b.id == 'TestCase') for b in n.bases)
                   for n in tree.body):
                names.append(path.stem)
    else:
        names = [n.removesuffix('.py') for n in args.modules]
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    counts = {}
    for name in names:
        if not re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*', name) or not (LIBRARY / (name+'.py')).is_file():
            parser.error('Unknown archived module: '+name)
        tests = loader.loadTestsFromModule(importlib.import_module(name))
        counts[name] = tests.countTestCases()
        suite.addTests(tests)
    start = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=1).run(suite)
    import numpy as np
    report = {'scope':'Archived unittest cases; no main functions or result writers invoked.',
              'python':platform.python_version(), 'numpy':np.__version__, 'modules':counts,
              'module_count':len(counts), 'tests_run':result.testsRun,
              'failures':len(result.failures), 'errors':len(result.errors),
              'skipped':len(result.skipped), 'seconds':round(time.perf_counter()-start, 3),
              'success':result.wasSuccessful()}
    if args.report:
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='modules'}, ensure_ascii=False, indent=2))
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__ == '__main__':
    main()
