"""Read-only replay beneath the reversible 1234 GR-response navigation.

Current 1234 scripts use actual paths. Closed evidence retains its original
logical paths and hash-checked bytes through the established closure reader.
The known historical paper-08 external edit is neither hidden nor rebound.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import types

import run_research_thermal_geometry_closed as previous
from research_layout import ROOT, ResearchLoader, no_archive_or_research_writes, relocation_original_bytes

THIS = Path(__file__).resolve()
STAGE = 'research_cognition_physics/archive_1234_'
LAYER = ROOT / 'research_cognition_physics/_migration/gr_response_active_1234/layer.json'
NAVIGATION = {'research_cognition_physics/' + name for name in
              ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')}


class GRResponseLayout(previous.ClosedThermalGeometryLayout):
    def __init__(self):
        raw = LAYER.read_bytes()
        self.response_layer, self.response_sha256 = json.loads(raw), previous.digest(raw)
        if self.response_layer.get('schema') != 'gr_response_active_navigation_v1':
            raise ValueError('Unexpected GR-response navigation schema')
        self.response_entries = {}
        for entry in self.response_layer['entries']:
            old, new = previous.previous.previous.previous.validate_entry(entry)
            if old not in NAVIGATION or new != old or old in self.response_entries:
                raise ValueError('GR-response navigation may change only four distinct root documents')
            self.response_entries[old] = entry
        if set(self.response_entries) != NAVIGATION:
            raise ValueError('All four root navigation documents are required')
        super().__init__()
        for name, entry in self.response_entries.items():
            boundary = self.geometry_closure_entries[name]
            if (entry['original_sha256'], entry['original_bytes']) != (boundary['current_sha256'], boundary['current_bytes']):
                raise ValueError('GR-response navigation does not join closed geometry boundary: ' + name)

    def before_geometry_closure(self, path):
        name = self.geometry_closure_name(path)
        if previous.inside(name, previous.OLD) and name not in self.geometry_closure_entries:
            raise FileNotFoundError('File absent from thermal-geometry closure: ' + str(name))
        raw = self.geometry_closure_actual(path).read_bytes()
        active = self.response_entries.get(name)
        if active:
            raw = relocation_original_bytes(raw, active)
        entry = self.geometry_closure_entries.get(name)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def script_key(self, path):
        name = previous.relative(path)
        if previous.inside(name, STAGE) and name.endswith('.py') and Path(str(path)).is_file():
            return name
        return super().script_key(path)

    def script_actual(self, key):
        return ROOT / key if previous.inside(key, STAGE) else super().script_actual(key)

    def verify_response_navigation(self):
        if previous.digest(LAYER.read_bytes()) != self.response_sha256:
            raise ValueError('GR-response layer changed during replay')
        for name, entry in self.response_entries.items():
            relocation_original_bytes((ROOT / name).read_bytes(), entry)
        for name, expected in self.response_layer['unchanged_support_sha256'].items():
            if previous.digest((ROOT / previous.safe_relative(name)).read_bytes()) != expected:
                raise ValueError('Frozen navigation support changed: ' + name)
        notice = ROOT / previous.safe_relative(self.response_layer['notice_path'])
        if previous.digest(notice.read_bytes()) != self.response_layer['notice_sha256']:
            raise ValueError('Active notice changed without a reversible navigation update')
        return dict(status='PASS', files=4, manifest_sha256=self.response_sha256,
                    closed_boundary_verified=True, original_navigation_recoverable=True,
                    frozen_archives_modified=False, stage_first_round=1234)

    def verify_layout(self):
        # The inherited exact historical check intentionally still fails on
        # paper 08. Scoped stage completion is separately checked by its code.
        return dict(response_navigation=self.verify_response_navigation(),
                    frozen_layers=super().verify_layout(), all_checks_passed=True)


class GRResponseLoader(ResearchLoader):
    def exec_module(self, module):
        if not previous.inside(self.key, STAGE):
            return previous.ClosedThermalGeometryLoader(self.runtime, self.key).exec_module(module)
        actual = self.runtime.store.script_actual(self.key)
        module.__file__ = str(actual)
        old_path = list(sys.path)
        sys.path.insert(0, str(actual.parent))
        try:
            exec(compile(actual.read_bytes(), str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = old_path


class GRResponseRuntime(previous.ClosedThermalGeometryRuntime):
    def __init__(self, store):
        super().__init__(store)
        old_scoped_import = self.scoped_builtins[previous.VIEW]['__import__']
        old_import_proxy = old_scoped_import('importlib')
        import_proxy = types.ModuleType('importlib')
        import_proxy.__dict__.update(vars(old_import_proxy))
        old_active = previous.previous
        active_proxy = types.ModuleType(old_active.__name__)
        active_proxy.__dict__.update(vars(old_active))
        active_proxy.ThermalGeometryLayout = GRResponseLayout

        def historical_import_module(name, package=None):
            return active_proxy if name == old_active.__name__ and package is None else old_import_proxy.import_module(name, package)

        import_proxy.import_module = historical_import_module

        def scoped_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == 'importlib' and level == 0:
                return import_proxy
            return old_scoped_import(name, globals, locals, fromlist, level)

        self.scoped_builtins[previous.VIEW]['__import__'] = scoped_import

    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, GRResponseLoader(self, key),
                                              origin=str(self.store.script_actual(key)))


@contextmanager
def frozen_runner_routes():
    old_run = subprocess.run

    def mapped(command, *args, **kwargs):
        if isinstance(command, (list, tuple)) and len(command) >= 6:
            same_python = os.path.normcase(os.path.abspath(str(command[0]))) == os.path.normcase(sys.executable)
            target = previous.relative(command[4])
            if same_python and list(command[1:4]) == ['-B', '-X', 'utf8']:
                frozen = ('scripts/run_research_thermal_geometry_closed.py',
                          'scripts/run_research_thermal_geometry_active.py',
                          'scripts/run_research_thermodynamics_closed.py')
                flag = command[5]
                if target in frozen and flag in ('--verify-navigation', '--verify-closure',
                                                  '--verify-layout', '--verify-prior-closure', '--script'):
                    if flag == '--verify-navigation':
                        flag = '--verify-previous-navigation'
                    elif flag == '--verify-closure' and target != frozen[0]:
                        flag = '--verify-prior-closure'
                    command = [*command[:4], str(THIS), flag, *command[6:]]
        return old_run(command, *args, **kwargs)

    subprocess.run = mapped
    try:
        yield
    finally:
        subprocess.run = old_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--script')
    group.add_argument('--verify-layout', action='store_true')
    group.add_argument('--verify-closure', action='store_true')
    group.add_argument('--verify-navigation', action='store_true')
    group.add_argument('--verify-previous-navigation', action='store_true', help=argparse.SUPPRESS)
    group.add_argument('--verify-prior-closure', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = GRResponseLayout()
    runtime = GRResponseRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    routes = previous.previous.previous.previous.combined.merged.previous.previous
    routes.THIS = THIS
    previous.previous.previous.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), routes.execution_routes(previous.ExecutionLayout(store), runtime), previous.previous.previous.gr_verification_route(), frozen_runner_routes():
        if args.script:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)
        else:
            result = (previous.previous.previous.ThermalLayout.verify_closure(store) if args.verify_prior_closure else
                      store.verify_active_navigation() if args.verify_previous_navigation else
                      store.verify_layout() if args.verify_layout else
                      store.verify_closure() if args.verify_closure else store.verify_response_navigation())
            print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
