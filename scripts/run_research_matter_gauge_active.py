"""Read-only replay beneath the reversible 1281 matter/gauge navigation.

New-stage scripts see current files. Closed-stage scripts retain their recorded
logical paths and bytes. The new navigation prefix is reversed in memory before
the unchanged 1234--1280 closure is checked. Earlier exceptions remain visible.
"""
from __future__ import annotations

import argparse
import builtins
from contextlib import contextmanager
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import types

import run_research_gr_response_closed as previous
from research_layout import ROOT, ResearchLoader, no_archive_or_research_writes
from research_layout import relocation_original_bytes

THIS = Path(__file__).resolve()
STAGE = 'research_cognition_physics/archive_1281_'
LAYER = ROOT / 'research_cognition_physics/_migration/matter_gauge_active_1281/layer.json'
VIEW = 'matter_gauge_active'
NAVIGATION = {'research_cognition_physics/' + name for name in
              ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')}
ENTRY_POINTS = [STAGE + '/README.md', STAGE + '/_shared/goal.md']
digest, relative, inside = previous.digest, previous.relative, previous.inside


class MatterGaugeLayout(previous.ClosedGRResponseLayout):
    def __init__(self):
        raw = LAYER.read_bytes()
        self.matter_layer, self.matter_sha256 = json.loads(raw), digest(raw)
        if (self.matter_layer.get('schema') != 'matter_gauge_active_navigation_v1'
                or self.matter_layer.get('stage_path') != STAGE
                or self.matter_layer.get('stage_first_round') != 1281
                or self.matter_layer.get('entry_points') != ENTRY_POINTS):
            raise ValueError('Unexpected matter/gauge navigation schema or stage')
        self.matter_entries = {}
        for entry in self.matter_layer['entries']:
            old, new = previous.gr.validate_entry(entry)
            if old not in NAVIGATION or new != old or old in self.matter_entries:
                raise ValueError('Only four distinct root navigation documents may change')
            if (len(entry['edits']) != 1 or entry['edits'][0]['before'] != ''
                    or entry['edits'][0]['start'] != entry['edits'][0]['end']):
                raise ValueError('Matter/gauge navigation must be a single prefix insertion')
            self.matter_entries[old] = entry
        if set(self.matter_entries) != NAVIGATION:
            raise ValueError('Exactly four root navigation documents are required')
        super().__init__()
        for name, entry in self.matter_entries.items():
            boundary = self.response_closure_entries[name]
            if (entry['original_sha256'], entry['original_bytes']) != (
                    boundary['current_sha256'], boundary['current_bytes']):
                raise ValueError('Navigation does not join the closed GR-response boundary: ' + name)

    def before_matter_navigation(self, path):
        name = self.response_closure_name(path)
        raw = self.response_closure_actual(path).read_bytes()
        entry = self.matter_entries.get(name)
        return relocation_original_bytes(raw, entry) if entry else raw

    def before_response_closure(self, path):
        name = self.response_closure_name(path)
        if inside(name, previous.OLD) and name not in self.response_closure_entries:
            raise FileNotFoundError('File absent from GR-response closure: ' + str(name))
        raw = self.before_matter_navigation(path)
        entry = self.response_closure_entries.get(name)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def script_key(self, path):
        name = relative(path)
        if inside(name, STAGE):
            return name if name.endswith('.py') and (ROOT / name).is_file() else None
        return super().script_key(path)

    def script_actual(self, key):
        return ROOT / key if inside(key, STAGE) else super().script_actual(key)

    def verify_matter_navigation(self):
        if digest(LAYER.read_bytes()) != self.matter_sha256:
            raise ValueError('Matter/gauge layer changed during replay')
        for name, entry in self.matter_entries.items():
            old = self.before_matter_navigation(ROOT / name)
            enc = 'utf-8-sig' if old.startswith(b'\xef\xbb\xbf') else 'utf8'
            if entry['edits'][0]['start'] != old.decode(enc).index('\n') + 1:
                raise ValueError('Navigation insertion must follow the original title: ' + name)
        support = self.matter_layer['unchanged_support_sha256']
        required = dict(self.response_closure_manifest['unchanged_support_sha256'])
        required[previous.MANIFEST.relative_to(ROOT).as_posix()] = self.response_closure_sha256
        required[previous.THIS.relative_to(ROOT).as_posix()] = digest(previous.THIS.read_bytes())
        if not isinstance(support, dict) or not set(required).issubset(support):
            raise ValueError('Missing inherited support or closure bindings')
        for name, expected in support.items():
            if digest((ROOT / previous.safe_relative(name)).read_bytes()) != expected:
                raise ValueError('Frozen navigation support changed: ' + name)
            if name in required and expected != required[name]:
                raise ValueError('Inherited support binding changed: ' + name)
        for name in ENTRY_POINTS:
            if not (ROOT / name).is_file():
                raise FileNotFoundError('New-stage navigation target is missing: ' + name)
        return dict(status='PASS', files=4, manifest_sha256=self.matter_sha256,
                    closed_boundary_verified=True, original_navigation_recoverable=True,
                    frozen_archives_modified=False, stage_first_round=1281,
                    entry_points=ENTRY_POINTS, frozen_support_files=len(support),
                    scientific_acceptance_checked=False)

    def verify_layout(self):
        return dict(matter_gauge_navigation=self.verify_matter_navigation(),
                    frozen_layers=super().verify_layout(), all_checks_passed=True)


class MatterGaugeLoader(ResearchLoader):
    def exec_module(self, module):
        if not inside(self.key, STAGE):
            return previous.ClosedGRResponseLoader(self.runtime, self.key).exec_module(module)
        actual = self.runtime.store.script_actual(self.key)
        module.__file__, module.__builtins__ = str(actual), self.runtime.scoped_builtins[VIEW]
        old_path = list(sys.path)
        sys.path.insert(0, str(actual.parent))
        try:
            exec(compile(actual.read_bytes(), str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = old_path


class MatterGaugeRuntime(previous.ClosedGRResponseRuntime):
    def __init__(self, store):
        super().__init__(store)
        # Current scripts keep native paths. Only historical reader imports are
        # adapted, and only inside modules loaded by this runtime.
        self.scoped_builtins[VIEW] = dict(vars(builtins))
        proxies = {}
        for module, layout_name, runtime_name in (
                (previous, 'ClosedGRResponseLayout', 'ClosedGRResponseRuntime'),
                (previous.previous, 'GRResponseLayout', 'GRResponseRuntime'),
                (previous.geometry.previous, 'ThermalGeometryLayout', 'ThermalGeometryRuntime')):
            proxy = types.ModuleType(module.__name__)
            proxy.__dict__.update(vars(module))
            setattr(proxy, layout_name, MatterGaugeLayout)
            setattr(proxy, runtime_name, MatterGaugeRuntime)
            proxies[module.__name__] = proxy
        for bindings in self.scoped_builtins.values():
            fallback = bindings['__import__']
            old_importlib = fallback('importlib')
            library = types.ModuleType('importlib')
            library.__dict__.update(vars(old_importlib))

            def import_module(name, package=None, original=old_importlib.import_module):
                return proxies[name] if package is None and name in proxies else original(name, package)

            library.import_module = import_module

            def scoped_import(name, globals=None, locals=None, fromlist=(), level=0,
                              original=fallback, importlib_proxy=library):
                if level == 0 and name in proxies:
                    return proxies[name]
                if level == 0 and name == 'importlib':
                    return importlib_proxy
                return original(name, globals, locals, fromlist, level)

            bindings['__import__'] = scoped_import

    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, MatterGaugeLoader(self, key),
                                              origin=str(self.store.script_actual(key)))


@contextmanager
def frozen_runner_routes():
    old_run = subprocess.run
    # Retain each older flag's scope rather than promoting it to the new layer.
    modes = {
        'scripts/run_research_gr_response_closed.py': {
            '--verify-navigation': '--verify-response-navigation',
            '--verify-closure': '--verify-closure',
            '--verify-geometry-closure': '--verify-geometry-closure',
            '--verify-previous-navigation': '--verify-previous-navigation',
            '--verify-prior-closure': '--verify-prior-closure',
            '--verify-gr-navigation': '--verify-gr-navigation',
            '--verify-local-closure': '--verify-local-closure'},
        'scripts/run_research_gr_response_active.py': {
            '--verify-navigation': '--verify-response-navigation',
            '--verify-closure': '--verify-geometry-closure',
            '--verify-previous-navigation': '--verify-previous-navigation',
            '--verify-prior-closure': '--verify-prior-closure'},
        'scripts/run_research_thermal_geometry_closed.py': {
            '--verify-navigation': '--verify-previous-navigation',
            '--verify-closure': '--verify-geometry-closure',
            '--verify-prior-closure': '--verify-prior-closure'},
        'scripts/run_research_thermal_geometry_active.py': {
            '--verify-navigation': '--verify-previous-navigation',
            '--verify-closure': '--verify-prior-closure'},
        'scripts/run_research_thermodynamics_closed.py': {
            '--verify-navigation': '--verify-gr-navigation',
            '--verify-closure': '--verify-prior-closure'},
        'scripts/run_research_gr_active.py': {
            '--verify-navigation': '--verify-gr-navigation',
            '--verify-closure': '--verify-local-closure'},
    }

    def mapped(command, *args, **kwargs):
        if isinstance(command, (list, tuple)) and len(command) >= 6:
            same_python = os.path.normcase(os.path.abspath(str(command[0]))) == os.path.normcase(sys.executable)
            target = relative(command[4])
            if same_python and list(command[1:4]) == ['-B', '-X', 'utf8'] and target in modes:
                flag = command[5]
                if flag in ('--script', '--verify-layout') or flag in modes[target]:
                    command = [*command[:4], str(THIS), modes[target].get(flag, flag), *command[6:]]
                    if kwargs.get('cwd') is not None:
                        name = relative(kwargs['cwd'])
                        if inside(name, previous.OLD):
                            kwargs['cwd'] = str(ROOT / (previous.NEW + name[len(previous.OLD):]))
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
    group.add_argument('--verify-navigation', action='store_true')
    group.add_argument('--verify-closure', action='store_true')
    group.add_argument('--verify-layout', action='store_true')
    for flag in ('--verify-response-navigation', '--verify-geometry-closure',
                 '--verify-previous-navigation', '--verify-prior-closure',
                 '--verify-gr-navigation', '--verify-local-closure'):
        group.add_argument(flag, action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = MatterGaugeLayout()
    runtime = MatterGaugeRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    previous.routes.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), previous.routes.execution_routes(previous.ExecutionLayout(store), runtime), frozen_runner_routes():
        if args.script:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)
        else:
            result = (store.verify_response_navigation() if args.verify_response_navigation else
                      previous.geometry.ClosedThermalGeometryLayout.verify_closure(store) if args.verify_geometry_closure else
                      store.verify_active_navigation() if args.verify_previous_navigation else
                      previous.thermal.ThermalLayout.verify_closure(store) if args.verify_prior_closure else
                      store.verify_navigation() if args.verify_gr_navigation else
                      store.verify_local_closure() if args.verify_local_closure else
                      store.verify_layout() if args.verify_layout else
                      store.verify_closure() if args.verify_closure else store.verify_matter_navigation())
            print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
