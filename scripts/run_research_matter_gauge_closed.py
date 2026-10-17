"""Read-only replay of the closed 1281--1287 matter/gauge stage.

Current closure bytes and inventory are checked before Markdown edits are
reversed in memory. Old/new paths, scoped imports and child executions use the
same historical view. No digest replacement, shadow tree, ZIP or symlink is used.
"""
from __future__ import annotations

import argparse
import builtins
from contextlib import contextmanager
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import types

import run_research_matter_gauge_active as previous
from research_layout import ROOT, Layout, ResearchLoader, _glob_match
from research_layout import no_archive_or_research_writes, relocation_original_bytes

response = previous.previous
geometry, thermal, gr, routes = response.geometry, response.thermal, response.gr, response.routes
THIS = Path(__file__).resolve()
BASE = 'research_cognition_physics'
OLD, NEW = BASE + '/archive_1281_', BASE + '/archive_1281_1287'
MANIFEST = ROOT / BASE / '_migration/closure_1281_1287/manifest.json'
VIEW = 'matter_gauge_closed'
digest, relative, inside = previous.digest, previous.relative, previous.inside
safe_relative = response.safe_relative


class ClosedMatterGaugeLayout(previous.MatterGaugeLayout):
    def __init__(self):
        raw = MANIFEST.read_bytes()
        self.matter_closure_sha256 = digest(raw)
        self.matter_closure_manifest = json.loads(raw)
        if self.matter_closure_manifest.get('schema') != 'research_closure_1281_1287_v1':
            raise ValueError('Unexpected matter/gauge closure schema')
        if self.matter_closure_manifest.get('stages') != [dict(
                source=OLD, destination=NEW, first_round=1281, last_round=1287)]:
            raise ValueError('Unexpected matter/gauge stage boundaries')
        self.matter_closure_entries, self.matter_closure_reverse = {}, {}
        old_names, new_names = set(), set()
        for entry in self.matter_closure_manifest['entries']:
            old, new = gr.validate_entry(entry)
            expected = NEW + old[len(OLD):] if inside(old, OLD) else old
            if new != expected or (not inside(old, OLD) and not old.endswith('.md')):
                raise ValueError('Unexpected matter/gauge closure destination: ' + old)
            if os.path.normcase(old) in old_names or os.path.normcase(new) in new_names:
                raise ValueError('Duplicate matter/gauge closure mapping')
            old_names.add(os.path.normcase(old))
            new_names.add(os.path.normcase(new))
            self.matter_closure_entries[old] = entry
            self.matter_closure_reverse[new] = old
        if not self.matter_closure_entries:
            raise ValueError('Empty matter/gauge closure inventory')
        self.matter_closure_directories = {
            p.as_posix() for name in self.matter_closure_entries
            for p in PurePosixPath(name).parents}
        super().__init__()
        self.current_files.update(self.matter_closure_entries)
        self.current_directories.update(self.matter_closure_directories)
        for name, boundary in self.matter_entries.items():
            entry = self.matter_closure_entries.get(name)
            if entry and (entry['original_sha256'], entry['original_bytes']) != (
                    boundary['current_sha256'], boundary['current_bytes']):
                raise ValueError('Closure does not join active matter navigation: ' + name)

    def matter_closure_name(self, path):
        name = relative(path)
        if name in self.matter_closure_reverse:
            return self.matter_closure_reverse[name]
        return OLD + name[len(NEW):] if inside(name, NEW) else name

    def matter_closure_logical(self, path):
        name = self.matter_closure_name(path)
        return ROOT / name if name is not None else Path(str(path))

    def matter_closure_actual(self, path):
        name = self.matter_closure_name(path)
        if name in self.matter_closure_entries:
            return ROOT / self.matter_closure_entries[name]['destination']
        return ROOT / (NEW + name[len(OLD):]) if inside(name, OLD) else Path(str(path))

    def before_matter_closure(self, path):
        name = self.matter_closure_name(path)
        if inside(name, OLD) and name not in self.matter_closure_entries:
            raise FileNotFoundError('File absent from matter closure: ' + str(name))
        raw = self.matter_closure_actual(path).read_bytes()
        entry = self.matter_closure_entries.get(name)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def before_matter_navigation(self, path):
        logical = self.matter_closure_logical(path)
        name = self.response_closure_name(logical)
        raw = self.before_matter_closure(self.response_closure_actual(logical))
        entry = self.matter_entries.get(name)
        return relocation_original_bytes(raw, entry) if entry else raw

    def key(self, path, view='legacy'):
        logical = self.matter_closure_logical(path)
        return relative(logical) if view == VIEW else super().key(logical, view)

    def actual(self, path, view='legacy'):
        if view == VIEW:
            return self.matter_closure_actual(path)
        return self.matter_closure_actual(super().actual(self.matter_closure_logical(path), view))

    def read(self, path, historical=True, view='legacy'):
        if view == VIEW:
            return (self.before_matter_closure(path) if historical
                    else self.matter_closure_actual(path).read_bytes())
        return super().read(self.matter_closure_logical(path), historical=historical, view=view)

    def glob_paths(self, directory, pattern, view):
        name = self.matter_closure_name(directory)
        if view != VIEW and not inside(name, OLD):
            yield from super().glob_paths(self.matter_closure_logical(directory), pattern, view)
            return
        prefix = name.rstrip('/') + '/' if name not in (None, '.') else ''
        found = {candidate: ROOT / candidate
                 for candidate in set(self.matter_closure_entries) | self.matter_closure_directories
                 if candidate.startswith(prefix) and _glob_match(candidate[len(prefix):], pattern)}
        if not inside(name, OLD):
            native = self.matter_closure_actual(directory)
            if native.is_dir():
                for actual in native.glob(pattern):
                    logical = self.matter_closure_name(actual)
                    if logical and logical.startswith(prefix) and _glob_match(logical[len(prefix):], pattern):
                        found[logical] = ROOT / logical
                    elif name is None:
                        found[str(actual)] = actual
        yield from (found[key] for key in sorted(found))

    def path_class(self, view='legacy'):
        if view != VIEW:
            return super().path_class(view)
        store = self
        base = Layout.path_class(self, view)

        class ClosedMatterPath(base):
            def exists(self, **kwargs):
                return store.key(self, VIEW) in store.matter_closure_directories or store.actual(self, VIEW).exists()

            def is_dir(self, **kwargs):
                return store.key(self, VIEW) in store.matter_closure_directories or store.actual(self, VIEW).is_dir()

        return ClosedMatterPath

    def script_key(self, path):
        name = self.matter_closure_name(path)
        if inside(name, OLD):
            return name if name in self.matter_closure_entries and name.endswith('.py') else None
        return super().script_key(self.matter_closure_logical(path))

    def script_actual(self, key):
        if inside(key, OLD):
            return self.matter_closure_actual(ROOT / key)
        return self.matter_closure_actual(super().script_actual(key))

    def verify_matter_navigation(self):
        if digest(previous.LAYER.read_bytes()) != self.matter_sha256:
            raise ValueError('Matter/gauge active layer changed during replay')
        for name, entry in self.matter_entries.items():
            old = self.before_matter_navigation(ROOT / name)
            enc = 'utf-8-sig' if old.startswith(b'\xef\xbb\xbf') else 'utf8'
            if entry['edits'][0]['start'] != old.decode(enc).index('\n') + 1:
                raise ValueError('Navigation insertion must follow the original title: ' + name)
        support = self.matter_layer['unchanged_support_sha256']
        required = dict(self.response_closure_manifest['unchanged_support_sha256'])
        required[response.MANIFEST.relative_to(ROOT).as_posix()] = self.response_closure_sha256
        required[response.THIS.relative_to(ROOT).as_posix()] = digest(response.THIS.read_bytes())
        if not isinstance(support, dict) or not set(required).issubset(support):
            raise ValueError('Missing inherited support or closure bindings')
        for name, expected in support.items():
            if digest(self.before_matter_closure(ROOT / safe_relative(name))) != expected:
                raise ValueError('Frozen navigation support changed: ' + name)
            if name in required and expected != required[name]:
                raise ValueError('Inherited support binding changed: ' + name)
        for name in previous.ENTRY_POINTS:
            self.before_matter_closure(ROOT / name)
        return dict(status='PASS', files=4, manifest_sha256=self.matter_sha256,
                    closed_boundary_verified=True, original_navigation_recoverable=True,
                    frozen_archives_modified=False, stage_first_round=1281,
                    entry_points=previous.ENTRY_POINTS, frozen_support_files=len(support),
                    scientific_acceptance_checked=False)

    def verify_closure(self):
        if digest(MANIFEST.read_bytes()) != self.matter_closure_sha256:
            raise ValueError('Matter/gauge closure changed during replay')
        expected, scientific, markdown = set(), 0, 0
        for name, entry in self.matter_closure_entries.items():
            raw = self.before_matter_closure(ROOT / name)
            if inside(name, OLD):
                expected.add(entry['destination'])
            scientific += Path(name).suffix.lower() in ('.py', '.json')
            markdown += bool(entry['edits'])
            if Path(name).suffix.lower() in ('.py', '.json') and digest(raw) != entry['current_sha256']:
                raise ValueError('Scientific evidence changed: ' + name)
        actual = {p.relative_to(ROOT).as_posix() for p in (ROOT / NEW).rglob('*') if p.is_file()}
        if (ROOT / OLD).exists() or not (ROOT / NEW).is_dir() or actual != expected:
            raise ValueError('Closed matter/gauge inventory differs from manifest')
        support = self.matter_closure_manifest['unchanged_support_sha256']
        if not isinstance(support, dict) or not support:
            raise ValueError('Frozen support bindings are required')
        for name, expected_hash in support.items():
            if digest((ROOT / safe_relative(name)).read_bytes()) != expected_hash:
                raise ValueError('Frozen support changed: ' + name)
        return dict(status='PASS', scope='1281--1287 directory closure only',
                    files=len(self.matter_closure_entries), moved_files=len(expected),
                    scientific_python_json_byte_identical=scientific,
                    reversible_markdown_files=markdown, frozen_support_files=len(support),
                    manifest_sha256=self.matter_closure_sha256,
                    old_directory_absent=True, current_and_original_hashes_verified=True,
                    directory_copy_used=False, zip_used=False,
                    prior_closure_external_edit_not_resolved=True)

    def verify_layout(self):
        return dict(matter_closure=self.verify_closure(),
                    recovered_matter_navigation=self.verify_matter_navigation(),
                    response_closure=response.ClosedGRResponseLayout.verify_closure(self),
                    recovered_response_navigation=self.verify_response_navigation(),
                    thermal_geometry_closure=geometry.ClosedThermalGeometryLayout.verify_closure(self),
                    recovered_previous_navigation=self.verify_active_navigation(),
                    prior_thermal_closure=thermal.ThermalLayout.verify_closure(self),
                    frozen_layers=gr.GRLayout.verify_layout(self), all_checks_passed=True)


class ClosedMatterGaugeLoader(ResearchLoader):
    def exec_module(self, module):
        if not inside(self.key, OLD):
            return previous.MatterGaugeLoader(self.runtime, self.key).exec_module(module)
        store = self.runtime.store
        logical, actual = ROOT / self.key, store.script_actual(self.key)
        module.__file__, module.__builtins__ = str(logical), self.runtime.scoped_builtins[VIEW]
        old_path = list(sys.path)
        sys.path.insert(0, str(logical.parent))
        try:
            exec(compile(store.read(logical, view=VIEW), str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = old_path


def adapt_reader_imports(runtime, layout_class, runtime_class, extra=()):
    """Replace only reader classes imported by replayed modules, not real modules."""
    proxies = {}
    readers = (
        (previous, 'MatterGaugeLayout', 'MatterGaugeRuntime'),
        (response, 'ClosedGRResponseLayout', 'ClosedGRResponseRuntime'),
        (response.previous, 'GRResponseLayout', 'GRResponseRuntime'),
        (geometry.previous, 'ThermalGeometryLayout', 'ThermalGeometryRuntime'),
        *extra)
    for module, layout_name, runtime_name in readers:
        proxy = types.ModuleType(module.__name__)
        proxy.__dict__.update(vars(module))
        setattr(proxy, layout_name, layout_class)
        setattr(proxy, runtime_name, runtime_class)
        proxies[module.__name__] = proxy
    for bindings in runtime.scoped_builtins.values():
        fallback = bindings['__import__']
        old_library = fallback('importlib')
        library = types.ModuleType('importlib')
        library.__dict__.update(vars(old_library))

        def import_module(name, package=None, original=old_library.import_module):
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


class ClosedMatterGaugeRuntime(previous.MatterGaugeRuntime):
    def __init__(self, store):
        super().__init__(store)
        import pathlib
        proxy = types.ModuleType('pathlib')
        proxy.__dict__.update(vars(pathlib))
        proxy.Path = store.path_class(VIEW)
        default_import, default_open = builtins.__import__, builtins.open

        def scoped_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == 'pathlib' and level == 0:
                return proxy
            return default_import(name, globals, locals, fromlist, level)

        def scoped_open(file, mode='r', buffering=-1, encoding=None, errors=None,
                        newline=None, closefd=True, opener=None):
            if isinstance(file, (str, bytes, os.PathLike)):
                name = os.fsdecode(file)
                key = store.key(name, VIEW)
                if key and (key.startswith(BASE + '/') or key in store.matter_closure_entries):
                    if opener is not None or not closefd:
                        raise NotImplementedError('Custom historical openers are unsupported')
                    return proxy.Path(name).open(mode, buffering, encoding, errors, newline)
            return default_open(file, mode, buffering, encoding, errors, newline, closefd, opener)

        self.scoped_builtins[VIEW] = dict(vars(builtins), __import__=scoped_import, open=scoped_open)
        adapt_reader_imports(self, ClosedMatterGaugeLayout, ClosedMatterGaugeRuntime)

    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, ClosedMatterGaugeLoader(self, key),
                                              origin=str(self.store.script_actual(key)))


class ExecutionLayout(response.ExecutionLayout):
    def actual(self, path, view='legacy'):
        return self.store.matter_closure_actual(super().actual(path, view))


RUNNER_MODES = {
    'scripts/run_research_matter_gauge_active.py': {
        '--verify-navigation': '--verify-navigation', '--verify-closure': '--verify-response-closure',
        '--verify-response-navigation': '--verify-response-navigation',
        '--verify-geometry-closure': '--verify-geometry-closure',
        '--verify-previous-navigation': '--verify-previous-navigation',
        '--verify-prior-closure': '--verify-prior-closure',
        '--verify-gr-navigation': '--verify-gr-navigation', '--verify-local-closure': '--verify-local-closure'},
    'scripts/run_research_gr_response_closed.py': {
        '--verify-navigation': '--verify-response-navigation', '--verify-closure': '--verify-response-closure',
        '--verify-geometry-closure': '--verify-geometry-closure',
        '--verify-previous-navigation': '--verify-previous-navigation',
        '--verify-prior-closure': '--verify-prior-closure',
        '--verify-gr-navigation': '--verify-gr-navigation', '--verify-local-closure': '--verify-local-closure'},
    'scripts/run_research_gr_response_active.py': {
        '--verify-navigation': '--verify-response-navigation', '--verify-closure': '--verify-geometry-closure',
        '--verify-previous-navigation': '--verify-previous-navigation', '--verify-prior-closure': '--verify-prior-closure'},
    'scripts/run_research_thermal_geometry_closed.py': {
        '--verify-navigation': '--verify-previous-navigation', '--verify-closure': '--verify-geometry-closure',
        '--verify-prior-closure': '--verify-prior-closure'},
    'scripts/run_research_thermal_geometry_active.py': {
        '--verify-navigation': '--verify-previous-navigation', '--verify-closure': '--verify-prior-closure'},
    'scripts/run_research_thermodynamics_closed.py': {
        '--verify-navigation': '--verify-gr-navigation', '--verify-closure': '--verify-prior-closure'},
    'scripts/run_research_gr_active.py': {
        '--verify-navigation': '--verify-gr-navigation', '--verify-closure': '--verify-local-closure'},
}
HISTORICAL_FLAGS = ('--verify-response-closure', '--verify-response-navigation', '--verify-geometry-closure',
                    '--verify-previous-navigation', '--verify-prior-closure', '--verify-gr-navigation',
                    '--verify-local-closure')


@contextmanager
def frozen_runner_routes(target=THIS, modes=None):
    old_run = subprocess.run
    modes = RUNNER_MODES if modes is None else modes

    def mapped(command, *args, **kwargs):
        if isinstance(command, (list, tuple)) and len(command) >= 6:
            same_python = os.path.normcase(os.path.abspath(str(command[0]))) == os.path.normcase(sys.executable)
            name = relative(command[4])
            if same_python and list(command[1:4]) == ['-B', '-X', 'utf8'] and name in modes:
                flag = command[5]
                if flag in ('--script', '--verify-layout') or flag in modes[name]:
                    command = [*command[:4], str(target), modes[name].get(flag, flag), *command[6:]]
                    if kwargs.get('cwd') is not None:
                        cwd = relative(kwargs['cwd'])
                        for old, new in ((OLD, NEW), (response.OLD, response.NEW)):
                            if inside(cwd, old):
                                kwargs['cwd'] = str(ROOT / (new + cwd[len(old):]))
                                break
        return old_run(command, *args, **kwargs)

    subprocess.run = mapped
    try:
        yield
    finally:
        subprocess.run = old_run


def historical_check(store, flag):
    return {
        '--verify-response-closure': lambda: response.ClosedGRResponseLayout.verify_closure(store),
        '--verify-response-navigation': store.verify_response_navigation,
        '--verify-geometry-closure': lambda: geometry.ClosedThermalGeometryLayout.verify_closure(store),
        '--verify-previous-navigation': store.verify_active_navigation,
        '--verify-prior-closure': lambda: thermal.ThermalLayout.verify_closure(store),
        '--verify-gr-navigation': store.verify_navigation,
        '--verify-local-closure': store.verify_local_closure,
    }[flag]()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--script')
    group.add_argument('--verify-closure', action='store_true')
    group.add_argument('--verify-navigation', action='store_true')
    group.add_argument('--verify-layout', action='store_true')
    for flag in HISTORICAL_FLAGS:
        group.add_argument(flag, action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = ClosedMatterGaugeLayout()
    runtime = ClosedMatterGaugeRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    routes.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), routes.execution_routes(ExecutionLayout(store), runtime), frozen_runner_routes():
        if args.script:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)
        else:
            selected = next((flag for flag in HISTORICAL_FLAGS if getattr(args, flag[2:].replace('-', '_'))), None)
            result = (historical_check(store, selected) if selected else
                      store.verify_layout() if args.verify_layout else
                      store.verify_closure() if args.verify_closure else store.verify_matter_navigation())
            print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
