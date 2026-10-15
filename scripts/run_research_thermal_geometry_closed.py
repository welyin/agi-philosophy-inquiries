"""Read-only replay of the closed 1219--1233 thermal--geometry stage.

Only recorded Markdown changes are reversed, in memory, after checking the
actual destination bytes. Old and new script paths share the original logical
view. No directory copy, ZIP, global pathlib replacement, or digest override
is used. --verify-closure verifies this migration, not the earlier closure
whose registered external edit to paper 08 remains unresolved.
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

import run_research_thermal_geometry_active as previous
from research_layout import ROOT, Layout, ResearchLoader, _glob_match
from research_layout import no_archive_or_research_writes, relocation_original_bytes

THIS = Path(__file__).resolve()
BASE = 'research_cognition_physics'
OLD = BASE + '/archive_1219_'
NEW = BASE + '/archive_1219_1233'
MANIFEST = ROOT / BASE / '_migration/closure_1219_1233/manifest.json'
VIEW = 'thermal_geometry_closed'
digest = previous.previous.digest
relative = previous.previous.relative
safe_relative = previous.previous.safe_relative
inside = previous.previous.inside


class ClosedThermalGeometryLayout(previous.ThermalGeometryLayout):
    def __init__(self):
        raw = MANIFEST.read_bytes()
        self.geometry_closure_sha256 = digest(raw)
        self.geometry_closure_manifest = json.loads(raw)
        if self.geometry_closure_manifest.get('schema') != 'research_closure_1219_1233_v1':
            raise ValueError('Unexpected thermal-geometry closure schema')
        if self.geometry_closure_manifest.get('stages') != [dict(
                source=OLD, destination=NEW, first_round=1219, last_round=1233)]:
            raise ValueError('Unexpected thermal-geometry stage boundaries')
        self.geometry_closure_entries, self.geometry_closure_reverse = {}, {}
        old_names, new_names = set(), set()
        for entry in self.geometry_closure_manifest['entries']:
            old, new = previous.previous.previous.validate_entry(entry)
            expected = NEW + old[len(OLD):] if inside(old, OLD) else old
            if new != expected or (not inside(old, OLD) and not old.endswith('.md')):
                raise ValueError('Unexpected thermal-geometry closure destination: ' + old)
            if os.path.normcase(old) in old_names or os.path.normcase(new) in new_names:
                raise ValueError('Duplicate thermal-geometry closure mapping')
            old_names.add(os.path.normcase(old)); new_names.add(os.path.normcase(new))
            self.geometry_closure_entries[old] = entry
            self.geometry_closure_reverse[new] = old
        self.geometry_closure_directories = {
            p.as_posix() for name in self.geometry_closure_entries for p in PurePosixPath(name).parents}
        super().__init__()
        self.current_files.update(self.geometry_closure_entries)
        self.current_directories.update(self.geometry_closure_directories)

    def geometry_closure_name(self, path):
        name = relative(path)
        if name in self.geometry_closure_reverse:
            return self.geometry_closure_reverse[name]
        return OLD + name[len(NEW):] if inside(name, NEW) else name

    def geometry_closure_logical(self, path):
        name = self.geometry_closure_name(path)
        return ROOT / name if name is not None else Path(str(path))

    def geometry_closure_actual(self, path):
        name = self.geometry_closure_name(path)
        if name in self.geometry_closure_entries:
            return ROOT / self.geometry_closure_entries[name]['destination']
        return ROOT / (NEW + name[len(OLD):]) if inside(name, OLD) else Path(str(path))

    def before_geometry_closure(self, path):
        name = self.geometry_closure_name(path)
        if inside(name, OLD) and name not in self.geometry_closure_entries:
            raise FileNotFoundError('File absent from thermal-geometry closure: ' + str(name))
        raw = self.geometry_closure_actual(path).read_bytes()
        entry = self.geometry_closure_entries.get(name)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def before_thermal(self, path):
        name = self.thermal_name(self.geometry_closure_logical(path))
        if previous.previous.moved(name) and name not in self.thermal_entries:
            raise FileNotFoundError('File absent from thermal closure inventory: ' + str(name))
        raw = self.before_geometry_closure(self.thermal_actual(path))
        active = self.geometry_entries.get(name)
        if active:
            raw = relocation_original_bytes(raw, active)
        entry = self.thermal_entries.get(name)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def key(self, path, view='legacy'):
        logical = self.geometry_closure_logical(path)
        return relative(logical) if view == VIEW else super().key(logical, view)

    def actual(self, path, view='legacy'):
        if view == VIEW:
            return self.geometry_closure_actual(path)
        return self.geometry_closure_actual(super().actual(self.geometry_closure_logical(path), view))

    def read(self, path, historical=True, view='legacy'):
        if view == VIEW:
            return self.before_geometry_closure(path) if historical else self.geometry_closure_actual(path).read_bytes()
        return super().read(self.geometry_closure_logical(path), historical=historical, view=view)

    def glob_paths(self, directory, pattern, view):
        name = self.geometry_closure_name(directory)
        if view != VIEW and not inside(name, OLD):
            yield from super().glob_paths(self.geometry_closure_logical(directory), pattern, view)
            return
        prefix = name.rstrip('/') + '/' if name not in (None, '.') else ''
        found = {candidate: ROOT / candidate
                 for candidate in set(self.geometry_closure_entries) | self.geometry_closure_directories
                 if candidate.startswith(prefix) and _glob_match(candidate[len(prefix):], pattern)}
        if not inside(name, OLD):
            native = self.geometry_closure_actual(directory)
            if native.is_dir():
                for actual in native.glob(pattern):
                    logical = self.geometry_closure_name(actual)
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

        class ClosedGeometryPath(base):
            # This view virtualizes only the newest moved archive. Earlier
            # stage source directories really are absent at this boundary;
            # inheriting all legacy directory aliases would conceal that fact.
            def exists(self, **kwargs):
                key = store.key(self, VIEW)
                return key in store.geometry_closure_directories or store.actual(self, VIEW).exists()

            def is_dir(self, **kwargs):
                key = store.key(self, VIEW)
                return key in store.geometry_closure_directories or store.actual(self, VIEW).is_dir()

        return ClosedGeometryPath

    def script_key(self, path):
        name = self.geometry_closure_name(path)
        if inside(name, OLD) and name in self.geometry_closure_entries and name.endswith('.py'):
            return name
        return super().script_key(self.geometry_closure_logical(path))

    def script_actual(self, key):
        if inside(key, OLD):
            return self.geometry_closure_actual(ROOT / key)
        return self.geometry_closure_actual(super().script_actual(key))

    def verify_active_navigation(self):
        if digest(previous.LAYER.read_bytes()) != self.geometry_sha256:
            raise ValueError('Active layer changed during replay')
        for name, entry in self.geometry_entries.items():
            relocation_original_bytes(self.before_geometry_closure(ROOT / name), entry)
        for name, expected in self.geometry_layer['unchanged_support_sha256'].items():
            if digest(self.before_geometry_closure(ROOT / safe_relative(name))) != expected:
                raise ValueError('Frozen support changed: ' + name)
        # Preserve the original runner's receipt fields and exact meaning.
        return dict(status='PASS', files=4, manifest_sha256=self.geometry_sha256,
                    closed_boundary_verified=True, original_navigation_recoverable=True,
                    frozen_archives_modified=False)

    def verify_closure(self):
        if digest(MANIFEST.read_bytes()) != self.geometry_closure_sha256:
            raise ValueError('Thermal-geometry closure changed during replay')
        expected, scientific, markdown = set(), 0, 0
        for name, entry in self.geometry_closure_entries.items():
            raw = self.before_geometry_closure(ROOT / name)
            if inside(name, OLD):
                expected.add(entry['destination'])
            scientific += Path(name).suffix.lower() in ('.py', '.json')
            markdown += bool(entry['edits'])
            if Path(name).suffix.lower() in ('.py', '.json') and digest(raw) != entry['current_sha256']:
                raise ValueError('Scientific evidence changed: ' + name)
        actual = {p.relative_to(ROOT).as_posix() for p in (ROOT / NEW).rglob('*') if p.is_file()}
        if (ROOT / OLD).exists() or not (ROOT / NEW).is_dir() or actual != expected:
            raise ValueError('Closed thermal-geometry inventory differs from manifest')
        support = self.geometry_closure_manifest['unchanged_support_sha256']
        if not isinstance(support, dict) or not support:
            raise ValueError('Frozen support bindings are required')
        for name, expected_hash in support.items():
            if digest((ROOT / safe_relative(name)).read_bytes()) != expected_hash:
                raise ValueError('Frozen support changed: ' + name)
        return dict(status='PASS', scope='1219--1233 directory closure only',
                    files=len(self.geometry_closure_entries), moved_files=len(expected),
                    scientific_python_json_byte_identical=scientific,
                    reversible_markdown_files=markdown, frozen_support_files=len(support),
                    manifest_sha256=self.geometry_closure_sha256,
                    old_directory_absent=True, current_and_original_hashes_verified=True,
                    directory_copy_used=False, zip_used=False,
                    prior_closure_external_edit_not_resolved=True)

    def verify_layout(self):
        # This deliberately does not promote the known old paper-08 exception.
        return dict(thermal_geometry_closure=self.verify_closure(),
                    recovered_active_navigation=self.verify_active_navigation(),
                    prior_thermal_closure=previous.previous.ThermalLayout.verify_closure(self),
                    frozen_layers=super().verify_layout(), all_checks_passed=True)


class ClosedThermalGeometryLoader(ResearchLoader):
    def exec_module(self, module):
        if not inside(self.key, OLD):
            return previous.ThermalGeometryLoader(self.runtime, self.key).exec_module(module)
        store = self.runtime.store
        logical, actual = ROOT / self.key, store.script_actual(self.key)
        module.__file__, module.__builtins__ = str(logical), self.runtime.scoped_builtins[VIEW]
        old_path = list(sys.path)
        sys.path.insert(0, str(logical.parent))
        try:
            exec(compile(store.read(logical, view=VIEW), str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = old_path


class ClosedThermalGeometryRuntime(previous.ThermalGeometryRuntime):
    def __init__(self, store):
        super().__init__(store)
        import pathlib
        path_proxy = types.ModuleType('pathlib')
        path_proxy.__dict__.update(vars(pathlib))
        path_proxy.Path = store.path_class(VIEW)
        import_proxy = types.ModuleType('importlib')
        import_proxy.__dict__.update(vars(importlib))
        active_proxy = types.ModuleType(previous.__name__)
        active_proxy.__dict__.update(vars(previous))
        active_proxy.ThermalGeometryLayout = ClosedThermalGeometryLayout
        default_import, default_open = builtins.__import__, builtins.open
        default_import_module = importlib.import_module

        def historical_import_module(name, package=None):
            # The frozen publisher creates its own old active-layout instance.
            # Give just that scoped import the outer reader; its old constants
            # and manifests still identify the old thermal closure and hashes.
            return active_proxy if name == previous.__name__ and package is None else default_import_module(name, package)

        import_proxy.import_module = historical_import_module

        def scoped_import(name, globals=None, locals=None, fromlist=(), level=0):
            if level == 0 and name == 'pathlib':
                return path_proxy
            if level == 0 and name == 'importlib':
                return import_proxy
            return default_import(name, globals, locals, fromlist, level)

        def scoped_open(file, mode='r', buffering=-1, encoding=None, errors=None,
                        newline=None, closefd=True, opener=None):
            if isinstance(file, (str, bytes, os.PathLike)):
                name = os.fsdecode(file)
                key = store.key(name, VIEW)
                if key and (key.startswith(BASE + '/') or key in store.geometry_closure_entries):
                    if opener is not None or not closefd:
                        raise NotImplementedError('Custom historical openers are unsupported')
                    return path_proxy.Path(name).open(mode, buffering, encoding, errors, newline)
            return default_open(file, mode, buffering, encoding, errors, newline, closefd, opener)
        self.scoped_builtins[VIEW] = dict(vars(builtins), __import__=scoped_import, open=scoped_open)

    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, ClosedThermalGeometryLoader(self, key),
                                              origin=str(self.store.script_actual(key)))


class ExecutionLayout(previous.previous.ExecutionLayout):
    def actual(self, path, view='legacy'):
        return self.store.geometry_closure_actual(super().actual(path, view))


@contextmanager
def frozen_runner_routes():
    old_run = subprocess.run

    def mapped(command, *args, **kwargs):
        if isinstance(command, (list, tuple)) and len(command) >= 6:
            same_python = os.path.normcase(os.path.abspath(str(command[0]))) == os.path.normcase(sys.executable)
            target = relative(command[4])
            if same_python and list(command[1:4]) == ['-B', '-X', 'utf8']:
                if target in ('scripts/run_research_thermal_geometry_active.py',
                              'scripts/run_research_thermodynamics_closed.py'):
                    flag = command[5]
                    if flag in ('--verify-closure', '--verify-navigation', '--verify-layout', '--script'):
                        # An old closure check must retain its old scope. In
                        # particular, paper 08 must still fail an exact replay.
                        flag = '--verify-prior-closure' if flag == '--verify-closure' else flag
                        command = [*command[:4], str(THIS), flag, *command[6:]]
        return old_run(command, *args, **kwargs)

    subprocess.run = mapped
    try:
        yield
    finally:
        subprocess.run = old_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--script')
    choice.add_argument('--verify-layout', action='store_true')
    choice.add_argument('--verify-closure', action='store_true')
    choice.add_argument('--verify-navigation', action='store_true')
    choice.add_argument('--verify-prior-closure', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = ClosedThermalGeometryLayout()
    runtime = ClosedThermalGeometryRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    routes = previous.previous.previous.combined.merged.previous.previous
    routes.THIS = THIS
    previous.previous.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), routes.execution_routes(ExecutionLayout(store), runtime), previous.previous.gr_verification_route(), frozen_runner_routes():
        if args.script:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)
        else:
            result = (previous.previous.ThermalLayout.verify_closure(store) if args.verify_prior_closure else
                      store.verify_layout() if args.verify_layout else
                      store.verify_closure() if args.verify_closure else store.verify_active_navigation())
            print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
