"""Read-only replay of the closed 1234--1280 GR-response stage.

The current closure is verified before recorded Markdown edits are reversed in
memory. Scientific Python/JSON remain byte-identical. Old and new script names
share the pre-closure logical view, including __file__, Path, open, imports and
child dispatch. No shadow checkout, archive extraction or digest override is
used. Earlier closure exceptions keep their original verification scope.
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

import run_research_gr_response_active as previous
from research_layout import ROOT, Layout, ResearchLoader, _glob_match
from research_layout import no_archive_or_research_writes, relocation_original_bytes

geometry = previous.previous
thermal = geometry.previous.previous
gr = thermal.previous
routes = gr.combined.merged.previous.previous
THIS = Path(__file__).resolve()
BASE = 'research_cognition_physics'
OLD = BASE + '/archive_1234_'
NEW = BASE + '/archive_1234_1280'
MANIFEST = ROOT / BASE / '_migration/closure_1234_1280/manifest.json'
VIEW = 'gr_response_closed'
digest, relative = geometry.digest, geometry.relative
safe_relative, inside = geometry.safe_relative, geometry.inside


class ClosedGRResponseLayout(previous.GRResponseLayout):
    def __init__(self):
        raw = MANIFEST.read_bytes()
        self.response_closure_sha256 = digest(raw)
        self.response_closure_manifest = json.loads(raw)
        if self.response_closure_manifest.get('schema') != 'research_closure_1234_1280_v1':
            raise ValueError('Unexpected GR-response closure schema')
        if self.response_closure_manifest.get('stages') != [dict(
                source=OLD, destination=NEW, first_round=1234, last_round=1280)]:
            raise ValueError('Unexpected GR-response stage boundaries')
        self.response_closure_entries, self.response_closure_reverse = {}, {}
        old_names, new_names = set(), set()
        for entry in self.response_closure_manifest['entries']:
            old, new = gr.validate_entry(entry)
            expected = NEW + old[len(OLD):] if inside(old, OLD) else old
            if new != expected or (not inside(old, OLD) and not old.endswith('.md')):
                raise ValueError('Unexpected GR-response closure destination: ' + old)
            if os.path.normcase(old) in old_names or os.path.normcase(new) in new_names:
                raise ValueError('Duplicate GR-response closure mapping')
            old_names.add(os.path.normcase(old))
            new_names.add(os.path.normcase(new))
            self.response_closure_entries[old] = entry
            self.response_closure_reverse[new] = old
        if not self.response_closure_entries:
            raise ValueError('Empty GR-response closure inventory')
        self.response_closure_directories = {
            p.as_posix() for name in self.response_closure_entries
            for p in PurePosixPath(name).parents}
        super().__init__()
        self.current_files.update(self.response_closure_entries)
        self.current_directories.update(self.response_closure_directories)
        # The old active layer stays frozen. Join it to this outer layer rather
        # than rebinding its hashes to the final navigation.
        for name, boundary in self.response_entries.items():
            entry = self.response_closure_entries.get(name)
            if entry and (entry['original_sha256'], entry['original_bytes']) != (
                    boundary['current_sha256'], boundary['current_bytes']):
                raise ValueError('Closure does not join active GR-response navigation: ' + name)

    def response_closure_name(self, path):
        name = relative(path)
        if name in self.response_closure_reverse:
            return self.response_closure_reverse[name]
        return OLD + name[len(NEW):] if inside(name, NEW) else name

    def response_closure_logical(self, path):
        name = self.response_closure_name(path)
        return ROOT / name if name is not None else Path(str(path))

    def response_closure_actual(self, path):
        name = self.response_closure_name(path)
        if name in self.response_closure_entries:
            return ROOT / self.response_closure_entries[name]['destination']
        return ROOT / (NEW + name[len(OLD):]) if inside(name, OLD) else Path(str(path))

    def before_response_closure(self, path):
        name = self.response_closure_name(path)
        if inside(name, OLD) and name not in self.response_closure_entries:
            raise FileNotFoundError('File absent from GR-response closure: ' + str(name))
        raw = self.response_closure_actual(path).read_bytes()
        entry = self.response_closure_entries.get(name)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def before_geometry_closure(self, path):
        logical = self.response_closure_logical(path)
        name = self.geometry_closure_name(logical)
        if inside(name, geometry.OLD) and name not in self.geometry_closure_entries:
            raise FileNotFoundError('File absent from thermal-geometry closure: ' + str(name))
        raw = self.before_response_closure(self.geometry_closure_actual(logical))
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

    def key(self, path, view='legacy'):
        logical = self.response_closure_logical(path)
        return relative(logical) if view == VIEW else super().key(logical, view)

    def actual(self, path, view='legacy'):
        if view == VIEW:
            return self.response_closure_actual(path)
        return self.response_closure_actual(super().actual(self.response_closure_logical(path), view))

    def read(self, path, historical=True, view='legacy'):
        if view == VIEW:
            return (self.before_response_closure(path) if historical
                    else self.response_closure_actual(path).read_bytes())
        return super().read(self.response_closure_logical(path), historical=historical, view=view)

    def glob_paths(self, directory, pattern, view):
        name = self.response_closure_name(directory)
        if view != VIEW and not inside(name, OLD):
            yield from super().glob_paths(self.response_closure_logical(directory), pattern, view)
            return
        prefix = name.rstrip('/') + '/' if name not in (None, '.') else ''
        found = {candidate: ROOT / candidate
                 for candidate in set(self.response_closure_entries) | self.response_closure_directories
                 if candidate.startswith(prefix) and _glob_match(candidate[len(prefix):], pattern)}
        if not inside(name, OLD):
            native = self.response_closure_actual(directory)
            if native.is_dir():
                for actual in native.glob(pattern):
                    logical = self.response_closure_name(actual)
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

        class ClosedResponsePath(base):
            def exists(self, **kwargs):
                key = store.key(self, VIEW)
                return key in store.response_closure_directories or store.actual(self, VIEW).exists()

            def is_dir(self, **kwargs):
                key = store.key(self, VIEW)
                return key in store.response_closure_directories or store.actual(self, VIEW).is_dir()

        return ClosedResponsePath

    def script_key(self, path):
        name = self.response_closure_name(path)
        if inside(name, OLD):
            return name if name in self.response_closure_entries and name.endswith('.py') else None
        return super().script_key(self.response_closure_logical(path))

    def script_actual(self, key):
        if inside(key, OLD):
            return self.response_closure_actual(ROOT / key)
        return self.response_closure_actual(super().script_actual(key))

    def verify_response_navigation(self):
        if digest(previous.LAYER.read_bytes()) != self.response_sha256:
            raise ValueError('GR-response active layer changed during replay')
        for name, entry in self.response_entries.items():
            relocation_original_bytes(self.before_response_closure(ROOT / name), entry)
        for name, expected in self.response_layer['unchanged_support_sha256'].items():
            if digest(self.before_response_closure(ROOT / safe_relative(name))) != expected:
                raise ValueError('Frozen navigation support changed: ' + name)
        notice = ROOT / safe_relative(self.response_layer['notice_path'])
        if digest(self.before_response_closure(notice)) != self.response_layer['notice_sha256']:
            raise ValueError('Active notice changed without a reversible navigation update')
        return dict(status='PASS', files=4, manifest_sha256=self.response_sha256,
                    closed_boundary_verified=True, original_navigation_recoverable=True,
                    frozen_archives_modified=False, stage_first_round=1234)

    def verify_closure(self):
        if digest(MANIFEST.read_bytes()) != self.response_closure_sha256:
            raise ValueError('GR-response closure changed during replay')
        expected, scientific, markdown = set(), 0, 0
        for name, entry in self.response_closure_entries.items():
            raw = self.before_response_closure(ROOT / name)
            if inside(name, OLD):
                expected.add(entry['destination'])
            scientific += Path(name).suffix.lower() in ('.py', '.json')
            markdown += bool(entry['edits'])
            if Path(name).suffix.lower() in ('.py', '.json') and digest(raw) != entry['current_sha256']:
                raise ValueError('Scientific evidence changed: ' + name)
        actual = {p.relative_to(ROOT).as_posix() for p in (ROOT / NEW).rglob('*') if p.is_file()}
        if (ROOT / OLD).exists() or not (ROOT / NEW).is_dir() or actual != expected:
            raise ValueError('Closed GR-response inventory differs from manifest')
        support = self.response_closure_manifest['unchanged_support_sha256']
        if not isinstance(support, dict) or not support:
            raise ValueError('Frozen support bindings are required')
        for name, expected_hash in support.items():
            if digest((ROOT / safe_relative(name)).read_bytes()) != expected_hash:
                raise ValueError('Frozen support changed: ' + name)
        return dict(status='PASS', scope='1234--1280 directory closure only',
                    files=len(self.response_closure_entries), moved_files=len(expected),
                    scientific_python_json_byte_identical=scientific,
                    reversible_markdown_files=markdown, frozen_support_files=len(support),
                    manifest_sha256=self.response_closure_sha256,
                    old_directory_absent=True, current_and_original_hashes_verified=True,
                    directory_copy_used=False, zip_used=False,
                    prior_closure_external_edit_not_resolved=True)

    def verify_layout(self):
        # Call each closure explicitly: a virtual self.verify_closure() in a
        # superclass would accidentally report the newest scope twice.
        return dict(response_closure=self.verify_closure(),
                    recovered_response_navigation=self.verify_response_navigation(),
                    thermal_geometry_closure=geometry.ClosedThermalGeometryLayout.verify_closure(self),
                    recovered_previous_navigation=self.verify_active_navigation(),
                    prior_thermal_closure=thermal.ThermalLayout.verify_closure(self),
                    frozen_layers=gr.GRLayout.verify_layout(self), all_checks_passed=True)


class ClosedGRResponseLoader(ResearchLoader):
    def exec_module(self, module):
        if not inside(self.key, OLD):
            return previous.GRResponseLoader(self.runtime, self.key).exec_module(module)
        store = self.runtime.store
        logical, actual = ROOT / self.key, store.script_actual(self.key)
        module.__file__, module.__builtins__ = str(logical), self.runtime.scoped_builtins[VIEW]
        old_path = list(sys.path)
        sys.path.insert(0, str(logical.parent))
        try:
            exec(compile(store.read(logical, view=VIEW), str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = old_path


class ClosedGRResponseRuntime(previous.GRResponseRuntime):
    def __init__(self, store):
        super().__init__(store)
        import pathlib
        path_proxy = types.ModuleType('pathlib')
        path_proxy.__dict__.update(vars(pathlib))
        path_proxy.Path = store.path_class(VIEW)
        default_import, default_open = builtins.__import__, builtins.open

        def scoped_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == 'pathlib' and level == 0:
                return path_proxy
            return default_import(name, globals, locals, fromlist, level)

        def scoped_open(file, mode='r', buffering=-1, encoding=None, errors=None,
                        newline=None, closefd=True, opener=None):
            if isinstance(file, (str, bytes, os.PathLike)):
                name = os.fsdecode(file)
                key = store.key(name, VIEW)
                if key and (key.startswith(BASE + '/') or key in store.response_closure_entries):
                    if opener is not None or not closefd:
                        raise NotImplementedError('Custom historical openers are unsupported')
                    return path_proxy.Path(name).open(mode, buffering, encoding, errors, newline)
            return default_open(file, mode, buffering, encoding, errors, newline, closefd, opener)

        self.scoped_builtins[VIEW] = dict(vars(builtins), __import__=scoped_import, open=scoped_open)
        # Frozen 1238 imports the active reader directly; an earlier publisher
        # uses importlib.import_module. Adapt only these scoped module objects.
        proxies = {}
        for module, layout_name, runtime_name in (
                (previous, 'GRResponseLayout', 'GRResponseRuntime'),
                (geometry.previous, 'ThermalGeometryLayout', 'ThermalGeometryRuntime')):
            proxy = types.ModuleType(module.__name__)
            proxy.__dict__.update(vars(module))
            setattr(proxy, layout_name, ClosedGRResponseLayout)
            setattr(proxy, runtime_name, ClosedGRResponseRuntime)
            proxies[module.__name__] = proxy
        for bindings in self.scoped_builtins.values():
            original_import = bindings['__import__']
            original_importlib = original_import('importlib')
            import_proxy = types.ModuleType('importlib')
            import_proxy.__dict__.update(vars(original_importlib))

            def import_module(name, package=None, fallback=original_importlib.import_module):
                return proxies[name] if package is None and name in proxies else fallback(name, package)

            import_proxy.import_module = import_module

            def historical_import(name, globals=None, locals=None, fromlist=(), level=0,
                                  fallback=original_import, library=import_proxy):
                if level == 0 and name in proxies:
                    return proxies[name]
                if level == 0 and name == 'importlib':
                    return library
                return fallback(name, globals, locals, fromlist, level)

            bindings['__import__'] = historical_import

    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, ClosedGRResponseLoader(self, key),
                                              origin=str(self.store.script_actual(key)))


class ExecutionLayout(geometry.ExecutionLayout):
    def actual(self, path, view='legacy'):
        return self.store.response_closure_actual(super().actual(path, view))


@contextmanager
def frozen_runner_routes():
    old_run = subprocess.run
    # The flag identifies the original check, not merely the newest manifest.
    modes = {
        'scripts/run_research_gr_response_active.py': {
            '--verify-navigation': '--verify-navigation',
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
                    flag = modes[target].get(flag, flag)
                    command = [*command[:4], str(THIS), flag, *command[6:]]
                    if kwargs.get('cwd') is not None:
                        # The manifest-backed store is supplied by the inner
                        # execution_routes adapter for actual script children.
                        name = relative(kwargs['cwd'])
                        if inside(name, OLD):
                            kwargs['cwd'] = str(ROOT / (NEW + name[len(OLD):]))
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
    group.add_argument('--verify-closure', action='store_true')
    group.add_argument('--verify-navigation', action='store_true')
    group.add_argument('--verify-layout', action='store_true')
    for flag in ('--verify-geometry-closure', '--verify-previous-navigation',
                 '--verify-prior-closure', '--verify-gr-navigation', '--verify-local-closure'):
        group.add_argument(flag, action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = ClosedGRResponseLayout()
    runtime = ClosedGRResponseRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    routes.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), routes.execution_routes(ExecutionLayout(store), runtime), frozen_runner_routes():
        if args.script:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)
        else:
            result = (geometry.ClosedThermalGeometryLayout.verify_closure(store) if args.verify_geometry_closure else
                      store.verify_active_navigation() if args.verify_previous_navigation else
                      thermal.ThermalLayout.verify_closure(store) if args.verify_prior_closure else
                      store.verify_navigation() if args.verify_gr_navigation else
                      store.verify_local_closure() if args.verify_local_closure else
                      store.verify_layout() if args.verify_layout else
                      store.verify_closure() if args.verify_closure else store.verify_response_navigation())
            print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
