"""Read-only replay after closing the 1161--1200 and 1201--1218 stages.

The outer manifest reverses only recorded Markdown changes. Python/JSON stay
byte-identical, old and new paths select the original logical module view,
and earlier frozen layout layers keep their own readers and manifests.
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

import run_research_gr_active as previous
from research_layout import ROOT, Layout, ResearchLoader, _glob_match
from research_layout import no_archive_or_research_writes, relocation_original_bytes

THIS = Path(__file__).resolve()
BASE = 'research_cognition_physics'
MANIFEST = ROOT / BASE / '_migration/closure_1161_1218/manifest.json'
STAGES = {
    BASE + '/archive_1161_': (BASE + '/archive_1161_1200', 1161, 1200),
    BASE + '/archive_1201_': (BASE + '/archive_1201_1218', 1201, 1218),
}
VIEW = 'thermodynamics_closed'
digest, relative, safe_relative = previous.digest, previous.relative, previous.safe_relative


def inside(name, prefix):
    return name is not None and (name == prefix or name.startswith(prefix + '/'))


def moved(name):
    return any(inside(name, prefix) for prefix in STAGES)


class ThermalLayout(previous.GRLayout):
    def __init__(self):
        raw = MANIFEST.read_bytes()
        self.thermal_sha256, self.thermal_manifest = digest(raw), json.loads(raw)
        if self.thermal_manifest.get('schema') != 'research_closure_1161_1218_v1':
            raise ValueError('Unsupported thermodynamics closure manifest')
        declared = {x['source']: (x['destination'], x['first_round'], x['last_round'])
                    for x in self.thermal_manifest['stages']}
        if declared != STAGES or len(self.thermal_manifest['stages']) != len(STAGES):
            raise ValueError('Unexpected closed-stage boundaries')
        self.thermal_entries, self.thermal_reverse = {}, {}
        old_names, new_names = set(), set()
        for entry in self.thermal_manifest['entries']:
            old, new = previous.validate_entry(entry)
            expected = next((dest + old[len(source):] for source, (dest, _, _) in STAGES.items()
                             if inside(old, source)), old)
            if new != expected or (not moved(old) and (old != new or not old.endswith('.md'))):
                raise ValueError('Unexpected thermal closure destination: ' + old)
            if os.path.normcase(old) in old_names or os.path.normcase(new) in new_names:
                raise ValueError('Duplicate thermal closure mapping')
            old_names.add(os.path.normcase(old)); new_names.add(os.path.normcase(new))
            self.thermal_entries[old], self.thermal_reverse[new] = entry, old
        self.thermal_directories = {p.as_posix() for name in self.thermal_entries
                                    for p in PurePosixPath(name).parents}
        super().__init__()
        self.current_files.update(self.thermal_entries)
        self.current_directories.update(self.thermal_directories)

    def thermal_name(self, path):
        name = relative(path)
        if name in self.thermal_reverse:
            return self.thermal_reverse[name]
        for source, (destination, _, _) in STAGES.items():
            if inside(name, destination):
                return source + name[len(destination):]
        return name

    def thermal_logical(self, path):
        name = self.thermal_name(path)
        return ROOT / name if name is not None else Path(str(path))

    def thermal_actual(self, path):
        name = self.thermal_name(path)
        if name in self.thermal_entries:
            return ROOT / self.thermal_entries[name]['destination']
        for source, (destination, _, _) in STAGES.items():
            if inside(name, source):
                return ROOT / (destination + name[len(source):])
        return Path(str(path))

    def before_thermal(self, path):
        name = self.thermal_name(path)
        if moved(name) and name not in self.thermal_entries:
            raise FileNotFoundError('File absent from thermal closure inventory: ' + name)
        raw = self.thermal_actual(path).read_bytes()
        entry = self.thermal_entries.get(name)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def before_gr(self, path):
        raw = self.before_thermal(path)
        entry = self.gr_entries.get(self.thermal_name(path))
        return relocation_original_bytes(raw, entry) if entry else raw

    def key(self, path, view='legacy'):
        logical = self.thermal_logical(path)
        return relative(logical) if view == VIEW else super().key(logical, view)

    def actual(self, path, view='legacy'):
        if view == VIEW:
            return self.thermal_actual(path)
        return self.thermal_actual(super().actual(self.thermal_logical(path), view))

    def read(self, path, historical=True, view='legacy'):
        if view == VIEW:
            return self.before_thermal(path) if historical else self.thermal_actual(path).read_bytes()
        return super().read(self.thermal_logical(path), historical=historical, view=view)

    def glob_paths(self, directory, pattern, view):
        name = self.thermal_name(directory)
        if view != VIEW and not moved(name):
            yield from super().glob_paths(self.thermal_logical(directory), pattern, view)
            return
        prefix = name.rstrip('/') + '/' if name not in (None, '.') else ''
        found = {candidate: ROOT / candidate
                 for candidate in set(self.thermal_entries) | self.thermal_directories
                 if candidate.startswith(prefix) and _glob_match(candidate[len(prefix):], pattern)}
        if not moved(name):
            native = self.thermal_actual(directory)
            if native.is_dir():
                for actual in native.glob(pattern):
                    logical = self.thermal_name(actual)
                    if logical and logical.startswith(prefix) and _glob_match(logical[len(prefix):], pattern):
                        found[logical] = ROOT / logical
                    elif name is None:
                        found[str(actual)] = actual
        yield from (found[key] for key in sorted(found))

    def path_class(self, view='legacy'):
        return Layout.path_class(self, view) if view == VIEW else super().path_class(view)

    def script_key(self, path):
        name = self.thermal_name(path)
        if moved(name) and name in self.thermal_entries and name.endswith('.py'):
            return name
        return super().script_key(self.thermal_logical(path))

    def script_actual(self, key):
        return self.thermal_actual(ROOT / key) if moved(key) else self.thermal_actual(super().script_actual(key))

    def verify_closure(self):
        if digest(MANIFEST.read_bytes()) != self.thermal_sha256:
            raise ValueError('Thermal manifest changed during replay')
        expected, scientific, markdown = set(), 0, 0
        for name, entry in self.thermal_entries.items():
            raw = self.before_thermal(ROOT / name)
            if moved(name):
                expected.add(entry['destination'])
            scientific += Path(name).suffix.lower() in ('.py', '.json')
            markdown += bool(entry['edits'])
            if Path(name).suffix.lower() in ('.py', '.json') and digest(raw) != entry['current_sha256']:
                raise ValueError('Scientific evidence changed: ' + name)
        actual = set()
        for source, (destination, _, _) in STAGES.items():
            if (ROOT / source).exists() or not (ROOT / destination).is_dir():
                raise ValueError('Old directory remains or new directory is absent: ' + source)
            actual.update(p.relative_to(ROOT).as_posix() for p in (ROOT / destination).rglob('*') if p.is_file())
        if actual != expected:
            raise ValueError('Thermal closure inventory differs from manifest')
        support = self.thermal_manifest['unchanged_support_sha256']
        if not isinstance(support, dict) or not support:
            raise ValueError('Frozen support bindings are required')
        for name, expected_hash in support.items():
            if digest((ROOT / safe_relative(name)).read_bytes()) != expected_hash:
                raise ValueError('Frozen support changed: ' + name)
        for name, expected_hash in self.thermal_manifest.get('new_support_sha256', {}).items():
            if digest((ROOT / safe_relative(name)).read_bytes()) != expected_hash:
                raise ValueError('New support changed: ' + name)
        return dict(status='PASS',files=len(self.thermal_entries),moved_files=len(expected),
                    scientific_python_json_byte_identical=scientific,reversible_markdown_files=markdown,
                    frozen_support_files=len(support),manifest_sha256=self.thermal_sha256,
                    old_directories_absent=True,current_and_original_hashes_verified=True,
                    directory_copy_used=False,zip_used=False)

    def verify_layout(self):
        return dict(thermal_closure=self.verify_closure(),frozen_layers=super().verify_layout(),
                    all_checks_passed=True)


class ThermalLoader(ResearchLoader):
    def exec_module(self, module):
        if not moved(self.key):
            return previous.GRLoader(self.runtime, self.key).exec_module(module)
        store = self.runtime.store
        logical, actual = ROOT / self.key, store.script_actual(self.key)
        module.__file__, module.__builtins__ = str(logical), self.runtime.scoped_builtins[VIEW]
        old_path = list(sys.path)
        sys.path.insert(0, str(logical.parent))
        try:
            exec(compile(store.read(logical, view=VIEW), str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = old_path


class ThermalRuntime(previous.GRRuntime):
    def __init__(self, store):
        super().__init__(store)
        import pathlib
        proxy = types.ModuleType('pathlib')
        proxy.__dict__.update(vars(pathlib))
        proxy.Path = store.path_class(VIEW)
        default_import, default_open = builtins.__import__, builtins.open

        def scoped_import(name, globals=None, locals=None, fromlist=(), level=0):
            return proxy if name == 'pathlib' and level == 0 else default_import(
                name, globals, locals, fromlist, level)

        def scoped_open(file, mode='r', buffering=-1, encoding=None, errors=None,
                        newline=None, closefd=True, opener=None):
            if isinstance(file, (str, bytes, os.PathLike)):
                name = os.fsdecode(file)
                key = store.key(name, VIEW)
                if key and (key.startswith(BASE + '/') or key in store.thermal_entries):
                    if opener is not None or not closefd:
                        raise NotImplementedError('Custom historical openers are unsupported')
                    return proxy.Path(name).open(mode, buffering, encoding, errors, newline)
            return default_open(file, mode, buffering, encoding, errors, newline, closefd, opener)
        self.scoped_builtins[VIEW] = dict(vars(builtins), __import__=scoped_import, open=scoped_open)

    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, ThermalLoader(self, key),
                                              origin=str(self.store.script_actual(key)))


class ExecutionLayout(previous.ExecutionLayout):
    def actual(self, path, view='legacy'):
        return self.store.thermal_actual(super().actual(path, view))


@contextmanager
def gr_verification_route():
    """Forward the frozen finalizer's known GR verification child command."""
    old_run = subprocess.run

    def mapped_run(command, *args, **kwargs):
        if isinstance(command, (list, tuple)) and len(command) == 6:
            same_python = os.path.normcase(os.path.abspath(str(command[0]))) == os.path.normcase(sys.executable)
            same_script = relative(command[4]) == 'scripts/run_research_gr_active.py'
            if same_python and same_script and list(command[1:4]) == ['-B', '-X', 'utf8']:
                if command[5] in ('--verify-navigation', '--verify-layout'):
                    command = [command[0], '-B', '-X', 'utf8', str(THIS), command[5]]
        return old_run(command, *args, **kwargs)
    subprocess.run = mapped_run
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
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store, runtime = ThermalLayout(), None
    runtime = ThermalRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    routes = previous.combined.merged.previous.previous
    routes.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), routes.execution_routes(ExecutionLayout(store), runtime), gr_verification_route():
        if args.script:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)
        else:
            result = (store.verify_layout() if args.verify_layout else
                      store.verify_closure() if args.verify_closure else store.verify_navigation())
            print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
