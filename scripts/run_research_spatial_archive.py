"""Read-only replay after physically closing archive_1063_ as 1063--1085.

Current files are verified before reversible Markdown path edits are undone.
The spatial-stage scripts retain their own pre-closure view, while earlier
scripts retain their existing historical views. No archived file is rewritten.
"""
from __future__ import annotations

import argparse
import builtins
from contextlib import contextmanager
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys
import types

import run_research_live as live_runtime
from research_layout import (Layout, ResearchRuntime, ResearchLoader, ROOT,
                             relocation_original_bytes, no_archive_or_research_writes)

THIS = Path(__file__).resolve()
OLD = 'research_cognition_physics/archive_1063_'
NEW = 'research_cognition_physics/archive_1063_1085'
MANIFEST = ROOT / 'research_cognition_physics/_migration/closure_1063_1085_20261010/manifest.json'


class SpatialLayout(live_runtime.LiveLayout):
    def __init__(self):
        self.spatial, self.spatial_reverse = {}, {}
        self.spatial_manifest = json.loads(MANIFEST.read_text(encoding='utf-8'))
        if self.spatial_manifest['schema'] != 'spatial_stage_closure_v1':
            raise ValueError('Unsupported spatial closure manifest')
        for entry in self.spatial_manifest['entries']:
            logical, destination = entry['original'], entry['destination']
            actual = ROOT / destination
            if not actual.resolve().is_relative_to(ROOT) or '..' in Path(logical).parts:
                raise ValueError('Closure path escapes project: ' + destination)
            moved = logical.startswith(OLD + '/')
            if moved and destination != NEW + logical[len(OLD):]:
                raise ValueError('Unexpected closure destination: ' + destination)
            if not moved and (logical != destination or not logical.endswith('.md')):
                raise ValueError('External closure entries must be in-place Markdown edits')
            physical = os.path.normcase(os.path.abspath(actual))
            if logical in self.spatial or physical in self.spatial_reverse:
                raise ValueError('Duplicate closure entry: ' + logical)
            self.spatial[logical] = (actual, entry)
            self.spatial_reverse[physical] = logical
        # Only this process redirects the old, byte-frozen runtime's layer
        # location. Its JSON and Python source remain byte-identical.
        previous = live_runtime.LAYER
        live_runtime.LAYER = ROOT / NEW / '_shared/navigation_layer.json'
        try:
            super().__init__()
        finally:
            live_runtime.LAYER = previous
        for logical in self.spatial:
            self.current_files.add(logical)
            self.current_directories.update(p.as_posix() for p in Path(logical).parents)

    @staticmethod
    def spatial_key(path):
        try:
            name = Path(os.path.abspath(str(path))).relative_to(ROOT).as_posix()
        except ValueError:
            return None
        for prefix in (OLD, NEW):
            if name == prefix or name.startswith(prefix + '/'):
                return OLD + name[len(prefix):]
        return None

    def key(self, path, view='legacy'):
        spatial = self.spatial_key(path)
        if spatial:
            return spatial
        if view == 'spatial':
            try:
                return Path(os.path.abspath(str(path))).relative_to(ROOT).as_posix()
            except ValueError:
                return None
        return super().key(path, view)

    def actual(self, path, view='legacy'):
        spatial = self.spatial_key(path)
        if spatial:
            return ROOT / (NEW + spatial[len(OLD):])
        return Path(str(path)) if view == 'spatial' else super().actual(path, view)

    def physical_bytes(self, path):
        """Outermost read boundary for a future independently verified format layer.

        An outer layer must verify actual current bytes and its inverse hash
        before returning the exact bytes represented by this closure manifest.
        """
        return Path(str(path)).read_bytes()

    def before_spatial(self, path):
        actual = self.actual(path, 'spatial')
        raw = self.physical_bytes(actual)
        logical = self.spatial_reverse.get(os.path.normcase(os.path.abspath(actual)))
        if logical is None:
            return raw
        original = relocation_original_bytes(raw, self.spatial[logical][1])
        if hasattr(self, 'read_paths'):
            name = actual.relative_to(ROOT).as_posix()
            self.read_paths.add(name)
            if raw != original:
                self.restored_markdown.add(name)
        return original

    def before_live(self, path):
        actual = self.actual(path, 'current')
        raw = self.before_spatial(actual)
        entry = self.live.get(os.path.normcase(os.path.abspath(actual)))
        if entry is None:
            return raw
        original = relocation_original_bytes(raw, entry)
        self.read_paths.add(entry['destination'])
        if raw != original:
            self.restored_markdown.add(entry['destination'])
        return original

    def read(self, path, historical=True, view='legacy'):
        if view == 'spatial':
            return self.before_spatial(path) if historical else self.physical_bytes(self.actual(path, view))
        return super().read(path, historical=historical, view=view)

    def script_key(self, path):
        key = self.spatial_key(path)
        if key in self.spatial and key.endswith('.py'):
            return key
        return super().script_key(path)

    def script_actual(self, key):
        return self.spatial[key][0] if key in self.spatial else super().script_actual(key)

    def lookup_script(self, name):
        for candidate in (ROOT / name, ROOT / 'research_cognition_physics' / name):
            key = self.script_key(candidate)
            if key:
                return key
        return super().lookup_script(name)

    def glob_paths(self, directory, pattern, view):
        if view == 'spatial' or self.spatial_key(directory):
            physical = self.actual(directory, view)
            if physical.is_dir():
                for path in sorted(physical.glob(pattern)):
                    key = self.spatial_key(path)
                    yield ROOT / key if key else path
            return
        yield from super().glob_paths(directory, pattern, view)

    def path_class(self, view='legacy'):
        inherited = super().path_class(view) if view != 'spatial' else Layout.path_class(self, view)
        if view != 'spatial':
            return inherited
        native = type(Path())

        class SpatialPath(inherited):
            def __init__(self, *parts):
                raw = native(*parts)
                if raw.is_absolute():
                    try:
                        relative = raw.relative_to(ROOT).as_posix()
                    except ValueError:
                        relative = ''
                    for prefix in (OLD, NEW):
                        if relative == prefix or relative.startswith(prefix + '/'):
                            native.__init__(self, ROOT / (OLD + relative[len(prefix):]))
                            return
                # Preserve normal pathlib spelling, including lexical '..',
                # for this stage's frozen path-derived dictionary keys.
                native.__init__(self, raw)

        return SpatialPath

    def verify_spatial_closure(self):
        identical = edited = scientific = moved = 0
        expected = set()
        for actual, entry in self.spatial.values():
            raw = self.physical_bytes(actual)
            original = relocation_original_bytes(raw, entry)
            identical += raw == original
            edited += raw != original
            if actual.suffix in ('.py', '.json'):
                if raw != original:
                    raise ValueError('Scientific bytes changed: ' + str(actual))
                scientific += 1
            if entry['original'].startswith(OLD + '/'):
                moved += 1
                expected.add(entry['destination'])
        if (ROOT / OLD).exists():
            raise ValueError('An old directory shell must not remain')
        actual_files = {p.relative_to(ROOT).as_posix()
                        for p in (ROOT / NEW).rglob('*') if p.is_file()}
        if actual_files != expected:
            raise ValueError('Closed-stage inventory differs from manifest')
        for name, expected_hash in self.spatial_manifest['unchanged_runtime_sha256'].items():
            if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected_hash:
                raise ValueError('Previous runtime changed: ' + name)
        return dict(files=len(self.spatial), moved_files=moved,
                    unchanged_files=identical, reversible_markdown_files=edited,
                    scientific_python_json_byte_identical=scientific,
                    old_directory_absent=True, zip_used=False, directory_copy_used=False,
                    current_and_original_hashes_verified=True)

    def verify_layout(self):
        spatial = self.verify_spatial_closure()
        previous = live_runtime.LAYER
        live_runtime.LAYER = ROOT / NEW / '_shared/navigation_layer.json'
        try:
            old = super().verify_layout()
        finally:
            live_runtime.LAYER = previous
        return dict(spatial_closure=spatial, previous_layout=old, all_checks_passed=True)


class SpatialLoader(ResearchLoader):
    def exec_module(self, module):
        if self.key not in self.runtime.store.spatial:
            return super().exec_module(module)
        store = self.runtime.store
        actual = store.script_actual(self.key)
        logical = ROOT / self.key
        raw = store.read(logical, view='spatial')
        module.__file__ = str(logical)
        module.__builtins__ = self.runtime.scoped_builtins['spatial']
        previous = list(sys.path)
        sys.path.insert(0, str(logical.parent))
        try:
            exec(compile(raw, str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = previous


class SpatialRuntime(ResearchRuntime):
    def __init__(self, store):
        super().__init__(store)
        import pathlib
        proxy = types.ModuleType('pathlib')
        proxy.__dict__.update(vars(pathlib))
        proxy.Path = store.path_class('spatial')
        default_import, default_open = builtins.__import__, builtins.open

        def scoped_import(name, globals=None, locals=None, fromlist=(), level=0):
            return proxy if name == 'pathlib' and level == 0 else default_import(name, globals, locals, fromlist, level)

        def scoped_open(file, mode='r', buffering=-1, encoding=None, errors=None,
                        newline=None, closefd=True, opener=None):
            if isinstance(file, (str, bytes, os.PathLike)):
                name = os.fsdecode(file)
                key = store.key(name, 'spatial')
                if key and (key.startswith('research_cognition_physics/') or key in store.spatial):
                    if opener is not None or not closefd:
                        raise NotImplementedError('Custom historical openers are unsupported')
                    return proxy.Path(name).open(mode, buffering, encoding, errors, newline)
            return default_open(file, mode, buffering, encoding, errors, newline, closefd, opener)

        self.scoped_builtins['spatial'] = dict(vars(builtins), __import__=scoped_import, open=scoped_open)

    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, SpatialLoader(self, key),
                                              origin=str(self.store.script_actual(key)))

    def run_script(self, name, arguments=()):
        key = self.store.lookup_script(name)
        previous_argv, previous_main = sys.argv, sys.modules.get('__main__')
        module = types.ModuleType('__main__')
        module.__spec__ = self.spec('__main__', key)
        sys.modules['__main__'] = module
        sys.argv = [str(ROOT / key), *arguments]
        try:
            SpatialLoader(self, key).exec_module(module)
        finally:
            sys.argv = previous_argv
            sys.modules['__main__'] = previous_main


@contextmanager
def execution_routes(store, runtime):
    old_path, old_run = runpy.run_path, subprocess.run

    def mapped_path(path_name, init_globals=None, run_name=None):
        key = store.script_key(path_name)
        if not key:
            return old_path(path_name, init_globals=init_globals, run_name=run_name)
        name = run_name or '<run_path>'
        module = types.ModuleType(name)
        if init_globals:
            module.__dict__.update(init_globals)
        module.__name__, module.__spec__ = name, runtime.spec(name, key)
        previous = sys.modules.get(name)
        sys.modules[name] = module
        try:
            SpatialLoader(runtime, key).exec_module(module)
            return module.__dict__.copy()
        finally:
            if previous is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = previous

    def mapped_run(command, *args, **kwargs):
        if isinstance(command, (list, tuple)) and len(command) >= 5:
            same_python = os.path.normcase(os.path.abspath(str(command[0]))) == os.path.normcase(os.path.abspath(sys.executable))
            if same_python and list(command[1:4]) == ['-B', '-X', 'utf8']:
                key = store.script_key(command[4])
                if key:
                    command = [command[0], '-B', '-X', 'utf8', str(THIS), '--script', key, '--', *command[5:]]
                    if kwargs.get('cwd') is not None:
                        kwargs['cwd'] = str(store.actual(kwargs['cwd'], 'spatial'))
        return old_run(command, *args, **kwargs)

    runpy.run_path, subprocess.run = mapped_path, mapped_run
    try:
        yield
    finally:
        runpy.run_path, subprocess.run = old_path, old_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--script')
    choice.add_argument('--verify-layout', action='store_true')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store, runtime = SpatialLayout(), None
    runtime = SpatialRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), execution_routes(store, runtime):
        if args.verify_layout:
            print(json.dumps(store.verify_layout(), ensure_ascii=False))
        else:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)


if __name__ == '__main__':
    main()
