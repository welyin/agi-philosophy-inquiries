"""Read-only replay after closing archive_1086_ as archive_1086_1105.

The outer boundary verifies real files and reverses only the recorded closure
edits. Earlier active, formula and archive layers remain byte-frozen. The SR
scripts receive their pre-closure paths and bytes; older scripts retain their
own historical views. No directory copies, ZIPs or old directory shells are used.
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

import run_research_active as active_runtime
from research_layout import ROOT, Layout, ResearchLoader, relocation_original_bytes, no_archive_or_research_writes

spatial = active_runtime.formatted.spatial
THIS = Path(__file__).resolve()
OLD = 'research_cognition_physics/archive_1086_'
NEW = 'research_cognition_physics/archive_1086_1105'
MAINTENANCE = ROOT / 'research_cognition_physics/_migration/sr_closure_1086_1105'
MANIFEST = MAINTENANCE / 'manifest.json'
PUBLICATION_RECEIPT = NEW + '/_shared/stage_paper_closure/publication_checks.json'


def path_key(path):
    return os.path.normcase(os.path.abspath(str(path)))


class SRLayout(active_runtime.ActiveLayout):
    def __init__(self):
        raw = MANIFEST.read_bytes()
        self.sr_closure_manifest_sha256 = hashlib.sha256(raw).hexdigest()
        self.sr_closure_manifest = json.loads(raw.decode('utf-8-sig'))
        if self.sr_closure_manifest.get('schema') != 'sr_stage_closure_v1':
            raise ValueError('Unsupported SR closure manifest')
        self.sr_closure, self.sr_closure_reverse = {}, {}
        for entry in self.sr_closure_manifest['entries']:
            logical, destination = entry['original'], entry['destination']
            actual = ROOT / destination
            if (not actual.resolve().is_relative_to(ROOT) or '..' in Path(logical).parts
                    or '\\' in logical or Path(logical).is_absolute()):
                raise ValueError('Unsafe closure path: ' + logical)
            moved = logical.startswith(OLD + '/')
            if moved and destination != NEW + logical[len(OLD):]:
                raise ValueError('Unexpected SR destination: ' + destination)
            if not moved and (logical != destination or not logical.endswith('.md')):
                raise ValueError('External closure entries must be Markdown edits')
            key = path_key(actual)
            if logical in self.sr_closure or key in self.sr_closure_reverse:
                raise ValueError('Duplicate closure entry: ' + logical)
            self.sr_closure[logical] = (actual, entry)
            self.sr_closure_reverse[key] = logical
        prior_nav = active_runtime.formatted.NAV
        active_runtime.formatted.NAV = ROOT / NEW / '_shared/navigation_layer.json'
        try:
            super().__init__()
        finally:
            active_runtime.formatted.NAV = prior_nav
        for logical in self.sr_closure:
            self.current_files.add(logical)
            self.current_directories.update(p.as_posix() for p in Path(logical).parents)

    @staticmethod
    def sr_key(path):
        try:
            name = Path(os.path.abspath(str(path))).relative_to(ROOT).as_posix()
        except ValueError:
            return None
        for prefix in (OLD, NEW):
            if name == prefix or name.startswith(prefix + '/'):
                return OLD + name[len(prefix):]
        return None

    def key(self, path, view='legacy'):
        own = self.sr_key(path)
        if own:
            return own
        if view == 'sr':
            try:
                return Path(os.path.abspath(str(path))).relative_to(ROOT).as_posix()
            except ValueError:
                return None
        return super().key(path, view)

    def actual(self, path, view='legacy'):
        own = self.sr_key(path)
        if own:
            return ROOT / (NEW + own[len(OLD):])
        return Path(str(path)) if view == 'sr' else super().actual(path, view)

    def before_closure(self, path):
        actual = self.actual(path, 'sr')
        raw = actual.read_bytes()
        logical = self.sr_closure_reverse.get(path_key(actual))
        if logical is None:
            return raw
        prior = relocation_original_bytes(raw, self.sr_closure[logical][1])
        if hasattr(self, 'read_paths'):
            self.read_paths.add(actual.relative_to(ROOT).as_posix())
            if raw != prior:
                self.restored_markdown.add(actual.relative_to(ROOT).as_posix())
        return prior

    def before_active(self, path):
        raw = self.before_closure(path)
        own = self.sr_key(path)
        logical_path = ROOT / own if own else Path(str(path))
        item = self.active.get(path_key(logical_path))
        return relocation_original_bytes(raw, item[1]) if item else raw

    def physical_bytes(self, path):
        raw = self.before_active(path)
        own = self.sr_key(path)
        key = path_key(ROOT / own if own else path)
        if key in self.formats:
            raw = active_runtime.formatted.formatting.restore_bytes(raw, self.formats[key])
        if key in self.sr:
            raw = relocation_original_bytes(raw, self.sr[key])
        return raw

    def read(self, path, historical=True, view='legacy'):
        if view == 'sr':
            return self.before_closure(path) if historical else self.actual(path, view).read_bytes()
        return super().read(path, historical=historical, view=view)

    def script_key(self, path):
        own = self.sr_key(path)
        if own in self.sr_closure and own.endswith('.py'):
            return own
        return super().script_key(path)

    def script_actual(self, key):
        return self.sr_closure[key][0] if key in self.sr_closure else super().script_actual(key)

    def lookup_script(self, name):
        for candidate in (ROOT / name, ROOT / 'research_cognition_physics' / name):
            key = self.script_key(candidate)
            if key:
                return key
        return super().lookup_script(name)

    def glob_paths(self, directory, pattern, view):
        if view == 'sr' or self.sr_key(directory):
            physical = self.actual(directory, view)
            if physical.is_dir():
                for path in sorted(physical.glob(pattern)):
                    key = self.sr_key(path)
                    yield ROOT / key if key else path
            return
        yield from super().glob_paths(directory, pattern, view)

    def path_class(self, view='legacy'):
        if view != 'sr':
            return super().path_class(view)
        inherited, native = Layout.path_class(self, view), type(Path())

        class SRPath(inherited):
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
                # Frozen source dictionaries sometimes retain lexical '..'.
                native.__init__(self, raw)
        return SRPath

    def verify_closure_layer(self):
        if hashlib.sha256(MANIFEST.read_bytes()).hexdigest() != self.sr_closure_manifest_sha256:
            raise ValueError('Closure manifest changed during verification')
        expected, scientific, changed = set(), 0, 0
        for logical, (actual, entry) in self.sr_closure.items():
            raw = actual.read_bytes()
            prior = relocation_original_bytes(raw, entry)
            changed += raw != prior
            if actual.suffix.lower() in ('.py', '.json'):
                if raw != prior:
                    raise ValueError('Frozen scientific bytes changed: ' + logical)
                scientific += 1
            if logical.startswith(OLD + '/'):
                expected.add(entry['destination'])
        inventory = {p.relative_to(ROOT).as_posix() for p in (ROOT / NEW).rglob('*') if p.is_file()}
        allowed = self.sr_closure_manifest.get('post_publication_files', [])
        if allowed != [PUBLICATION_RECEIPT] or PUBLICATION_RECEIPT in expected:
            raise ValueError('Unexpected post-publication allowance')
        if (ROOT / OLD).exists() or not expected <= inventory or inventory - expected - set(allowed):
            raise ValueError('SR archive inventory differs from the manifest')
        for name, digest in self.sr_closure_manifest['unchanged_support_sha256'].items():
            if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
                raise ValueError('Frozen support changed: ' + name)
        return dict(files=len(self.sr_closure), moved_files=len(expected),
                    scientific_python_json_byte_identical=scientific,
                    reversible_markdown_files=changed, old_directory_absent=True,
                    zip_used=False, directory_copy_used=False,
                    post_publication_receipt_sha256={name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                                                    for name in allowed if name in inventory},
                    manifest_sha256=self.sr_closure_manifest_sha256,
                    current_and_original_hashes_verified=True)

    def verify_layout(self):
        closure = self.verify_closure_layer()
        prior = super().verify_layout()
        return {'sr_closure': closure, 'active_and_earlier_layers': prior, 'all_checks_passed': True}


class SRLoader(ResearchLoader):
    def exec_module(self, module):
        store = self.runtime.store
        if not self.key.startswith(OLD + '/'):
            return spatial.SpatialLoader(self.runtime, self.key).exec_module(module)
        actual, logical = store.script_actual(self.key), ROOT / self.key
        module.__file__ = str(logical)
        module.__builtins__ = self.runtime.scoped_builtins['sr']
        previous = list(sys.path)
        sys.path.insert(0, str(logical.parent))
        try:
            exec(compile(store.read(logical, view='sr'), str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = previous


class SRRuntime(spatial.SpatialRuntime):
    def __init__(self, store):
        super().__init__(store)
        import pathlib
        proxy = types.ModuleType('pathlib')
        proxy.__dict__.update(vars(pathlib))
        proxy.Path = store.path_class('sr')
        original_import, original_open = builtins.__import__, builtins.open

        def scoped_import(name, globals=None, locals=None, fromlist=(), level=0):
            return proxy if name == 'pathlib' and level == 0 else original_import(name, globals, locals, fromlist, level)

        def scoped_open(file, mode='r', buffering=-1, encoding=None, errors=None,
                        newline=None, closefd=True, opener=None):
            if isinstance(file, (str, bytes, os.PathLike)):
                name = os.fsdecode(file)
                key = store.key(name, 'sr')
                if key and (key.startswith('research_cognition_physics/') or key in store.sr_closure):
                    if opener is not None or not closefd:
                        raise NotImplementedError('Custom historical openers are unsupported')
                    return proxy.Path(name).open(mode, buffering, encoding, errors, newline)
            return original_open(file, mode, buffering, encoding, errors, newline, closefd, opener)
        self.scoped_builtins['sr'] = dict(vars(builtins), __import__=scoped_import, open=scoped_open)

    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, SRLoader(self, key), origin=str(self.store.script_actual(key)))

    def run_script(self, name, arguments=()):
        key = self.store.lookup_script(name)
        old_argv, old_main = sys.argv, sys.modules.get('__main__')
        module = types.ModuleType('__main__')
        module.__spec__ = self.spec('__main__', key)
        sys.modules['__main__'] = module
        sys.argv = [str(ROOT / key), *arguments]
        try:
            module.__spec__.loader.exec_module(module)
        finally:
            sys.argv, sys.modules['__main__'] = old_argv, old_main


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
            module.__spec__.loader.exec_module(module)
            return module.__dict__.copy()
        finally:
            if previous is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = previous

    def mapped_run(command, *args, **kwargs):
        if isinstance(command, (list, tuple)) and len(command) >= 5:
            same_python = path_key(command[0]) == path_key(sys.executable)
            if same_python and list(command[1:4]) == ['-B', '-X', 'utf8']:
                key = store.script_key(command[4])
                if key:
                    command = [command[0], '-B', '-X', 'utf8', str(THIS), '--script', key, '--', *command[5:]]
                    if kwargs.get('cwd') is not None:
                        kwargs['cwd'] = str(store.actual(kwargs['cwd'], 'sr'))
        return old_run(command, *args, **kwargs)
    runpy.run_path, subprocess.run = mapped_path, mapped_run
    try:
        yield
    finally:
        runpy.run_path, subprocess.run = old_path, old_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--script')
    group.add_argument('--verify-layout', action='store_true')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = SRLayout()
    runtime = SRRuntime(store)
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
