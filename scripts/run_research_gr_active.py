"""Read-only replay after closing archive_1159_ as archive_1159_1160.

Old/new closure paths both select each module's original file/Markdown view.
Optional 1161 navigation is a separate reversible outer layer. Current 1161
scripts retain actual paths. No global pathlib patch or archive copy is used.
"""
from __future__ import annotations

import argparse
import builtins
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
import re
import sys
import types

import run_research_local_dynamics_active as previous
from research_layout import ROOT, Layout, ResearchLoader, _glob_match
from research_layout import no_archive_or_research_writes, relocation_original_bytes

combined = previous.previous
digest, relative, safe_relative = combined.digest, combined.relative, combined.safe_relative
THIS = Path(__file__).resolve()
BASE = 'research_cognition_physics'
OLD = BASE + '/archive_1159_'
NEW = BASE + '/archive_1159_1160'
CURRENT = BASE + '/archive_1161_'
MANIFEST = ROOT / BASE / '_migration/closure_1159_1160/manifest.json'
GR_LAYER = ROOT / BASE / '_migration/gr_active_1161/layer.json'
NAVIGATION = previous.NAVIGATION


def in_stage(name, prefix=OLD):
    return name is not None and (name == prefix or name.startswith(prefix + '/'))


def validate_entry(entry):
    old, new = safe_relative(entry['original']), safe_relative(entry['destination'])
    for key in ('original_sha256', 'current_sha256'):
        if not isinstance(entry.get(key), str) or not re.fullmatch(r'[0-9a-f]{64}', entry[key]):
            raise ValueError('Invalid closure hash: ' + old)
    for key in ('original_bytes', 'current_bytes'):
        if type(entry.get(key)) is not int or entry[key] < 0:
            raise ValueError('Invalid closure size: ' + old)
    if not isinstance(entry.get('edits'), list):
        raise ValueError('Missing closure edits: ' + old)
    for edit in entry['edits']:
        if (type(edit.get('start')) is not int or type(edit.get('end')) is not int
                or edit['start'] < 0 or edit['end'] < edit['start']
                or not isinstance(edit.get('before'), str) or not isinstance(edit.get('after'), str)):
            raise ValueError('Invalid closure edit: ' + old)
    if not new.endswith('.md') and (entry['edits']
            or entry['original_sha256'] != entry['current_sha256']
            or entry['original_bytes'] != entry['current_bytes']):
        raise ValueError('Non-Markdown evidence changed: ' + old)
    if Path(new).suffix.lower() in ('.zip', '.pyc', '.pyo'):
        raise ValueError('Unsupported closure asset: ' + new)
    return old, new


class GRLayout(previous.LocalDynamicsLayout):
    def __init__(self):
        raw = MANIFEST.read_bytes()
        self.closure_sha256 = digest(raw)
        self.local_closure = json.loads(raw)
        if self.local_closure.get('schema') != 'research_closure_1159_1160_v1':
            raise ValueError('Unsupported 1159--1160 closure manifest')
        self.closure_entries, self.closure_reverse = {}, {}
        old_names, new_names = set(), set()
        for entry in self.local_closure['entries']:
            old, new = validate_entry(entry)
            if in_stage(old):
                if new != NEW + old[len(OLD):]:
                    raise ValueError('Unexpected closure destination: ' + old)
            elif old != new or not old.endswith('.md'):
                raise ValueError('External closure changes must be in-place Markdown')
            if os.path.normcase(old) in old_names or os.path.normcase(new) in new_names:
                raise ValueError('Duplicate closure mapping')
            old_names.add(os.path.normcase(old)); new_names.add(os.path.normcase(new))
            self.closure_entries[old], self.closure_reverse[new] = entry, old
        self.closure_directories = {p.as_posix() for name in self.closure_entries
                                    for p in PurePosixPath(name).parents}
        self.gr_entries, self.gr_layer, self.gr_sha256 = {}, None, None
        if GR_LAYER.exists():
            raw = GR_LAYER.read_bytes()
            self.gr_layer, self.gr_sha256 = json.loads(raw), digest(raw)
            if self.gr_layer.get('schema') != 'gr_active_navigation_v1':
                raise ValueError('Unsupported GR navigation layer')
            for entry in self.gr_layer['entries']:
                name, new = validate_entry(entry)
                if name not in NAVIGATION or new != name or name in self.gr_entries:
                    raise ValueError('GR navigation may change exactly four root files')
                boundary = self.closure_entries[name]
                if (entry['original_sha256'] != boundary['current_sha256']
                        or entry['original_bytes'] != boundary['current_bytes']):
                    raise ValueError('GR navigation does not join closure boundary')
                self.gr_entries[name] = entry
            if set(self.gr_entries) != NAVIGATION:
                raise ValueError('Incomplete GR navigation')
        super().__init__()
        self.current_files.update(self.closure_entries)
        self.current_directories.update(self.closure_directories)
        for name in NAVIGATION:
            entry, boundary = self.closure_entries[name], self.dynamics_entries[name]
            if (entry['original_sha256'] != boundary['current_sha256']
                    or entry['original_bytes'] != boundary['current_bytes']):
                raise ValueError('Closure does not join frozen local-dynamics navigation')

    def closure_name(self, path):
        name = relative(path)
        if name in self.closure_reverse:
            return self.closure_reverse[name]
        return OLD + name[len(NEW):] if in_stage(name, NEW) else name

    def closure_logical(self, path):
        name = self.closure_name(path)
        return ROOT / name if name is not None else Path(str(path))

    def closure_actual(self, path):
        name = self.closure_name(path)
        if name in self.closure_entries:
            return ROOT / self.closure_entries[name]['destination']
        return ROOT / (NEW + name[len(OLD):]) if in_stage(name) else Path(str(path))

    def before_gr(self, path):
        actual = Path(str(path))
        raw = actual.read_bytes()
        entry = self.gr_entries.get(relative(actual))
        return relocation_original_bytes(raw, entry) if entry else raw

    def before_local_closure(self, path):
        name = self.closure_name(path)
        if in_stage(name) and name not in self.closure_entries:
            raise FileNotFoundError('Historical file absent from closure manifest: ' + name)
        actual = self.closure_actual(path)
        raw = self.before_gr(actual)
        entry = self.closure_entries.get(name)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def before_dynamics(self, path):
        raw = self.before_local_closure(path)
        entry = self.dynamics_entries.get(self.closure_name(path))
        return relocation_original_bytes(raw, entry) if entry else raw

    def key(self, path, view='legacy'):
        logical = self.closure_logical(path)
        return relative(logical) if view == 'dynamics_closed' else super().key(logical, view)

    def actual(self, path, view='legacy'):
        if view == 'dynamics_closed':
            return self.closure_actual(path)
        return self.closure_actual(super().actual(self.closure_logical(path), view))

    def read(self, path, historical=True, view='legacy'):
        if view == 'dynamics_closed':
            return self.before_local_closure(path) if historical else self.closure_actual(path).read_bytes()
        return super().read(self.closure_logical(path), historical=historical, view=view)

    def glob_paths(self, directory, pattern, view):
        name = self.closure_name(directory)
        if view != 'dynamics_closed' and not in_stage(name):
            yield from super().glob_paths(self.closure_logical(directory), pattern, view)
            return
        prefix = name.rstrip('/') + '/' if name not in (None, '.') else ''
        found = {candidate: ROOT / candidate
                 for candidate in set(self.closure_entries) | self.closure_directories
                 if candidate.startswith(prefix) and _glob_match(candidate[len(prefix):], pattern)}
        # Within the closed archive the manifest is the complete historical
        # inventory, not an unrestricted scan of any later administrative files.
        if not in_stage(name):
            native = Path(str(self.closure_logical(directory)))
            if native.is_dir():
                for actual in native.glob(pattern):
                    physical = relative(actual)
                    if in_stage(physical, NEW):
                        continue
                    logical = self.closure_name(actual)
                    if logical and logical.startswith(prefix) and _glob_match(logical[len(prefix):], pattern):
                        found[logical] = ROOT / logical
                    elif name is None:
                        found[str(actual)] = actual
        yield from (found[key] for key in sorted(found))

    def path_class(self, view='legacy'):
        return Layout.path_class(self, view) if view == 'dynamics_closed' else super().path_class(view)

    def script_key(self, path):
        name = self.closure_name(path)
        if in_stage(name) and name in self.closure_entries and name.endswith('.py'):
            return name
        if in_stage(name, CURRENT) and name.endswith('.py') and Path(str(path)).is_file():
            return name
        return super().script_key(self.closure_logical(path))

    def script_actual(self, key):
        if in_stage(key):
            return self.closure_actual(ROOT / key)
        if in_stage(key, CURRENT):
            return ROOT / key
        return self.closure_actual(super().script_actual(key))

    def verify_navigation(self):
        if self.gr_layer is None:
            return dict(status='NOT_INSTALLED', files=0, closure_navigation_active=True)
        if digest(GR_LAYER.read_bytes()) != self.gr_sha256:
            raise ValueError('GR navigation changed during verification')
        for name in self.gr_entries:
            self.before_gr(ROOT / name)
        for name, expected in self.gr_layer['unchanged_support_sha256'].items():
            if digest((ROOT / safe_relative(name)).read_bytes()) != expected:
                raise ValueError('GR frozen support changed: ' + name)
        return dict(status='PASS', files=4, manifest_sha256=self.gr_sha256,
                    original_navigation_recoverable=True, closure_boundary_verified=True)

    def verify_local_closure(self):
        if digest(MANIFEST.read_bytes()) != self.closure_sha256:
            raise ValueError('Closure manifest changed during replay')
        expected, scientific, markdown = set(), 0, 0
        for name, entry in self.closure_entries.items():
            raw = self.before_local_closure(ROOT / name)
            if in_stage(name):
                expected.add(entry['destination'])
            scientific += Path(name).suffix.lower() in ('.py', '.json')
            markdown += bool(entry['edits'])
            if Path(name).suffix.lower() in ('.py', '.json') and digest(raw) != entry['current_sha256']:
                raise ValueError('Scientific evidence changed: ' + name)
        inventory = {p.relative_to(ROOT).as_posix() for p in (ROOT / NEW).rglob('*') if p.is_file()}
        if inventory != expected or (ROOT / OLD).exists():
            raise ValueError('Closed archive inventory differs from manifest')
        for name, expected_hash in self.local_closure['unchanged_support_sha256'].items():
            if digest((ROOT / safe_relative(name)).read_bytes()) != expected_hash:
                raise ValueError('Frozen closure support changed: ' + name)
        return dict(status='PASS', files=len(self.closure_entries), moved_files=len(expected),
                    scientific_python_json_byte_identical=scientific, reversible_markdown_files=markdown,
                    frozen_support_files=len(self.local_closure['unchanged_support_sha256']),
                    manifest_sha256=self.closure_sha256, old_directory_absent=True,
                    current_and_original_hashes_verified=True, directory_copy_used=False, zip_used=False)

    def verify_layout(self):
        return dict(gr_navigation=self.verify_navigation(), closed_local_dynamics=self.verify_local_closure(),
                    frozen_layers=super().verify_layout(), all_checks_passed=True)


class GRLoader(ResearchLoader):
    def exec_module(self, module):
        if not in_stage(self.key) and not in_stage(self.key, CURRENT):
            return previous.LocalDynamicsLoader(self.runtime, self.key).exec_module(module)
        store = self.runtime.store
        logical, actual = ROOT / self.key, store.script_actual(self.key)
        if in_stage(self.key):
            module.__builtins__ = self.runtime.scoped_builtins['dynamics_closed']
            raw = store.read(logical, view='dynamics_closed')
        else:
            raw = actual.read_bytes()
            logical = actual
        module.__file__ = str(logical)
        old_path = list(sys.path)
        sys.path.insert(0, str(logical.parent))
        try:
            exec(compile(raw, str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = old_path


class GRRuntime(previous.LocalDynamicsRuntime):
    def __init__(self, store):
        super().__init__(store)
        import pathlib
        proxy = types.ModuleType('pathlib')
        proxy.__dict__.update(vars(pathlib))
        proxy.Path = store.path_class('dynamics_closed')
        default_import, default_open = builtins.__import__, builtins.open

        def scoped_import(name, globals=None, locals=None, fromlist=(), level=0):
            return proxy if name == 'pathlib' and level == 0 else default_import(
                name, globals, locals, fromlist, level)

        def scoped_open(file, mode='r', buffering=-1, encoding=None, errors=None,
                        newline=None, closefd=True, opener=None):
            if isinstance(file, (str, bytes, os.PathLike)):
                name = os.fsdecode(file)
                key = store.key(name, 'dynamics_closed')
                if key and (key.startswith(BASE + '/') or key in store.closure_entries):
                    if opener is not None or not closefd:
                        raise NotImplementedError('Custom historical openers are unsupported')
                    return proxy.Path(name).open(mode, buffering, encoding, errors, newline)
            return default_open(file, mode, buffering, encoding, errors, newline, closefd, opener)
        self.scoped_builtins['dynamics_closed'] = dict(vars(builtins), __import__=scoped_import, open=scoped_open)

    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, GRLoader(self, key),
                                              origin=str(self.store.script_actual(key)))


class ExecutionLayout(combined.ExecutionLayout):
    def actual(self, path, view='legacy'):
        return self.store.closure_actual(super().actual(path, view))


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
    store = GRLayout()
    runtime = GRRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    routes = combined.merged.previous.previous
    routes.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), routes.execution_routes(ExecutionLayout(store), runtime):
        if args.script:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)
        else:
            result = (store.verify_layout() if args.verify_layout else
                      store.verify_local_closure() if args.verify_closure else store.verify_navigation())
            print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
