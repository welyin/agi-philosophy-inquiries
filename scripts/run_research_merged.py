"""Read-only replay through the 1106--1119 merge and paper-numbering layer.

The manifest maps real files to their exact pre-merge names and bytes. Earlier
runtime code and manifests stay frozen. Each module retains its own historical
view and __file__; no copied tree, legacy directory shell or ZIP is needed.
"""
from __future__ import annotations

import argparse
import builtins
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import sys
import types

import run_research_motion_active as previous
from research_layout import (
    ROOT, Layout, ResearchLoader, _glob_match,
    no_archive_or_research_writes, relocation_original_bytes,
)

THIS = Path(__file__).resolve()
BASE = 'research_cognition_physics'
NEW = BASE + '/archive_1106_1119'
STAGES = {
    BASE + '/archive_1106_': '1106_1115',
    BASE + '/archive_1116_': '1116_1119',
}
MANIFEST = ROOT / BASE / '_migration/merge_1106_1119/manifest.json'


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def relative(path):
    try:
        return Path(os.path.abspath(str(path))).relative_to(ROOT).as_posix()
    except ValueError:
        return None


def safe_relative(name):
    if not isinstance(name, str):
        raise ValueError('Manifest paths must be strings')
    path = PurePosixPath(name)
    if (not name or '\\' in name or ':' in name or path.is_absolute()
            or '..' in path.parts or path.as_posix() != name
            or set(path.parts) & {'.git', '.codex', '.agents', '.research_runtime'}
            or not (ROOT / name).resolve().is_relative_to(ROOT)):
        raise ValueError('Unsafe merge path: ' + name)
    return name


def moved_stage(name):
    return next((old for old in STAGES if name == old or name.startswith(old + '/')), None)


def stage_destination(name):
    old = moved_stage(name)
    if old is None:
        return name
    tail = name[len(old):].lstrip('/')
    if tail in ('README.md', '文件索引.md'):
        return NEW + '/_shared/' + STAGES[old] + '/' + tail
    if tail == '_shared' or tail.startswith('_shared/'):
        return NEW + '/_shared/' + STAGES[old] + tail[len('_shared'):]
    return NEW + ('/' + tail if tail else '')


class MergedLayout(previous.MotionLayout):
    def __init__(self):
        raw = MANIFEST.read_bytes()
        self.merge_manifest_sha256 = digest(raw)
        self.merge_manifest = json.loads(raw.decode('utf-8-sig'))
        if self.merge_manifest.get('schema') != 'research_merge_1106_1119_v1':
            raise ValueError('Unsupported research merge manifest')
        self.merge_entries, self.merge_reverse = {}, {}
        normalized_originals, normalized_destinations = set(), set()
        for entry in self.merge_manifest['entries']:
            old, new = safe_relative(entry['original']), safe_relative(entry['destination'])
            if moved_stage(old):
                if new != stage_destination(old):
                    raise ValueError('Unexpected merged-stage destination: ' + new)
            elif old != new:
                before, after = PurePosixPath(old), PurePosixPath(new)
                prefix = after.name[:-len(before.name)] if after.name.endswith(before.name) else ''
                if (before.parent.as_posix() != BASE or after.parent != before.parent
                        or not before.name.endswith('_阶段论文.md')
                        or not re.fullmatch(r'0[1-5][. _、—-]+', prefix)):
                    raise ValueError('Only five numbered root stage papers may be renamed: ' + old)
            elif not old.endswith('.md'):
                raise ValueError('External merge changes must be Markdown: ' + old)
            old_key, new_key = os.path.normcase(old), os.path.normcase(new)
            if old_key in normalized_originals or new_key in normalized_destinations:
                raise ValueError('Duplicate merge source or destination: ' + old)
            normalized_originals.add(old_key)
            normalized_destinations.add(new_key)
            for field in ('original_sha256', 'current_sha256'):
                if not isinstance(entry.get(field), str) or not re.fullmatch(r'[0-9a-f]{64}', entry[field]):
                    raise ValueError('Invalid merge hash: ' + old)
            for field in ('original_bytes', 'current_bytes'):
                if type(entry.get(field)) is not int or entry[field] < 0:
                    raise ValueError('Invalid merge byte count: ' + old)
            if not isinstance(entry.get('edits'), list):
                raise ValueError('Merge edits must be a list: ' + old)
            for edit in entry['edits']:
                if (type(edit.get('start')) is not int or type(edit.get('end')) is not int
                        or edit['start'] < 0 or edit['end'] < edit['start']
                        or not isinstance(edit.get('before'), str)
                        or not isinstance(edit.get('after'), str)):
                    raise ValueError('Invalid merge edit: ' + old)
            if not new.endswith('.md') and (entry['edits']
                    or entry['original_sha256'] != entry['current_sha256']
                    or entry['original_bytes'] != entry['current_bytes']):
                raise ValueError('Non-Markdown merge assets must be byte-identical: ' + old)
            if Path(new).suffix.lower() in ('.zip', '.pyc', '.pyo'):
                raise ValueError('Caches and ZIPs are not replay assets: ' + new)
            self.merge_entries[old] = entry
            self.merge_reverse[new] = old
        self.merge_directories = {
            parent.as_posix() for name in self.merge_entries
            for parent in PurePosixPath(name).parents
        }
        created = self.merge_manifest.get('created_files', {})
        self.created_sizes = {}
        if isinstance(created, list):
            self.created_files = {entry['path']: entry['sha256'] for entry in created}
            self.created_sizes = {entry['path']: entry['bytes'] for entry in created}
            if len(self.created_files) != len(created):
                raise ValueError('Duplicate new navigation file')
        else:
            self.created_files = created
        if not isinstance(self.created_files, dict) or set(self.created_files) != {
                NEW + '/README.md', NEW + '/文件索引.md'}:
            raise ValueError('Exactly the merged README and index must be new files')
        self.premerge_path = self._premerge_path_class()
        super().__init__()
        self.current_files.update(self.merge_entries)
        self.current_directories.update(self.merge_directories)
        # Frozen verifiers occasionally read the paths in their inventories
        # directly. Adapt those file objects at the same outer byte boundary;
        # do not alter entries, hashes, source code or global pathlib behavior.
        for attribute in ('entries', 'relocations', 'recent', 'spatial', 'sr_closure', 'active'):
            inventory = getattr(self, attribute)
            setattr(self, attribute, {
                key: (self.premerge_path(path), entry)
                for key, (path, entry) in inventory.items()
            })

    def original_name(self, path):
        name = relative(path)
        if name in self.merge_reverse:
            return self.merge_reverse[name]
        # A shared subdirectory has an unambiguous old name; the merged root
        # itself has two predecessors and is deliberately never aliased.
        if name is not None:
            for old, group in STAGES.items():
                prefix = NEW + '/_shared/' + group
                if name == prefix or name.startswith(prefix + '/'):
                    return old + '/_shared' + name[len(prefix):]
            prefix = NEW + '/'
            if name.startswith(prefix):
                first = name[len(prefix):].split('/', 1)[0]
                if first.isdigit() and 1106 <= int(first) <= 1119:
                    old = next(old for old, group in STAGES.items()
                               if int(group.split('_')[0]) <= int(first) <= int(group.split('_')[1]))
                    return old + name[len(NEW):]
        return name

    def premerge_name_path(self, path):
        name = self.original_name(path)
        return ROOT / name if name is not None else Path(str(path))

    def merge_actual(self, path):
        old = self.original_name(path)
        if old in self.merge_entries:
            return ROOT / self.merge_entries[old]['destination']
        if old and moved_stage(old):
            return ROOT / stage_destination(old)
        return Path(str(path))

    def before_merge(self, path):
        old = self.original_name(path)
        actual = self.merge_actual(path)
        raw = actual.read_bytes()
        entry = self.merge_entries.get(old)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def before_motion(self, path):
        old = self.original_name(path)
        raw = self.before_merge(path)
        entry = self.motion_entries.get(old)
        return relocation_original_bytes(raw, entry) if entry else raw

    def key(self, path, view='legacy'):
        path = self.premerge_name_path(path)
        return relative(path) if view == 'merged' else super().key(path, view)

    def actual(self, path, view='legacy'):
        if view == 'merged':
            return self.merge_actual(path)
        return self.premerge_path(super().actual(self.premerge_name_path(path), view))

    def read(self, path, historical=True, view='legacy'):
        if not historical:
            prior_location = (self.premerge_name_path(path) if view == 'merged'
                              else super().actual(self.premerge_name_path(path), view))
            return self.merge_actual(prior_location).read_bytes()
        if view == 'merged':
            return self.before_merge(path)
        return super().read(self.premerge_name_path(path), historical=True, view=view)

    def _premerge_path_class(self):
        store, native = self, type(Path())

        class BeforeMergePath(native):
            def read_bytes(self):
                return store.before_merge(self)

            def read_text(self, encoding=None, errors=None, newline=None):
                with self.open('r', encoding=encoding, errors=errors, newline=newline) as stream:
                    return stream.read()

            def open(self, mode='r', buffering=-1, encoding=None, errors=None, newline=None):
                if any(char in mode for char in 'wax+'):
                    raise PermissionError('Historical replay is read-only: ' + str(self))
                stream = io.BytesIO(self.read_bytes())
                return stream if 'b' in mode else io.TextIOWrapper(
                    stream, encoding=encoding or 'utf8', errors=errors, newline=newline)

            def exists(self, **kwargs):
                return (store.original_name(self) in store.merge_directories
                        or store.merge_actual(self).exists(**kwargs))

            def is_file(self, **kwargs):
                return store.merge_actual(self).is_file(**kwargs)

            def is_dir(self, **kwargs):
                return (store.original_name(self) in store.merge_directories
                        or store.merge_actual(self).is_dir(**kwargs))

            def glob(self, pattern, **kwargs):
                if kwargs:
                    raise NotImplementedError('Historical glob keyword options are unsupported')
                yield from (type(self)(p) for p in store.premerge_glob(self, pattern))

            def rglob(self, pattern, **kwargs):
                yield from self.glob('**/' + pattern, **kwargs)

            def iterdir(self):
                yield from self.glob('*')

        return BeforeMergePath

    def premerge_glob(self, directory, pattern):
        name = self.original_name(directory)
        prefix = name.rstrip('/') + '/' if name not in (None, '.') else ''
        found = {}
        for candidate in set(self.merge_entries) | self.merge_directories:
            if candidate.startswith(prefix) and _glob_match(candidate[len(prefix):], pattern):
                found[candidate] = ROOT / candidate
        native = Path(str(self.premerge_name_path(directory)))
        if native.is_dir():
            for actual in native.glob(pattern):
                physical = relative(actual)
                if physical == NEW or (physical and physical.startswith(NEW + '/')):
                    continue
                logical = self.original_name(actual)
                if logical and logical.startswith(prefix) and _glob_match(logical[len(prefix):], pattern):
                    found[logical] = ROOT / logical
                elif name is None:
                    found[str(actual)] = actual
        yield from (found[key] for key in sorted(found))

    def glob_paths(self, directory, pattern, view):
        if view == 'merged' or moved_stage(self.original_name(directory) or ''):
            yield from self.premerge_glob(directory, pattern)
        else:
            yield from super().glob_paths(self.premerge_name_path(directory), pattern, view)

    def path_class(self, view='legacy'):
        return Layout.path_class(self, view) if view == 'merged' else super().path_class(view)

    def script_key(self, path):
        old = self.original_name(path)
        if old in self.merge_entries and old.endswith('.py') and moved_stage(old):
            return old
        return super().script_key(self.premerge_name_path(path))

    def script_actual(self, key):
        if key in self.merge_entries and moved_stage(key):
            return ROOT / self.merge_entries[key]['destination']
        return self.merge_actual(super().script_actual(key))

    def lookup_script(self, name):
        for candidate in (ROOT / name, ROOT / BASE / name):
            key = self.script_key(candidate)
            if key:
                return key
        return super().lookup_script(name)

    def verify_closure_layer(self):
        # The previous adapter covers only the five Motion-layer documents.
        # Numbered papers and external link repairs also need the same reader.
        class ClosureReader:
            def __init__(adapter):
                adapter.sr_closure = {
                    key: (previous.BeforeMotionFile(self, path), entry)
                    for key, (path, entry) in self.sr_closure.items()
                }

            def __getattr__(adapter, name):
                return getattr(self, name)

        return previous.previous.SRLayout.verify_closure_layer(ClosureReader())

    def verify_merge(self):
        if digest(MANIFEST.read_bytes()) != self.merge_manifest_sha256:
            raise ValueError('Merge manifest changed during replay')
        moved, scientific, markdown = 0, 0, 0
        expected = set(self.created_files)
        for old, entry in self.merge_entries.items():
            prior = self.before_merge(ROOT / old)
            moved += old != entry['destination']
            markdown += entry['original_sha256'] != entry['current_sha256']
            if Path(old).suffix.lower() in ('.py', '.json'):
                if digest(prior) != entry['current_sha256']:
                    raise ValueError('Frozen scientific file changed: ' + old)
                scientific += 1
            if moved_stage(old):
                expected.add(entry['destination'])
            if old != entry['destination'] and (ROOT / old).exists():
                raise ValueError('An old file remains after relocation: ' + old)
        for name, expected_hash in self.created_files.items():
            raw = (ROOT / safe_relative(name)).read_bytes()
            if (digest(raw) != expected_hash
                    or (name in self.created_sizes and len(raw) != self.created_sizes[name])):
                raise ValueError('New navigation file changed: ' + name)
        inventory = {p.relative_to(ROOT).as_posix() for p in (ROOT / NEW).rglob('*') if p.is_file()}
        if inventory != expected or any((ROOT / old).exists() for old in STAGES):
            raise ValueError('Merged archive inventory differs from the manifest')
        support = self.merge_manifest.get('unchanged_support_sha256')
        if not isinstance(support, dict) or not support:
            raise ValueError('Frozen-support hashes are required')
        for name, expected_hash in support.items():
            if digest((ROOT / safe_relative(name)).read_bytes()) != expected_hash:
                raise ValueError('Frozen support changed: ' + name)
        return dict(files=len(self.merge_entries), moved_files=moved,
                    reversible_markdown_files=markdown,
                    scientific_python_json_byte_identical=scientific,
                    new_navigation_files=len(self.created_files),
                    frozen_support_files=len(support),
                    manifest_sha256=self.merge_manifest_sha256,
                    current_and_original_hashes_verified=True,
                    old_directories_absent=True, directory_copy_used=False, zip_used=False)

    def verify_layout(self):
        merged = self.verify_merge()
        frozen = super().verify_layout()
        return dict(merged_layout=merged, frozen_layers=frozen, all_checks_passed=True)


class MergedLoader(ResearchLoader):
    def exec_module(self, module):
        if not moved_stage(self.key):
            return previous.previous.SRLoader(self.runtime, self.key).exec_module(module)
        store = self.runtime.store
        logical, actual = ROOT / self.key, store.script_actual(self.key)
        module.__file__ = str(logical)
        module.__builtins__ = self.runtime.scoped_builtins['merged']
        old_path = list(sys.path)
        sys.path.insert(0, str(logical.parent))
        try:
            exec(compile(store.read(logical, view='merged'), str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = old_path


class MergedRuntime(previous.previous.SRRuntime):
    def __init__(self, store):
        super().__init__(store)
        import pathlib
        proxy = types.ModuleType('pathlib')
        proxy.__dict__.update(vars(pathlib))
        proxy.Path = store.path_class('merged')
        default_import, default_open = builtins.__import__, builtins.open

        def scoped_import(name, globals=None, locals=None, fromlist=(), level=0):
            return proxy if name == 'pathlib' and level == 0 else default_import(
                name, globals, locals, fromlist, level)

        def scoped_open(file, mode='r', buffering=-1, encoding=None, errors=None,
                        newline=None, closefd=True, opener=None):
            if isinstance(file, (str, bytes, os.PathLike)):
                name = os.fsdecode(file)
                key = store.key(name, 'merged')
                if key and (key.startswith(BASE + '/') or key in store.merge_entries):
                    if opener is not None or not closefd:
                        raise NotImplementedError('Custom historical openers are unsupported')
                    return proxy.Path(name).open(mode, buffering, encoding, errors, newline)
            return default_open(file, mode, buffering, encoding, errors, newline, closefd, opener)

        self.scoped_builtins['merged'] = dict(vars(builtins), __import__=scoped_import, open=scoped_open)

    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, MergedLoader(self, key),
                                              origin=str(self.store.script_actual(key)))


class ExecutionLayout:
    """Give child-process cwd a physical path while retaining logical __file__."""

    def __init__(self, store):
        self.store = store

    def __getattr__(self, name):
        return getattr(self.store, name)

    def actual(self, path, view='legacy'):
        return self.store.merge_actual(self.store.actual(path, view))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--script')
    group.add_argument('--verify-layout', action='store_true')
    group.add_argument('--verify-merge', action='store_true')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = MergedLayout()
    runtime = MergedRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    # The inherited router calls runtime.spec, retaining this loader for runpy
    # and subprocess replays, including nested old-stage constructors.
    previous.previous.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), previous.previous.execution_routes(ExecutionLayout(store), runtime):
        if args.verify_layout or args.verify_merge:
            print(json.dumps(store.verify_layout() if args.verify_layout else store.verify_merge(),
                             ensure_ascii=False))
        else:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)


if __name__ == '__main__':
    main()
