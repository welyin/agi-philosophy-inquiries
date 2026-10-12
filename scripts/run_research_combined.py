"""Read-only replay through the 1120--1158 archive merge.

The outer manifest reverses navigation edits and maps real files back to each
module's pre-merge path. Older replay layers stay frozen. No directory copy,
ZIP, global pathlib replacement, or substitute scientific evidence is used.
"""
from __future__ import annotations

import argparse
import builtins
import importlib.util
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import sys
import types

import run_research_field_active as previous
from research_layout import (
    ROOT, Layout, ResearchLoader, _glob_match,
    no_archive_or_research_writes, relocation_original_bytes,
)

merged = previous.previous
THIS = Path(__file__).resolve()
BASE = 'research_cognition_physics'
NEW = BASE + '/archive_1120_1158'
STAGES = {
    BASE + '/archive_1120_': (1120, 1156),
    BASE + '/archive_1157_1158': (1157, 1158),
}
MANIFEST = ROOT / BASE / '_migration/merge_1120_1158/manifest.json'
digest, relative, safe_relative = merged.digest, merged.relative, merged.safe_relative


def moved_stage(name):
    return next((old for old in STAGES if name == old or name.startswith(old + '/')), None)


def stage_destination(name):
    old = moved_stage(name)
    if old is None:
        return name
    tail = name[len(old):].lstrip('/')
    group = '_'.join(map(str, STAGES[old]))
    if tail in ('README.md', '文件索引.md'):
        return NEW + '/_shared/' + group + '/' + tail
    if tail == '_shared' or tail.startswith('_shared/'):
        return NEW + '/_shared/' + group + tail[len('_shared'):]
    return NEW + ('/' + tail if tail else '')


class CombinedLayout(previous.FieldLayout):
    def __init__(self):
        raw = MANIFEST.read_bytes()
        self.combined_manifest_sha256 = digest(raw)
        self.combined_manifest = json.loads(raw.decode('utf-8-sig'))
        if self.combined_manifest.get('schema') != 'research_merge_1120_1158_v1':
            raise ValueError('Unsupported 1120--1158 merge manifest')
        entries = self.combined_manifest.get('entries')
        if not isinstance(entries, list):
            raise ValueError('Combined entries must be a list')
        self.combined_entries, self.combined_reverse = {}, {}
        sources, destinations = set(), set()
        for entry in entries:
            old, new = safe_relative(entry['original']), safe_relative(entry['destination'])
            if moved_stage(old):
                if new != stage_destination(old):
                    raise ValueError('Unexpected combined destination: ' + new)
            elif old != new or not old.endswith('.md'):
                raise ValueError('External changes must be in-place Markdown: ' + old)
            old_key, new_key = os.path.normcase(old), os.path.normcase(new)
            if old_key in sources or new_key in destinations:
                raise ValueError('Duplicate combined source or destination: ' + old)
            sources.add(old_key)
            destinations.add(new_key)
            for field in ('original_sha256', 'current_sha256'):
                if not isinstance(entry.get(field), str) or not re.fullmatch(r'[0-9a-f]{64}', entry[field]):
                    raise ValueError('Invalid combined hash: ' + old)
            for field in ('original_bytes', 'current_bytes'):
                if type(entry.get(field)) is not int or entry[field] < 0:
                    raise ValueError('Invalid combined size: ' + old)
            if not isinstance(entry.get('edits'), list):
                raise ValueError('Combined edits must be a list: ' + old)
            for edit in entry['edits']:
                if (type(edit.get('start')) is not int or type(edit.get('end')) is not int
                        or edit['start'] < 0 or edit['end'] < edit['start']
                        or not isinstance(edit.get('before'), str)
                        or not isinstance(edit.get('after'), str)):
                    raise ValueError('Invalid combined edit: ' + old)
            if not new.endswith('.md') and (entry['edits']
                    or entry['original_sha256'] != entry['current_sha256']
                    or entry['original_bytes'] != entry['current_bytes']):
                raise ValueError('Non-Markdown assets must stay byte-identical: ' + old)
            if Path(new).suffix.lower() == '.zip':
                raise ValueError('ZIPs are not replay assets: ' + new)
            if 'optional_local_cache' in entry:
                if type(entry['optional_local_cache']) is not bool:
                    raise ValueError('Optional-cache marker must be Boolean: ' + old)
                if entry['optional_local_cache'] and (
                        Path(old).suffix.lower() not in ('.pyc', '.pyo')
                        or '__pycache__' not in PurePosixPath(old).parts):
                    raise ValueError('Only local bytecode caches may be optional: ' + old)
            # A pre-existing local bytecode cache may be moved intact. It is
            # inventoried only; script lookup and loaders accept source .py.
            self.combined_entries[old] = entry
            self.combined_reverse[new] = old
        self.combined_directories = {
            parent.as_posix() for name in self.combined_entries
            for parent in PurePosixPath(name).parents
        }
        created = self.combined_manifest.get('created_files', {})
        self.combined_created_sizes = {}
        if isinstance(created, list):
            self.combined_created = {e['path']: e['sha256'] for e in created}
            self.combined_created_sizes = {e['path']: e['bytes'] for e in created}
            if len(self.combined_created) != len(created):
                raise ValueError('Duplicate combined navigation file')
        else:
            self.combined_created = created
        if not isinstance(self.combined_created, dict) or set(self.combined_created) != {
                NEW + '/README.md', NEW + '/文件索引.md'}:
            raise ValueError('Exactly the combined README and index must be new')
        super().__init__()
        self.current_files.update(self.combined_entries)
        self.current_directories.update(self.combined_directories)

    def combined_original_name(self, path):
        name = relative(path)
        if name in self.combined_reverse:
            return self.combined_reverse[name]
        if name is not None:
            for old, bounds in STAGES.items():
                prefix = NEW + '/_shared/' + '_'.join(map(str, bounds))
                if name == prefix or name.startswith(prefix + '/'):
                    return old + '/_shared' + name[len(prefix):]
            if name.startswith(NEW + '/'):
                first = name[len(NEW) + 1:].split('/', 1)[0]
                if first.isdigit():
                    for old, (start, end) in STAGES.items():
                        if start <= int(first) <= end:
                            return old + name[len(NEW):]
        return name

    def combined_name_path(self, path):
        name = self.combined_original_name(path)
        return ROOT / name if name is not None else Path(str(path))

    def combined_actual(self, path):
        old = self.combined_original_name(path)
        if old in self.combined_entries:
            return ROOT / self.combined_entries[old]['destination']
        if old and moved_stage(old):
            return ROOT / stage_destination(old)
        return Path(str(path))

    def before_combined(self, path):
        old = self.combined_original_name(path)
        actual = self.combined_actual(path)
        if old and moved_stage(old) and old not in self.combined_entries:
            raise FileNotFoundError('Historical file missing from combined manifest: ' + old)
        raw = actual.read_bytes()
        entry = self.combined_entries.get(old)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def before_field(self, path):
        # The single physical-read boundary for all earlier inverse layers.
        raw = self.before_combined(path)
        entry = self.field_entries.get(relative(self.combined_name_path(path)))
        return relocation_original_bytes(raw, entry) if entry else raw

    def key(self, path, view='legacy'):
        logical = self.combined_name_path(path)
        return relative(logical) if view == 'combined' else super().key(logical, view)

    def actual(self, path, view='legacy'):
        if view == 'combined':
            return self.combined_actual(path)
        return super().actual(self.combined_name_path(path), view)

    def read(self, path, historical=True, view='legacy'):
        old = self.combined_original_name(path)
        entry = self.combined_entries.get(old)
        if entry and entry.get('optional_local_cache'):
            raise PermissionError('A local bytecode cache is not a research input: ' + old)
        if view == 'combined':
            return self.before_combined(path) if historical else self.combined_actual(path).read_bytes()
        return super().read(self.combined_name_path(path), historical=historical, view=view)

    def combined_glob(self, directory, pattern):
        name = self.combined_original_name(directory)
        prefix = name.rstrip('/') + '/' if name not in (None, '.') else ''
        found = {}
        for candidate in set(self.combined_entries) | self.combined_directories:
            if candidate.startswith(prefix) and _glob_match(candidate[len(prefix):], pattern):
                found[candidate] = ROOT / candidate
        native = Path(str(self.combined_name_path(directory)))
        if native.is_dir():
            for actual in native.glob(pattern):
                physical = relative(actual)
                if physical == NEW or (physical and physical.startswith(NEW + '/')):
                    continue
                logical = self.combined_original_name(actual)
                if logical and logical.startswith(prefix) and _glob_match(logical[len(prefix):], pattern):
                    found[logical] = ROOT / logical
                elif name is None:
                    found[str(actual)] = actual
        yield from (found[key] for key in sorted(found))

    def glob_paths(self, directory, pattern, view):
        if view == 'combined' or moved_stage(self.combined_original_name(directory) or ''):
            yield from self.combined_glob(directory, pattern)
        else:
            yield from super().glob_paths(self.combined_name_path(directory), pattern, view)

    def path_class(self, view='legacy'):
        return Layout.path_class(self, view) if view == 'combined' else super().path_class(view)

    def script_key(self, path):
        old = self.combined_original_name(path)
        if old in self.combined_entries and moved_stage(old) and old.endswith('.py'):
            return old
        return super().script_key(self.combined_name_path(path))

    def script_actual(self, key):
        if key in self.combined_entries and moved_stage(key):
            return ROOT / self.combined_entries[key]['destination']
        return self.combined_actual(super().script_actual(key))

    def lookup_script(self, name):
        for candidate in (ROOT / name, ROOT / BASE / name):
            key = self.script_key(candidate)
            if key:
                return key
        return super().lookup_script(name)

    def verify_combined(self):
        if digest(MANIFEST.read_bytes()) != self.combined_manifest_sha256:
            raise ValueError('Combined manifest changed during replay')
        expected = set(self.combined_created)
        moved = scientific = markdown = caches = absent_optional_caches = 0
        for old, entry in self.combined_entries.items():
            moved += old != entry['destination']
            caches += Path(old).suffix.lower() in ('.pyc', '.pyo')
            if old != entry['destination'] and (ROOT / old).exists():
                raise ValueError('Old file remains after combined relocation: ' + old)
            if (entry.get('optional_local_cache')
                    and not (ROOT / entry['destination']).exists()):
                absent_optional_caches += 1
                continue
            prior = self.before_combined(ROOT / old)
            markdown += entry['original_sha256'] != entry['current_sha256']
            if Path(old).suffix.lower() in ('.py', '.json'):
                if digest(prior) != entry['current_sha256']:
                    raise ValueError('Scientific bytes changed: ' + old)
                scientific += 1
            if moved_stage(old):
                expected.add(entry['destination'])
        for name, expected_hash in self.combined_created.items():
            raw = (ROOT / safe_relative(name)).read_bytes()
            if digest(raw) != expected_hash or (name in self.combined_created_sizes
                    and len(raw) != self.combined_created_sizes[name]):
                raise ValueError('Combined navigation file changed: ' + name)
        inventory = {p.relative_to(ROOT).as_posix() for p in (ROOT / NEW).rglob('*') if p.is_file()}
        if inventory != expected or any((ROOT / old).exists() for old in STAGES):
            raise ValueError('Combined archive inventory differs from manifest')
        support = self.combined_manifest.get('unchanged_support_sha256')
        if not isinstance(support, dict) or not support:
            raise ValueError('Frozen support hashes are required')
        for name, expected_hash in support.items():
            if digest((ROOT / safe_relative(name)).read_bytes()) != expected_hash:
                raise ValueError('Frozen support changed: ' + name)
        return dict(files=len(self.combined_entries), moved_files=moved,
                    scientific_python_json_byte_identical=scientific,
                    reversible_markdown_files=markdown, relocated_cache_files=caches,
                    absent_optional_local_caches=absent_optional_caches,
                    cache_files_executed=0, new_navigation_files=len(self.combined_created),
                    frozen_support_files=len(support), manifest_sha256=self.combined_manifest_sha256,
                    current_and_original_hashes_verified=True, old_directories_absent=True,
                    directory_copy_used=False, zip_used=False)

    def verify_layout(self):
        combined = self.verify_combined()
        frozen = super().verify_layout()
        return dict(combined_layout=combined, frozen_layers=frozen, all_checks_passed=True)


class CombinedLoader(ResearchLoader):
    def exec_module(self, module):
        if not moved_stage(self.key):
            return merged.MergedLoader(self.runtime, self.key).exec_module(module)
        store = self.runtime.store
        logical, actual = ROOT / self.key, store.script_actual(self.key)
        module.__file__ = str(logical)
        module.__builtins__ = self.runtime.scoped_builtins['combined']
        old_path = list(sys.path)
        sys.path.insert(0, str(logical.parent))
        try:
            exec(compile(store.read(logical, view='combined'), str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = old_path


class CombinedRuntime(merged.MergedRuntime):
    def __init__(self, store):
        super().__init__(store)
        import pathlib
        proxy = types.ModuleType('pathlib')
        proxy.__dict__.update(vars(pathlib))
        proxy.Path = store.path_class('combined')
        default_import, default_open = builtins.__import__, builtins.open

        def scoped_import(name, globals=None, locals=None, fromlist=(), level=0):
            return proxy if name == 'pathlib' and level == 0 else default_import(
                name, globals, locals, fromlist, level)

        def scoped_open(file, mode='r', buffering=-1, encoding=None, errors=None,
                        newline=None, closefd=True, opener=None):
            if isinstance(file, (str, bytes, os.PathLike)):
                name = os.fsdecode(file)
                key = store.key(name, 'combined')
                if key and (key.startswith(BASE + '/') or key in store.combined_entries):
                    if opener is not None or not closefd:
                        raise NotImplementedError('Custom historical openers are unsupported')
                    return proxy.Path(name).open(mode, buffering, encoding, errors, newline)
            return default_open(file, mode, buffering, encoding, errors, newline, closefd, opener)
        self.scoped_builtins['combined'] = dict(vars(builtins), __import__=scoped_import, open=scoped_open)

    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, CombinedLoader(self, key),
                                              origin=str(self.store.script_actual(key)))


class ExecutionLayout(merged.ExecutionLayout):
    def actual(self, path, view='legacy'):
        return self.store.combined_actual(super().actual(path, view))


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
    store = CombinedLayout()
    runtime = CombinedRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    routes = merged.previous.previous
    routes.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), routes.execution_routes(ExecutionLayout(store), runtime):
        if args.verify_layout or args.verify_merge:
            print(json.dumps(store.verify_layout() if args.verify_layout else store.verify_combined(),
                             ensure_ascii=False))
        else:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)


if __name__ == '__main__':
    main()
