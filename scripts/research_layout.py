"""Resolve archived research dependencies to their current, physical locations.

There is no archive extraction or shadow directory. Historical Markdown hashes
can be checked by reversing the recorded link-only edits in memory; every read
first verifies the current file. Scientific Python and JSON bytes are unchanged.
"""
from __future__ import annotations

import ast
import builtins
from contextlib import contextmanager
import hashlib
import importlib.abc
import importlib.util
import io
import json
import os
from pathlib import Path, PurePosixPath
import sys
import types

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'research_cognition_physics'
MIG = BASE / 'archive_764_/_migration'
EARLY = 'research_cognition_physics/archive_001_222/research_process'
SECOND = 'research_cognition_physics/archive_223_230'
LATE = 'research_cognition_physics/archive_231_'


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def local_backup(name):
    p = PurePosixPath(name.replace('\\', '/'))
    return (any(t.startswith('navigation_before_') for t in p.parts)
            or 'runtime_cache' in p.parts or '__pycache__' in p.parts
            or p.suffix in ('.pyc', '.pyo', '.zip')
            or p.name == 'talk_20260629_kikukawa.pdf')


def original_bytes(raw, entry):
    """Prove both hashes; never substitute a stored digest for actual content."""
    if sha(raw) != entry['current_sha256']:
        raise ValueError('Current file changed: ' + entry['destination'])
    if entry['original_sha256'] == entry['current_sha256']:
        return raw
    if not entry['rewrite_markdown_links'] or not entry['destination'].endswith('.md'):
        raise ValueError('Unaccounted migration change: ' + entry['destination'])
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    text, pieces, offset, previous = raw.decode(enc), [], 0, 0
    for edit in entry['link_edits']:
        start = edit['start'] + offset
        end = start + len(edit['after'])
        if start < previous or text[start:end] != edit['after']:
            raise ValueError('Recorded link edit no longer matches: ' + entry['destination'])
        pieces.extend((text[previous:start], edit['before']))
        previous = end
        offset += len(edit['after']) - (edit['end'] - edit['start'])
    pieces.append(text[previous:])
    restored = ''.join(pieces).encode(enc)
    if sha(restored) != entry['original_sha256']:
        raise ValueError('Original hash not recovered: ' + entry['destination'])
    return restored


class Layout:
    def __init__(self, root=ROOT):
        self.root = Path(root).resolve()
        self.base = self.root / 'research_cognition_physics'
        self.migration = self.base / 'archive_764_/_migration'
        self.entries, self.reverse, self.directories = {}, {}, set()
        self.read_paths, self.restored_markdown = set(), set()
        for source, old_prefix, new_base in (
            ('manifest.json', LATE + '/', self.base),
            ('layout_001_230_20261004/manifest.json', '', self.root),
        ):
            data = json.loads((self.migration / source).read_text('utf8'))
            for entry in data['entries']:
                logical = old_prefix + entry['original']
                actual = new_base / entry['destination']
                if not actual.resolve().is_relative_to(self.root):
                    raise ValueError('Destination outside project: ' + str(actual))
                self.entries[logical] = (actual, entry)
                self.reverse.setdefault(os.path.normcase(str(actual)), logical)
                for parent in PurePosixPath(logical).parents:
                    self.directories.add(str(parent))

    def key(self, path):
        raw = os.path.abspath(str(path))
        if os.path.normcase(raw) in self.reverse:
            return self.reverse[os.path.normcase(raw)]
        try:
            return Path(raw).relative_to(self.root).as_posix()
        except ValueError:
            return None

    def actual(self, path):
        key = self.key(path)
        return self.entries[key][0] if key in self.entries else Path(str(path))

    def read(self, path, historical=True):
        key = self.key(path)
        if local_backup(str(path)) and not str(path).endswith('talk_20260629_kikukawa.pdf'):
            raise FileNotFoundError('Local backup is not a runtime dependency: ' + str(path))
        if key not in self.entries:
            return Path(str(path)).read_bytes()
        actual, entry = self.entries[key]
        raw = actual.read_bytes()
        original = original_bytes(raw, entry)
        self.read_paths.add(actual.relative_to(self.root).as_posix())
        if raw != original:
            self.restored_markdown.add(actual.relative_to(self.root).as_posix())
        return original if historical else raw

    def lookup_script(self, name):
        candidates = [name, LATE + '/' + name, EARLY + '/' + name, SECOND + '/' + name]
        current = self.key(self.root / name)
        if current:
            candidates.insert(0, current)
        current = self.key(self.base / name)
        if current:
            candidates.insert(0, current)
        for key in candidates:
            if key in self.entries and key.endswith('.py'):
                return key
        raise ValueError('Unknown research script: ' + name)

    def verify(self, include_local=False):
        unchanged = links = skipped = 0
        for key, (path, entry) in self.entries.items():
            if not include_local and local_backup(key):
                skipped += 1
                continue
            raw = path.read_bytes()
            original = original_bytes(raw, entry)
            if raw == original:
                unchanged += 1
            else:
                links += 1
        return dict(mapped_files=len(self.entries), checked_files=unchanged + links,
                    unchanged_files=unchanged, verified_link_only_changes=links,
                    optional_local_files_skipped=skipped, zip_files_read=0)

    def path_class(self):
        store = self
        native = type(Path())

        class ResearchPath(native):
            """Read-only legacy name that resolves to a current physical file."""
            def __new__(cls, *parts, **kwargs):
                return super().__new__(cls, *parts, **kwargs)

            def __init__(self, *parts):
                raw = native(*parts)
                key = store.reverse.get(os.path.normcase(os.path.abspath(str(raw))))
                super().__init__(store.root / key if key else raw)

            def resolve(self, strict=False):
                result = type(self)(os.path.abspath(str(self)))
                if strict and not result.exists():
                    raise FileNotFoundError(str(self))
                return result

            def exists(self, **kwargs):
                key = store.key(self)
                if key in store.entries:
                    return store.entries[key][0].is_file()
                return key in store.directories or Path(str(self)).exists()

            def is_file(self, **kwargs):
                return store.actual(self).is_file()

            def is_dir(self, **kwargs):
                return store.key(self) in store.directories or Path(str(self)).is_dir()

            def read_bytes(self):
                return store.read(self)

            def open(self, mode='r', buffering=-1, encoding=None, errors=None, newline=None):
                if any(c in mode for c in 'wax+'):
                    raise PermissionError('Historical replay is read-only: ' + str(self))
                stream = io.BytesIO(self.read_bytes())
                return stream if 'b' in mode else io.TextIOWrapper(
                    stream, encoding=encoding or 'utf8', errors=errors, newline=newline)

            def glob(self, pattern, **kwargs):
                key = store.key(self)
                if key not in store.directories:
                    yield from (type(self)(p) for p in Path(str(self)).glob(pattern))
                    return
                prefix = key.rstrip('/') + '/'
                for logical in sorted(store.entries):
                    if logical.startswith(prefix):
                        rest = logical[len(prefix):]
                        if PurePosixPath(rest).match(pattern) and ('/' in pattern or '/' not in rest):
                            yield type(self)(store.root / logical)

            def rglob(self, pattern, **kwargs):
                yield from self.glob('**/' + pattern)

        return ResearchPath


class ResearchLoader(importlib.abc.Loader):
    def __init__(self, runtime, key):
        self.runtime, self.key = runtime, key

    def create_module(self, spec):
        return None

    def exec_module(self, module):
        store = self.runtime.store
        actual, entry = store.entries[self.key]
        raw = store.read(store.root / self.key)
        module.__file__ = str(actual)
        module.__builtins__ = self.runtime.builtins
        exec(compile(raw, str(actual), 'exec'), module.__dict__)


class ResearchRuntime(importlib.abc.MetaPathFinder):
    def __init__(self, store, scope='late'):
        self.store = store
        self.path = store.path_class()
        self.search = ([EARLY, SECOND, LATE] if scope == 'early' else [LATE, SECOND, EARLY])
        proxy = types.ModuleType('pathlib')
        import pathlib
        proxy.__dict__.update(vars(pathlib))
        proxy.Path = self.path
        default_import = builtins.__import__

        def scoped_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == 'pathlib' and level == 0:
                return proxy
            return default_import(name, globals, locals, fromlist, level)

        self.builtins = dict(vars(builtins), __import__=scoped_import)

    def spec(self, name, key):
        actual = self.store.entries[key][0]
        return importlib.util.spec_from_loader(name, ResearchLoader(self, key), origin=str(actual))

    def find_spec(self, fullname, path=None, target=None):
        name = fullname.rsplit('.', 1)[-1]
        directories = list(path) if path is not None else [
            *sys.path, *(str(self.store.root / p) for p in self.search)]
        for directory in directories:
            key = self.store.key(Path(directory) / (name + '.py'))
            if key in self.store.entries:
                return self.spec(fullname, key)
            package = self.store.key(Path(directory) / name)
            if package in self.store.directories:
                spec = importlib.util.spec_from_loader(fullname, loader=None, is_package=True)
                spec.submodule_search_locations = [str(self.store.root / package)]
                return spec
        return None

    @contextmanager
    def installed(self):
        old_spec = importlib.util.spec_from_file_location

        def mapped_spec(name, location=None, *args, **kwargs):
            key = self.store.key(location) if location is not None else None
            if key in self.store.entries and key.endswith('.py'):
                return self.spec(name, key)
            return old_spec(name, location, *args, **kwargs)

        sys.meta_path.insert(0, self)
        importlib.util.spec_from_file_location = mapped_spec
        try:
            yield self
        finally:
            sys.meta_path.remove(self)
            importlib.util.spec_from_file_location = old_spec

    def run_script(self, name, arguments=()):
        key = self.store.lookup_script(name)
        old_argv = sys.argv
        old_main = sys.modules.get('__main__')
        module = types.ModuleType('__main__')
        module.__spec__ = self.spec('__main__', key)
        sys.modules['__main__'] = module
        sys.argv = [str(self.store.entries[key][0]), *arguments]
        try:
            ResearchLoader(self, key).exec_module(module)
        finally:
            sys.argv = old_argv
            sys.modules['__main__'] = old_main


def no_archive_or_research_writes(event, args):
    """Process-local guard, also covering reads that bypass the path adapter."""
    if event in ('os.remove', 'os.rmdir', 'os.mkdir', 'os.rename', 'os.link', 'os.symlink'):
        for item in args[:2]:
            if isinstance(item, (str, bytes, os.PathLike)):
                if Path(os.fsdecode(item)).absolute().is_relative_to(BASE):
                    raise PermissionError('Research evidence is read-only during replay')
    if event != 'open' or not isinstance(args[0], (str, bytes, os.PathLike)):
        return
    path = os.fsdecode(args[0])
    if path.lower().endswith('.zip'):
        raise PermissionError('ZIP reads are forbidden during current-layout replay')
    if any(p.startswith('navigation_before_') for p in Path(path).parts):
        raise PermissionError('Navigation backups are not runtime dependencies')
    mode, flags = args[1:3]
    writing = (isinstance(mode, str) and any(c in mode for c in 'wax+')) or (
        isinstance(flags, int) and flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC))
    if writing and Path(path).absolute().is_relative_to(BASE):
        raise PermissionError('Research evidence is read-only during replay: ' + path)


def test_modules(store, stage):
    prefix = EARLY if stage == '1' else SECOND
    result = []
    for key, (actual, entry) in store.entries.items():
        relative = key.removeprefix(prefix + '/')
        if key.startswith(prefix + '/') and '/' not in relative and key.endswith('.py'):
            tree = ast.parse(actual.read_text('utf-8-sig'))
            if any(isinstance(n, ast.ClassDef) and any(
                (isinstance(b, ast.Attribute) and b.attr == 'TestCase') or
                (isinstance(b, ast.Name) and b.id == 'TestCase') for b in n.bases
                ) for n in tree.body):
                result.append(PurePosixPath(key).stem)
    return sorted(result)
