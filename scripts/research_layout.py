"""Resolve archived research dependencies to their current, physical locations.

There is no archive extraction or shadow directory. Historical Markdown hashes
can be checked by reversing the recorded link-only edits in memory; every read
first verifies the current file. Scientific Python and JSON bytes are unchanged.
"""
from __future__ import annotations

import ast
import builtins
from contextlib import contextmanager
import fnmatch
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
MIG = BASE / '_migration'
SPLIT_MANIFEST = 'layout_764_1008_20261008/manifest.json'
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


def relocation_original_bytes(raw, entry):
    """Verify and reverse the newest recorded UTF-8 text edits in memory."""
    if sha(raw) != entry['current_sha256'] or len(raw) != entry['current_bytes']:
        raise ValueError('Relocated file changed: ' + entry['destination'])
    edits = entry['edits']
    if edits or entry['original_sha256'] != entry['current_sha256']:
        if not entry['destination'].endswith('.md'):
            raise ValueError('Non-Markdown migration edit: ' + entry['destination'])
        enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
        text, pieces, offset, previous, original_end = raw.decode(enc), [], 0, 0, 0
        for edit in edits:
            if edit['start'] < original_end or edit['end'] - edit['start'] != len(edit['before']):
                raise ValueError('Invalid relocation edit interval: ' + entry['destination'])
            start = edit['start'] + offset
            end = start + len(edit['after'])
            if start < previous or text[start:end] != edit['after']:
                raise ValueError('Relocation edit no longer matches: ' + entry['destination'])
            pieces.extend((text[previous:start], edit['before']))
            previous, original_end = end, edit['end']
            offset += len(edit['after']) - len(edit['before'])
        pieces.append(text[previous:])
        restored = ''.join(pieces).encode(enc)
    else:
        restored = raw
    if sha(restored) != entry['original_sha256'] or len(restored) != entry['original_bytes']:
        raise ValueError('Pre-split hash not recovered: ' + entry['destination'])
    return restored


def _glob_match(name, pattern):
    """Path-segment glob matching, including zero components for **."""
    parts, rules = tuple(PurePosixPath(name).parts), tuple(PurePosixPath(pattern).parts)
    memo = {}

    def match(i, j):
        if (i, j) not in memo:
            if j == len(rules):
                result = i == len(parts)
            elif rules[j] == '**':
                result = match(i, j + 1) or (i < len(parts) and match(i + 1, j))
            else:
                result = i < len(parts) and fnmatch.fnmatchcase(parts[i], rules[j]) and match(i + 1, j + 1)
            memo[i, j] = result
        return memo[i, j]

    return match(0, 0)


class Layout:
    def __init__(self, root=ROOT):
        self.root = Path(root).resolve()
        self.base = self.root / 'research_cognition_physics'
        central = self.base / '_migration'
        self.migration = central if (central / 'manifest.json').is_file() else self.base / 'archive_764_/_migration'
        self.entries, self.reverse, self.directories = {}, {}, set()
        self.relocations, self.relocation_reverse = {}, {}
        self.current_directories, self.current_files = set(), set()
        self.legacy_destinations, self.legacy_aliases = {}, {}
        self.directory_aliases = {'legacy': {}, 'current': {}}
        self.read_paths, self.restored_markdown = set(), set()
        self.relocation_manifest = central / SPLIT_MANIFEST
        if self.relocation_manifest.is_file():
            data = json.loads(self.relocation_manifest.read_text('utf8'))
            if data['version'] != 1:
                raise ValueError('Unsupported relocation manifest version')
            for entry in data['entries']:
                logical = PurePosixPath(entry['original'].replace('\\', '/')).as_posix()
                actual = self.root / entry['destination']
                if not actual.resolve().is_relative_to(self.root) or logical.startswith('../'):
                    raise ValueError('Relocation outside project: ' + str(actual))
                if logical in self.relocations:
                    raise ValueError('Duplicate relocation name: ' + logical)
                physical = os.path.normcase(os.path.abspath(str(actual)))
                if physical in self.relocation_reverse:
                    raise ValueError('Duplicate relocation destination: ' + str(actual))
                self.relocations[logical] = (actual, entry)
                self.relocation_reverse[physical] = logical
                self.current_files.add(logical)
        for source, old_prefix, new_base in (
            ('manifest.json', LATE + '/', self.base),
            ('layout_001_230_20261004/manifest.json', '', self.root),
        ):
            data = json.loads(self.read_pre_split(self.migration / source))
            for entry in data['entries']:
                logical = old_prefix + entry['original']
                prior = new_base / entry['destination']
                prior_key = prior.relative_to(self.root).as_posix()
                actual = self.relocations.get(prior_key, (prior, None))[0]
                if not actual.resolve().is_relative_to(self.root):
                    raise ValueError('Destination outside project: ' + str(actual))
                self.entries[logical] = (actual, entry)
                self.legacy_destinations[logical] = prior
                self.legacy_aliases[os.path.normcase(os.path.abspath(str(prior)))] = logical
                self.reverse.setdefault(os.path.normcase(str(actual)), logical)
                self.current_files.add(prior_key)
                for parent in PurePosixPath(logical).parents:
                    self.directories.add(str(parent))
        for logical in self.current_files:
            self.current_directories.update(str(parent) for parent in PurePosixPath(logical).parents)
        # A moved round has one old directory; a split stage may span several
        # physical directories. These aliases never require an on-disk shell.
        for view, items in (('legacy', self.entries), ('current', self.relocations)):
            candidates = {}
            for logical, (actual, _) in items.items():
                for old_parent, new_parent in zip((self.root / logical).parents, actual.parents):
                    if (not old_parent.is_relative_to(self.base) or old_parent == self.base
                            or not new_parent.is_relative_to(self.base) or new_parent == self.base):
                        break
                    candidates.setdefault(os.path.normcase(str(new_parent)), set()).add(str(old_parent))
            self.directory_aliases[view] = {
                physical: next(iter(names)) for physical, names in candidates.items() if len(names) == 1}

    def key(self, path, view='legacy'):
        if view not in ('legacy', 'current'):
            raise ValueError('Unknown historical view: ' + view)
        raw = os.path.abspath(str(path))
        normalized = os.path.normcase(raw)
        try:
            direct = Path(raw).relative_to(self.root).as_posix()
        except ValueError:
            return None
        if view == 'legacy':
            if direct in self.entries:
                return direct
            if normalized in self.reverse:
                return self.reverse[normalized]
            if normalized in self.legacy_aliases:
                return self.legacy_aliases[normalized]
        if direct in self.relocations:
            return direct
        if normalized in self.relocation_reverse:
            return self.relocation_reverse[normalized]
        alias = self.directory_aliases[view].get(normalized)
        return Path(alias).relative_to(self.root).as_posix() if alias else direct

    def actual(self, path, view='legacy'):
        key = self.key(path, view)
        if key in self.entries:
            return self.entries[key][0]
        return self.relocations[key][0] if key in self.relocations else Path(str(path))

    def read_pre_split(self, path):
        """Read the real file, verifying the new layer before returning its prior bytes."""
        key = self.key(path, 'current')
        if key not in self.relocations:
            return Path(str(path)).read_bytes()
        actual, entry = self.relocations[key]
        raw = actual.read_bytes()
        original = relocation_original_bytes(raw, entry)
        self.read_paths.add(actual.relative_to(self.root).as_posix())
        if raw != original:
            self.restored_markdown.add(actual.relative_to(self.root).as_posix())
        return original

    def read(self, path, historical=True, view='legacy'):
        key = self.key(path, view)
        if local_backup(str(path)) and not str(path).endswith('talk_20260629_kikukawa.pdf'):
            raise FileNotFoundError('Local backup is not a runtime dependency: ' + str(path))
        if key not in self.entries:
            prior = self.read_pre_split(path)
            return prior if historical else self.actual(path, 'current').read_bytes()
        actual, entry = self.entries[key]
        prior = self.read_pre_split(self.legacy_destinations[key])
        original = original_bytes(prior, entry)
        raw = actual.read_bytes()
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
        for candidate in (self.root / name, self.base / name):
            key = self.script_key(candidate)
            if key:
                return key
        matches = [key for key in self.relocations if key.endswith('/' + name) and key.endswith('.py')]
        if len(matches) == 1:
            return matches[0]
        raise ValueError('Unknown research script: ' + name)

    def script_key(self, path):
        legacy = self.key(path, 'legacy')
        if legacy in self.entries and legacy.endswith('.py'):
            return legacy
        current = self.key(path, 'current')
        if current in self.relocations and current.endswith('.py'):
            return current
        return None

    def script_actual(self, key):
        return self.entries[key][0] if key in self.entries else self.relocations[key][0]

    def verify(self, include_local=False):
        unchanged = links = skipped = 0
        for key, (path, entry) in self.entries.items():
            if not include_local and local_backup(key):
                skipped += 1
                continue
            raw = self.read_pre_split(self.legacy_destinations[key])
            original = original_bytes(raw, entry)
            if raw == original:
                unchanged += 1
            else:
                links += 1
        return dict(mapped_files=len(self.entries), checked_files=unchanged + links,
                    unchanged_files=unchanged, verified_link_only_changes=links,
                    optional_local_files_skipped=skipped, zip_files_read=0)

    def verify_relocation(self):
        """New migration statistics, separate from the frozen 001--775 statistics."""
        identical = rewritten = edits = 0
        kinds = {}
        for actual, entry in self.relocations.values():
            if actual.suffix.lower() == '.zip':
                raise ValueError('ZIP is not a relocation runtime dependency: ' + str(actual))
            raw = actual.read_bytes()
            prior = relocation_original_bytes(raw, entry)
            identical += raw == prior
            rewritten += raw != prior
            edits += len(entry['edits'])
            for edit in entry['edits']:
                kinds[edit['kind']] = kinds.get(edit['kind'], 0) + 1
        return dict(mapped_files=len(self.relocations), checked_files=identical + rewritten,
                    byte_identical_files=identical, rewritten_markdown_files=rewritten,
                    recorded_text_edits=edits, edit_kinds=kinds,
                    current_and_pre_split_hashes_verified=True, zip_files_read=0)

    def glob_paths(self, directory, pattern, view):
        key = self.key(directory, view)
        if view == 'current' and key in self.directories and key not in self.current_directories:
            view = 'legacy'
        prefix = key.rstrip('/') + '/' if key not in (None, '.') else ''
        files = set(self.entries) if view == 'legacy' else self.current_files
        dirs = self.directories if view == 'legacy' else self.current_directories
        found = {}
        for logical in files | dirs:
            if logical.startswith(prefix) and _glob_match(logical[len(prefix):], pattern):
                found[logical] = self.root / logical
        native = Path(str(directory))
        if native.is_dir():
            for actual in native.glob(pattern):
                logical = self.key(actual, view)
                if logical and logical.startswith(prefix) and _glob_match(logical[len(prefix):], pattern):
                    found[logical] = self.root / logical
                elif key is None:
                    found[str(actual)] = actual
        yield from (found[name] for name in sorted(found))

    def path_class(self, view='legacy'):
        store = self
        native = type(Path())

        class ResearchPath(native):
            """Read-only legacy name that resolves to a current physical file."""
            def __new__(cls, *parts, **kwargs):
                return super().__new__(cls, *parts, **kwargs)

            def __init__(self, *parts):
                raw = native(*parts)
                key = store.key(raw, view)
                super().__init__(store.root / key if raw.is_absolute() and key is not None else raw)

            def resolve(self, strict=False):
                result = type(self)(os.path.abspath(str(self)))
                if strict and not result.exists():
                    raise FileNotFoundError(str(self))
                return result

            def exists(self, **kwargs):
                key = store.key(self, view)
                dirs = store.directories if view == 'legacy' else store.current_directories
                return key in dirs or key in store.directories or store.actual(self, view).exists()

            def is_file(self, **kwargs):
                return store.actual(self, view).is_file()

            def is_dir(self, **kwargs):
                dirs = store.directories if view == 'legacy' else store.current_directories
                key = store.key(self, view)
                return key in dirs or key in store.directories or Path(str(self)).is_dir()

            def read_bytes(self):
                return store.read(self, view=view)

            def read_text(self, encoding=None, errors=None, newline=None):
                with self.open('r', encoding=encoding, errors=errors, newline=newline) as stream:
                    return stream.read()

            def open(self, mode='r', buffering=-1, encoding=None, errors=None, newline=None):
                if any(c in mode for c in 'wax+'):
                    raise PermissionError('Historical replay is read-only: ' + str(self))
                stream = io.BytesIO(self.read_bytes())
                return stream if 'b' in mode else io.TextIOWrapper(
                    stream, encoding=encoding or 'utf8', errors=errors, newline=newline)

            def glob(self, pattern, **kwargs):
                if kwargs:
                    raise NotImplementedError('Historical glob keyword options are not supported')
                yield from (type(self)(p) for p in store.glob_paths(self, pattern, view))

            def rglob(self, pattern, **kwargs):
                yield from self.glob('**/' + pattern, **kwargs)

            def iterdir(self):
                yield from self.glob('*')

        return ResearchPath


class ResearchLoader(importlib.abc.Loader):
    def __init__(self, runtime, key):
        self.runtime, self.key = runtime, key

    def create_module(self, spec):
        return None

    def exec_module(self, module):
        store = self.runtime.store
        actual = store.script_actual(self.key)
        view = 'legacy' if self.key in store.entries else 'current'
        logical = store.root / self.key
        raw = store.read(logical, view=view)
        module.__file__ = str(logical)
        module.__builtins__ = self.runtime.scoped_builtins[view]
        previous_path = list(sys.path)
        sys.path.insert(0, str(logical.parent))
        try:
            exec(compile(raw, str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = previous_path


class ResearchRuntime(importlib.abc.MetaPathFinder):
    def __init__(self, store, scope='late'):
        self.store = store
        self.path = store.path_class()
        self.search = ([EARLY, SECOND, LATE] if scope == 'early' else [LATE, SECOND, EARLY])
        import pathlib
        default_import = builtins.__import__
        default_open = builtins.open
        self.scoped_builtins = {}

        def make_builtins(view):
            proxy = types.ModuleType('pathlib')
            proxy.__dict__.update(vars(pathlib))
            proxy.Path = store.path_class(view)

            def scoped_import(name, globals=None, locals=None, fromlist=(), level=0):
                if name == 'pathlib' and level == 0:
                    return proxy
                return default_import(name, globals, locals, fromlist, level)

            def scoped_open(file, mode='r', buffering=-1, encoding=None, errors=None,
                            newline=None, closefd=True, opener=None):
                if isinstance(file, (str, bytes, os.PathLike)):
                    path = os.fsdecode(file)
                    key = store.key(path, view)
                    if key and (key.startswith('research_cognition_physics/')
                                or key in store.entries or key in store.relocations):
                        if opener is not None or not closefd:
                            raise NotImplementedError('Custom historical file openers are not supported')
                        return proxy.Path(path).open(mode, buffering, encoding, errors, newline)
                return default_open(file, mode, buffering, encoding, errors, newline, closefd, opener)

            return dict(vars(builtins), __import__=scoped_import, open=scoped_open)

        self.scoped_builtins = {view: make_builtins(view) for view in ('legacy', 'current')}
        self.builtins = self.scoped_builtins['legacy']

    def spec(self, name, key):
        actual = self.store.script_actual(key)
        return importlib.util.spec_from_loader(name, ResearchLoader(self, key), origin=str(actual))

    def find_spec(self, fullname, path=None, target=None):
        name = fullname.rsplit('.', 1)[-1]
        directories = list(path) if path is not None else [
            *sys.path, *(str(self.store.root / p) for p in self.search)]
        for directory in directories:
            key = self.store.script_key(Path(directory) / (name + '.py'))
            if key:
                return self.spec(fullname, key)
            package = self.store.key(Path(directory) / name, 'current')
            if package in self.store.directories or package in self.store.current_directories:
                spec = importlib.util.spec_from_loader(fullname, loader=None, is_package=True)
                spec.submodule_search_locations = [str(self.store.root / package)]
                return spec
        return None

    @contextmanager
    def installed(self):
        old_spec = importlib.util.spec_from_file_location

        def mapped_spec(name, location=None, *args, **kwargs):
            key = self.store.script_key(location) if location is not None else None
            if key:
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
        sys.argv = [str(self.store.root / key), *arguments]
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
