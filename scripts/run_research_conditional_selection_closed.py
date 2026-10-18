"""Read-only replay of closed conditional-selection rounds 1308--1312.

The newest migration is reversed in memory before the unchanged active layer
and older closures are inspected. Python/JSON evidence remains byte-identical;
old and current paths resolve to one physical file, without shadow directories.
"""
from __future__ import annotations

import argparse
import builtins
from contextlib import contextmanager
import importlib.util
import json
import os
import re
from pathlib import Path, PurePosixPath
import subprocess
import sys
import types

import run_research_conditional_selection_active as previous
from research_layout import ROOT, Layout, ResearchLoader, _glob_match
from research_layout import no_archive_or_research_writes, relocation_original_bytes

structure = previous.previous
matter = structure.matter
THIS = Path(__file__).resolve()
BASE = 'research_cognition_physics'
STAGES = [
    dict(source=BASE + '/archive_1308_', destination=BASE + '/archive_1308_1312',
         first_round=1308, last_round=1312),
]
PREFIXES = {s['source']: s['destination'] for s in STAGES}
MANIFEST = ROOT / BASE / '_migration/closure_1308_1312/manifest.json'
VIEW = 'conditional_selection_closed'
digest, relative, inside = previous.digest, previous.relative, previous.inside
safe_relative = matter.safe_relative


def closed_stage(name):
    return any(inside(name, old) for old in PREFIXES)


def validate_entry(entry):
    # Preserve pre-existing caches as opaque inventory, never executable inputs.
    # All source execution below still requires a .py logical script key.
    old, new = safe_relative(entry['original']), safe_relative(entry['destination'])
    if Path(new).suffix.lower() not in ('.pyc', '.pyo'):
        return matter.gr.validate_entry(entry)
    if ('__pycache__' not in PurePosixPath(new).parts or not closed_stage(old)
            or entry.get('edits') != []
            or not re.fullmatch(r'[0-9a-f]{64}', entry.get('original_sha256', ''))
            or entry['current_sha256'] != entry['original_sha256']
            or type(entry.get('original_bytes')) is not int or entry['original_bytes'] < 0
            or entry['current_bytes'] != entry['original_bytes']):
        raise ValueError('Invalid opaque cache inventory: ' + new)
    return old, new


class ClosedConditionalSelectionLayout(previous.ConditionalSelectionLayout):
    def __init__(self):
        raw = MANIFEST.read_bytes()
        self.conditional_closure_sha256 = digest(raw)
        self.conditional_closure_manifest = json.loads(raw)
        if self.conditional_closure_manifest.get('schema') != 'research_closure_1308_1312_v1':
            raise ValueError('Unexpected conditional-selection closure schema')
        if self.conditional_closure_manifest.get('stages') != STAGES:
            raise ValueError('Unexpected structure-selection stage boundaries')
        self.conditional_closure_entries, self.conditional_closure_reverse = {}, {}
        old_names, new_names = set(), set()
        for entry in self.conditional_closure_manifest['entries']:
            old, new = validate_entry(entry)
            expected = next((dest + old[len(src):] for src, dest in PREFIXES.items()
                             if inside(old, src)), old)
            if new != expected or (not closed_stage(old) and not old.endswith('.md')):
                raise ValueError('Unexpected structure-selection destination: ' + old)
            if os.path.normcase(old) in old_names or os.path.normcase(new) in new_names:
                raise ValueError('Duplicate conditional-selection closure mapping')
            if Path(old).suffix.lower() in ('.py', '.json') and (
                    entry['edits'] or entry['original_sha256'] != entry['current_sha256']
                    or entry['original_bytes'] != entry['current_bytes']):
                raise ValueError('Scientific Python/JSON must remain byte-identical: ' + old)
            old_names.add(os.path.normcase(old))
            new_names.add(os.path.normcase(new))
            self.conditional_closure_entries[old] = entry
            self.conditional_closure_reverse[new] = old
        if not self.conditional_closure_entries:
            raise ValueError('Empty conditional-selection closure inventory')
        self.conditional_closure_directories = {
            p.as_posix() for name in self.conditional_closure_entries
            for p in PurePosixPath(name).parents}
        super().__init__()
        self.current_files.update(self.conditional_closure_entries)
        self.current_directories.update(self.conditional_closure_directories)
        for name, boundary in self.conditional_entries.items():
            entry = self.conditional_closure_entries.get(name)
            if entry is None or (entry['original_sha256'], entry['original_bytes']) != (
                    boundary['current_sha256'], boundary['current_bytes']):
                raise ValueError('Closure does not join active conditional navigation: ' + name)

    def conditional_closure_name(self, path):
        name = relative(path)
        if name in self.conditional_closure_reverse:
            return self.conditional_closure_reverse[name]
        for old, new in PREFIXES.items():
            if inside(name, new):
                return old + name[len(new):]
        return name

    def conditional_closure_logical(self, path):
        name = self.conditional_closure_name(path)
        return ROOT / name if name is not None else Path(str(path))

    def conditional_closure_actual(self, path):
        name = self.conditional_closure_name(path)
        if name in self.conditional_closure_entries:
            return ROOT / self.conditional_closure_entries[name]['destination']
        for old, new in PREFIXES.items():
            if inside(name, old):
                return ROOT / (new + name[len(old):])
        return Path(str(path))

    def before_conditional_closure(self, path):
        name = self.conditional_closure_name(path)
        if closed_stage(name) and name not in self.conditional_closure_entries:
            raise FileNotFoundError('File absent from conditional-selection closure: ' + str(name))
        raw = self.conditional_closure_actual(path).read_bytes()
        entry = self.conditional_closure_entries.get(name)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def before_conditional_selection_navigation(self, path):
        logical = self.conditional_closure_logical(path)
        name = self.structure_closure_name(logical)
        raw = self.before_conditional_closure(self.structure_closure_actual(logical))
        entry = self.conditional_entries.get(name)
        return relocation_original_bytes(raw, entry) if entry else raw

    def key(self, path, view='legacy'):
        logical = self.conditional_closure_logical(path)
        return relative(logical) if view == VIEW else super().key(logical, view)

    def actual(self, path, view='legacy'):
        if view == VIEW:
            return self.conditional_closure_actual(path)
        return self.conditional_closure_actual(super().actual(self.conditional_closure_logical(path), view))

    def read(self, path, historical=True, view='legacy'):
        if view == VIEW:
            return (self.before_conditional_closure(path) if historical
                    else self.conditional_closure_actual(path).read_bytes())
        return super().read(self.conditional_closure_logical(path), historical=historical, view=view)

    def glob_paths(self, directory, pattern, view):
        name = self.conditional_closure_name(directory)
        if view != VIEW and not closed_stage(name):
            yield from super().glob_paths(self.conditional_closure_logical(directory), pattern, view)
            return
        prefix = name.rstrip('/') + '/' if name not in (None, '.') else ''
        found = {candidate: ROOT / candidate
                 for candidate in set(self.conditional_closure_entries) | self.conditional_closure_directories
                 if candidate.startswith(prefix) and _glob_match(candidate[len(prefix):], pattern)}
        if not closed_stage(name):
            native = self.conditional_closure_actual(directory)
            if native.is_dir():
                for actual in native.glob(pattern):
                    logical = self.conditional_closure_name(actual)
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

        class ClosedStructurePath(base):
            def exists(self, **kwargs):
                return store.key(self, VIEW) in store.conditional_closure_directories or store.actual(self, VIEW).exists()

            def is_dir(self, **kwargs):
                return store.key(self, VIEW) in store.conditional_closure_directories or store.actual(self, VIEW).is_dir()

        return ClosedStructurePath

    def script_key(self, path):
        name = self.conditional_closure_name(path)
        if closed_stage(name):
            return name if name in self.conditional_closure_entries and name.endswith('.py') else None
        return super().script_key(self.conditional_closure_logical(path))

    def script_actual(self, key):
        if closed_stage(key):
            return self.conditional_closure_actual(ROOT / key)
        return self.conditional_closure_actual(super().script_actual(key))

    def verify_conditional_navigation(self):
        if digest(previous.LAYER.read_bytes()) != self.conditional_sha256:
            raise ValueError('Conditional-selection active layer changed during replay')
        for name, entry in self.conditional_entries.items():
            old = self.before_conditional_selection_navigation(ROOT / name)
            enc = 'utf-8-sig' if old.startswith(b'\xef\xbb\xbf') else 'utf8'
            if entry['edits'][0]['start'] != old.decode(enc).index('\n') + 1:
                raise ValueError('Navigation insertion must follow the original title: ' + name)
        support = self.conditional_layer['unchanged_support_sha256']
        required = dict(self.structure_closure_manifest['unchanged_support_sha256'])
        required[structure.MANIFEST.relative_to(ROOT).as_posix()] = self.structure_closure_sha256
        required[structure.THIS.relative_to(ROOT).as_posix()] = digest(structure.THIS.read_bytes())
        if not isinstance(support, dict) or not set(required).issubset(support):
            raise ValueError('Missing inherited support or closure bindings')
        for name, expected in support.items():
            if digest(self.before_conditional_closure(ROOT / safe_relative(name))) != expected:
                raise ValueError('Frozen conditional-navigation support changed: ' + name)
            if name in required and expected != required[name]:
                raise ValueError('Inherited support binding changed: ' + name)
        for name in previous.ENTRY_POINTS:
            self.before_conditional_closure(ROOT / name)
        return dict(status='PASS', files=4, manifest_sha256=self.conditional_sha256,
                    closed_boundary_verified=True, original_navigation_recoverable=True,
                    stage_first_round=1308, entry_points=previous.ENTRY_POINTS,
                    frozen_support_files=len(support), scientific_acceptance_checked=False)

    def verify_closure(self):
        if digest(MANIFEST.read_bytes()) != self.conditional_closure_sha256:
            raise ValueError('Conditional-selection closure changed during replay')
        expected = {new: set() for new in PREFIXES.values()}
        scientific = markdown = 0
        for name, entry in self.conditional_closure_entries.items():
            raw = self.before_conditional_closure(ROOT / name)
            for old, new in PREFIXES.items():
                if inside(name, old):
                    expected[new].add(entry['destination'])
            scientific += Path(name).suffix.lower() in ('.py', '.json')
            markdown += bool(entry['edits'])
            if Path(name).suffix.lower() in ('.py', '.json') and digest(raw) != entry['current_sha256']:
                raise ValueError('Scientific evidence changed: ' + name)
        for old, new in PREFIXES.items():
            actual = {p.relative_to(ROOT).as_posix() for p in (ROOT / new).rglob('*') if p.is_file()}
            if (ROOT / old).exists() or not (ROOT / new).is_dir() or actual != expected[new]:
                raise ValueError('Closed structure-selection inventory differs from manifest: ' + new)
        support = self.conditional_closure_manifest['unchanged_support_sha256']
        required = dict(self.conditional_layer['unchanged_support_sha256'])
        required[previous.LAYER.relative_to(ROOT).as_posix()] = self.conditional_sha256
        required[previous.THIS.relative_to(ROOT).as_posix()] = digest(previous.THIS.read_bytes())
        if not isinstance(support, dict) or not set(required).issubset(support):
            raise ValueError('Missing inherited frozen support bindings')
        for name, expected_hash in support.items():
            if digest((ROOT / safe_relative(name)).read_bytes()) != expected_hash:
                raise ValueError('Frozen support changed: ' + name)
            if name in required and expected_hash != required[name]:
                raise ValueError('Inherited support binding changed: ' + name)
        return dict(status='PASS', scope='1308--1312 directory closure only',
                    files=len(self.conditional_closure_entries), moved_files=sum(map(len, expected.values())),
                    scientific_python_json_byte_identical=scientific,
                    reversible_markdown_files=markdown, frozen_support_files=len(support),
                    manifest_sha256=self.conditional_closure_sha256,
                    old_directories_absent=True, current_and_original_hashes_verified=True,
                    directory_copy_used=False, zip_used=False)

    def verify_layout(self):
        return dict(conditional_selection_closure=self.verify_closure(),
                    recovered_conditional_navigation=self.verify_conditional_navigation(),
                    structure_closure=structure.ClosedStructureSelectionLayout.verify_closure(self),
                    recovered_structure_navigation=self.verify_structure_navigation(),
                    matter_closure=matter.ClosedMatterGaugeLayout.verify_closure(self),
                    recovered_matter_navigation=self.verify_matter_navigation(),
                    response_closure=matter.response.ClosedGRResponseLayout.verify_closure(self),
                    recovered_response_navigation=self.verify_response_navigation(),
                    thermal_geometry_closure=matter.geometry.ClosedThermalGeometryLayout.verify_closure(self),
                    recovered_previous_navigation=self.verify_active_navigation(),
                    prior_thermal_closure=matter.thermal.ThermalLayout.verify_closure(self),
                    frozen_layers=matter.gr.GRLayout.verify_layout(self), all_checks_passed=True)


class ClosedConditionalSelectionLoader(ResearchLoader):
    def exec_module(self, module):
        if not closed_stage(self.key):
            return previous.ConditionalSelectionLoader(self.runtime, self.key).exec_module(module)
        store = self.runtime.store
        logical, actual = ROOT / self.key, store.script_actual(self.key)
        module.__file__, module.__builtins__ = str(logical), self.runtime.scoped_builtins[VIEW]
        old_path = list(sys.path)
        sys.path.insert(0, str(logical.parent))
        try:
            exec(compile(store.read(logical, view=VIEW), str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = old_path


class ClosedConditionalSelectionRuntime(previous.ConditionalSelectionRuntime):
    def __init__(self, store):
        super().__init__(store)
        import pathlib
        proxy = types.ModuleType('pathlib')
        proxy.__dict__.update(vars(pathlib))
        proxy.Path = store.path_class(VIEW)
        default_import, default_open = builtins.__import__, builtins.open

        def scoped_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == 'pathlib' and level == 0:
                return proxy
            return default_import(name, globals, locals, fromlist, level)

        def scoped_open(file, mode='r', buffering=-1, encoding=None, errors=None,
                        newline=None, closefd=True, opener=None):
            if isinstance(file, (str, bytes, os.PathLike)):
                name = os.fsdecode(file)
                key = store.key(name, VIEW)
                if key and (key.startswith(BASE + '/') or key in store.conditional_closure_entries):
                    if opener is not None or not closefd:
                        raise NotImplementedError('Custom historical openers are unsupported')
                    return proxy.Path(name).open(mode, buffering, encoding, errors, newline)
            return default_open(file, mode, buffering, encoding, errors, newline, closefd, opener)

        self.scoped_builtins[VIEW] = dict(vars(builtins), __import__=scoped_import, open=scoped_open)
        matter.adapt_reader_imports(self, ClosedConditionalSelectionLayout, ClosedConditionalSelectionRuntime,
            ((matter, 'ClosedMatterGaugeLayout', 'ClosedMatterGaugeRuntime'),
             (previous, 'ConditionalSelectionLayout', 'ConditionalSelectionRuntime'),
             (structure, 'ClosedStructureSelectionLayout', 'ClosedStructureSelectionRuntime')))

    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, ClosedConditionalSelectionLoader(self, key),
                                              origin=str(self.store.script_actual(key)))


class ExecutionLayout(structure.ExecutionLayout):
    def actual(self, path, view='legacy'):
        return self.store.conditional_closure_actual(super().actual(path, view))


HISTORICAL_FLAGS = ('--verify-conditional-navigation', '--verify-structure-closure', *previous.HISTORICAL_FLAGS)


def runner_modes():
    modes = previous.runner_modes()
    modes['scripts/run_research_structure_selection_closed.py']['--verify-closure'] = '--verify-structure-closure'
    modes['scripts/run_research_conditional_selection_active.py'] = {
        '--verify-navigation': '--verify-conditional-navigation',
        '--verify-closure': '--verify-structure-closure',
        **{flag: flag for flag in previous.HISTORICAL_FLAGS}}
    return modes


@contextmanager
def frozen_runner_routes(target=THIS, modes=None):
    """Preserve each old flag's scope and map only recognized Python dispatch."""
    old_run = subprocess.run
    modes = runner_modes() if modes is None else modes

    def mapped(command, *args, **kwargs):
        if isinstance(command, (list, tuple)) and len(command) >= 6:
            same_python = os.path.normcase(os.path.abspath(str(command[0]))) == os.path.normcase(sys.executable)
            name = relative(command[4])
            if same_python and list(command[1:4]) == ['-B', '-X', 'utf8'] and name in modes:
                flag = command[5]
                if flag in ('--script', '--verify-layout') or flag in modes[name]:
                    command = [*command[:4], str(target), modes[name].get(flag, flag), *command[6:]]
                    if kwargs.get('cwd') is not None:
                        cwd = relative(kwargs['cwd'])
                        for old, new in (*PREFIXES.items(), *structure.PREFIXES.items(), (matter.OLD, matter.NEW),
                                         (matter.response.OLD, matter.response.NEW)):
                            if inside(cwd, old):
                                kwargs['cwd'] = str(ROOT / (new + cwd[len(old):]))
                                break
        return old_run(command, *args, **kwargs)

    subprocess.run = mapped
    try:
        yield
    finally:
        subprocess.run = old_run


def historical_check(store, flag):
    if flag == '--verify-conditional-navigation':
        return store.verify_conditional_navigation()
    if flag == '--verify-structure-closure':
        return structure.ClosedStructureSelectionLayout.verify_closure(store)
    if flag == '--verify-structure-navigation':
        return store.verify_structure_navigation()
    return structure.historical_check(store, flag)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--script')
    group.add_argument('--verify-closure', action='store_true')
    group.add_argument('--verify-navigation', action='store_true')
    group.add_argument('--verify-layout', action='store_true')
    for flag in HISTORICAL_FLAGS:
        group.add_argument(flag, action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = ClosedConditionalSelectionLayout()
    runtime = ClosedConditionalSelectionRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    matter.routes.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), matter.routes.execution_routes(ExecutionLayout(store), runtime), frozen_runner_routes():
        if args.script:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)
        else:
            selected = next((flag for flag in HISTORICAL_FLAGS if getattr(args, flag[2:].replace('-', '_'))), None)
            result = (historical_check(store, selected) if selected else
                      store.verify_layout() if args.verify_layout else
                      store.verify_closure() if args.verify_closure else store.verify_conditional_navigation())
            print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
