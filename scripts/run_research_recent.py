"""Read-only replay after the 1044--1062 archive directory migration.

The newest manifest records physical moves and reversible Markdown edits.
Every read verifies real current bytes before reconstructing the earlier view;
the older runtime files, scientific files and frozen receipts remain untouched.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import sys

import run_research_closed as closed_runtime
from research_layout import (Layout, ResearchRuntime, ROOT,
                             relocation_original_bytes,
                             no_archive_or_research_writes)

THIS = Path(__file__).resolve()
MANIFEST = ROOT / 'research_cognition_physics/_migration/layout_1044_1062_20261009/manifest.json'
PREFIXES = {
    'research_cognition_physics/archive_1044_':
        'research_cognition_physics/archive_1044_1045',
    'research_cognition_physics/archive_1046_':
        'research_cognition_physics/archive_1046_1062',
}


class RecentLayout(closed_runtime.ClosedLayout):
    """Add a separate newest layer without replacing any older manifest entry."""

    def __init__(self):
        # Parent construction reads older manifests through read_pre_split.
        self.recent = {}
        self.recent_reverse = {}
        super().__init__()
        data = json.loads(MANIFEST.read_text(encoding='utf-8'))
        if data['schema'] != 'research_round_archive_layout_v1':
            raise ValueError('Unsupported recent archive manifest')
        self.recent_manifest = data
        for entry in data['entries']:
            logical = self._safe_relative(entry['original'])
            destination = self._safe_relative(entry['destination'])
            actual = ROOT / destination
            if not actual.resolve().is_relative_to(ROOT):
                raise ValueError('Destination escapes project: ' + destination)
            moved = self.recent_key(ROOT / logical)
            if moved:
                if destination != self.actual(ROOT / logical, 'current').relative_to(ROOT).as_posix():
                    raise ValueError('Unexpected stage destination: ' + destination)
            elif logical != destination or not destination.endswith('.md'):
                raise ValueError('External entry must be an in-place Markdown edit: ' + logical)
            physical = os.path.normcase(os.path.abspath(actual))
            if logical in self.recent or physical in self.recent_reverse:
                raise ValueError('Duplicate recent migration entry: ' + logical)
            self.recent[logical] = (actual, entry)
            self.recent_reverse[physical] = logical
            self.current_files.add(logical)
            self.current_directories.update(p.as_posix() for p in PurePosixPath(logical).parents)

    @staticmethod
    def _safe_relative(value):
        path = PurePosixPath(value.replace('\\', '/'))
        if path.is_absolute() or '..' in path.parts or not path.parts or ':' in path.parts[0]:
            raise ValueError('Unsafe manifest path: ' + value)
        if path.parts[0] not in ('research_cognition_physics', '猜想'):
            raise ValueError('Recent archive entry outside permitted trees: ' + value)
        if path.parts[0] == '猜想' and path.suffix != '.md':
            raise ValueError('Only in-place methodology Markdown edits are permitted: ' + value)
        return path.as_posix()

    @staticmethod
    def recent_key(path):
        try:
            relative = Path(os.path.abspath(str(path))).relative_to(ROOT).as_posix()
        except ValueError:
            return None
        for old, new in PREFIXES.items():
            for prefix in (old, new):
                if relative == prefix or relative.startswith(prefix + '/'):
                    return old + relative[len(prefix):]
        return None

    def key(self, path, view='legacy'):
        if view == 'current':
            closed_current = self.current_closed_name(path)
            if closed_current:
                return closed_current
        return self.recent_key(path) or super().key(path, view)

    @staticmethod
    def current_closed_name(path):
        """A post-1043 caller can explicitly retain the already-closed name."""
        try:
            relative = Path(os.path.abspath(str(path))).relative_to(ROOT).as_posix()
        except ValueError:
            return None
        prefix = closed_runtime.NEW
        return relative if relative == prefix or relative.startswith(prefix + '/') else None

    def actual(self, path, view='legacy'):
        key = self.recent_key(path)
        if key:
            for old, new in PREFIXES.items():
                if key == old or key.startswith(old + '/'):
                    return ROOT / (new + key[len(old):])
        return super().actual(path, view)

    def before_recent(self, path):
        """Read a physical file and prove its bytes before the newest migration."""
        actual = self.actual(path, 'current')
        physical = os.path.normcase(os.path.abspath(actual))
        logical = self.recent_reverse.get(physical)
        raw = actual.read_bytes()
        if logical is None:
            return raw
        entry = self.recent[logical][1]
        original = relocation_original_bytes(raw, entry)
        name = actual.relative_to(ROOT).as_posix()
        self.read_paths.add(name)
        if original != raw:
            self.restored_markdown.add(name)
        return original

    def read_pre_split(self, path):
        # Reconstruct the newest layer first, then let the frozen earlier
        # manifest prove its own pre-migration bytes. In-place edits to the
        # 1009--1043 README therefore cannot conceal an older hash mismatch.
        key = self.key(path, 'current')
        if key in self.relocations:
            actual, entry = self.relocations[key]
            raw = self.before_recent(actual)
            original = relocation_original_bytes(raw, entry)
            name = actual.relative_to(ROOT).as_posix()
            self.read_paths.add(name)
            if raw != original:
                self.restored_markdown.add(name)
            return original
        return self.before_recent(path)

    def read(self, path, historical=True, view='legacy'):
        if view == 'current' and (not self.closed_key(path) or self.current_closed_name(path)):
            return self.before_recent(path) if historical else self.actual(path, view).read_bytes()
        # Skip ClosedLayout's direct physical read of unrelated current files:
        # those files can now include root navigation covered by this layer.
        return Layout.read(self, path, historical=historical, view=view)

    def script_key(self, path):
        key = self.recent_key(path)
        if key in self.recent and key.endswith('.py'):
            return key
        # Executing an older module still gives that module its own original
        # __file__ and path view, even when its caller used the closed name.
        closed = self.closed_key(path)
        if closed in self.relocations and closed.endswith('.py'):
            return closed
        return super().script_key(path)

    def script_actual(self, key):
        if key in self.recent and key.endswith('.py'):
            return self.recent[key][0]
        return super().script_actual(key)

    def lookup_script(self, name):
        try:
            return super().lookup_script(name)
        except ValueError:
            matches = [key for key in self.recent
                       if key.endswith('/' + name.replace('\\', '/')) and key.endswith('.py')]
            if len(matches) == 1:
                return matches[0]
            raise

    def glob_paths(self, directory, pattern, view):
        if view == 'current' and self.current_closed_name(directory):
            physical = self.actual(directory, view)
            if physical.is_dir():
                yield from sorted(physical.glob(pattern))
            return
        if self.recent_key(directory):
            physical = self.actual(directory, view)
            if physical.is_dir():
                yield from (ROOT / self.recent_key(p) for p in sorted(physical.glob(pattern)))
            return
        yield from super().glob_paths(directory, pattern, view)

    def path_class(self, view='legacy'):
        inherited = super().path_class(view)
        if view != 'current':
            return inherited
        native = type(Path())

        class RecentPath(inherited):
            def __init__(self, *parts):
                raw = native(*parts)
                if raw.is_absolute():
                    try:
                        relative = raw.relative_to(ROOT).as_posix()
                    except ValueError:
                        relative = ''
                    for old, new in PREFIXES.items():
                        for prefix in (old, new):
                            if relative == prefix or relative.startswith(prefix + '/'):
                                # pathlib retains lexical '..' until resolve().
                                # Frozen receipts sometimes include that exact
                                # spelling as a dictionary key. Normalize only
                                # during actual file lookup, not construction.
                                native.__init__(self, ROOT / (old + relative[len(prefix):]))
                                return
                super().__init__(*parts)

        return RecentPath

    def verify_affected_prior_layers(self):
        # Check only overlaps with this migration. These three live root
        # navigation documents had already evolved since the 764--1008
        # migration; the newest manifest proves their actual pre-move bytes,
        # but it does not recreate that unrelated older navigation snapshot.
        evolving_navigation = {
            'research_cognition_physics/README.md',
            'research_cognition_physics/research_direction.md',
            'research_cognition_physics/RESEARCH_STATE.md',
        }
        verified, prior_navigation_differences = [], []
        for actual, entry in self.relocations.values():
            physical = os.path.normcase(os.path.abspath(actual))
            if physical not in self.recent_reverse:
                continue
            raw = self.before_recent(actual)
            digest = hashlib.sha256(raw).hexdigest()
            if digest != entry['current_sha256'] or len(raw) != entry['current_bytes']:
                recent_entry = self.recent[self.recent_reverse[physical]][1]
                if (entry['destination'] not in evolving_navigation
                        or digest != recent_entry['original_sha256']):
                    raise ValueError('Affected prior entry changed: ' + entry['destination'])
                prior_navigation_differences.append(dict(
                    path=entry['destination'], older_layout_sha256=entry['current_sha256'],
                    before_recent_migration_sha256=digest,
                    newest_layer_original_bytes_verified=True))
                continue
            original = relocation_original_bytes(raw, entry)
            if actual.suffix in ('.py', '.json') and raw != original:
                raise ValueError('Prior scientific bytes changed: ' + str(actual))
            verified.append(entry['destination'])
        return dict(scope='Only prior manifest entries affected by this migration',
                    verified_prior_entries=verified,
                    preexisting_live_navigation_differences=prior_navigation_differences,
                    whole_older_layout_reverified=False)

    def verify_closure(self):
        if (ROOT / closed_runtime.OLD).exists():
            raise ValueError('Old 1009 directory must not be retained')
        changed_markdown = 0
        for entry in self.closure['entries']:
            actual = ROOT / entry['destination']
            raw = self.before_recent(actual)
            original = relocation_original_bytes(raw, entry)
            if actual.suffix in ('.py', '.json') and raw != original:
                raise ValueError('Scientific bytes changed: ' + str(actual))
            changed_markdown += raw != original
        return dict(files=len(self.closure['entries']),
                    reversible_markdown_files=changed_markdown,
                    scientific_python_json_bytes_unchanged=True,
                    frozen_runtime_files_unchanged=True,
                    zip_used=False, old_directory_retained=False)

    def verify_layout(self):
        identical = rewritten = moved = 0
        expected = {new: set() for new in PREFIXES.values()}
        for actual, entry in self.recent.values():
            raw = actual.read_bytes()
            original = relocation_original_bytes(raw, entry)
            if actual.suffix in ('.py', '.json') and raw != original:
                raise ValueError('Scientific bytes changed: ' + str(actual))
            identical += raw == original
            rewritten += raw != original
            if entry['original'] != entry['destination']:
                moved += 1
            for new in expected:
                if entry['destination'].startswith(new + '/'):
                    expected[new].add(entry['destination'])
        for old, new in PREFIXES.items():
            if (ROOT / old).exists() or not (ROOT / new).is_dir():
                raise ValueError('Archive move incomplete: ' + old)
            actual_files = {p.relative_to(ROOT).as_posix()
                            for p in (ROOT / new).rglob('*') if p.is_file()}
            if actual_files != expected[new]:
                raise ValueError('Archive inventory differs from manifest: ' + new)
        for name, digest in self.recent_manifest.get('unchanged_runtime_sha256', {}).items():
            if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
                raise ValueError('Earlier runtime changed: ' + name)
        return dict(
            recent_layout=dict(files=len(self.recent), moved_files=moved,
                               byte_identical_files=identical,
                               reversible_markdown_files=rewritten,
                               scientific_python_json_bytes_unchanged=True,
                               original_and_current_hashes_verified=True,
                               old_directories_retained=False),
            affected_prior_layers=self.verify_affected_prior_layers(),
            prior_1009_closure=self.verify_closure(),
            all_checks_passed=True,
        )


@contextmanager
def execution_routes(store, runtime):
    # Reuse the existing, narrow runpy/subprocess adapter. Its child dispatch
    # target is a process-local module variable; its file remains unchanged.
    previous = closed_runtime.THIS
    closed_runtime.THIS = THIS
    try:
        with closed_runtime.execution_routes(store, runtime):
            yield
    finally:
        closed_runtime.THIS = previous


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--script')
    choice.add_argument('--verify-layout', action='store_true')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    # Preserve the caller's numerical thread environment: frozen strict
    # floating-point comparisons can depend on the original BLAS settings.
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = RecentLayout()
    runtime = ResearchRuntime(store)
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
