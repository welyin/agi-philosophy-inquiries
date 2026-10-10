"""Replay frozen research with a reversible layer for four live root documents.

Only archive_1063_/_shared/navigation_layer.json owns the new navigation edits.
The archived manifests, research bytes and earlier runtime files stay frozen.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import sys

import run_research_recent as recent_runtime
from research_layout import (ResearchRuntime, ROOT, relocation_original_bytes,
                             no_archive_or_research_writes)

THIS = Path(__file__).resolve()
LAYER = ROOT / 'research_cognition_physics/archive_1063_/_shared/navigation_layer.json'
NAVIGATION = {
    'research_cognition_physics/README.md',
    'research_cognition_physics/research_direction.md',
    'research_cognition_physics/RESEARCH_STATE.md',
    'research_cognition_physics/ROADMAP.md',
}


class LiveLayout(recent_runtime.RecentLayout):
    def __init__(self):
        self.live = {}
        super().__init__()
        self.navigation_layer = json.loads(LAYER.read_text(encoding='utf-8'))
        if self.navigation_layer['schema'] != 'research_live_navigation_layer_v1':
            raise ValueError('Unsupported live navigation layer')
        seen = set()
        for entry in self.navigation_layer['entries']:
            logical = entry['original']
            if logical not in NAVIGATION or entry['destination'] != logical:
                raise ValueError('Live layer is restricted to four in-place root documents')
            if logical in seen:
                raise ValueError('Duplicate live navigation entry: ' + logical)
            prior = self.recent[logical][1]
            if (entry['original_sha256'] != prior['current_sha256']
                    or entry['original_bytes'] != prior['current_bytes']):
                raise ValueError('Live layer does not join the frozen recent manifest: ' + logical)
            actual = ROOT / logical
            if not actual.resolve().is_relative_to(ROOT):
                raise ValueError('Live navigation escapes project: ' + logical)
            self.live[os.path.normcase(os.path.abspath(actual))] = entry
            seen.add(logical)
        if seen != NAVIGATION:
            raise ValueError('The live layer must cover exactly the four root documents')

    def before_live(self, path):
        """Verify live physical bytes, then return the exact archived navigation."""
        actual = self.actual(path, 'current')
        raw = actual.read_bytes()
        entry = self.live.get(os.path.normcase(os.path.abspath(actual)))
        if entry is None:
            return raw
        original = relocation_original_bytes(raw, entry)
        self.read_paths.add(entry['destination'])
        if original != raw:
            self.restored_markdown.add(entry['destination'])
        return original

    def before_recent(self, path):
        # The parent reader has a direct native read; repeat only that small
        # layer boundary so live bytes are reversed before the frozen layer.
        actual = self.actual(path, 'current')
        raw = self.before_live(actual)
        logical = self.recent_reverse.get(os.path.normcase(os.path.abspath(actual)))
        if logical is None:
            return raw
        original = relocation_original_bytes(raw, self.recent[logical][1])
        name = actual.relative_to(ROOT).as_posix()
        self.read_paths.add(name)
        if original != raw:
            self.restored_markdown.add(name)
        return original

    def verify_layout(self):
        # This is the previous verifier with before_live at its physical-read
        # boundary. All 604 old entries, hashes and inventory checks remain.
        identical = rewritten = moved = 0
        expected = {new: set() for new in recent_runtime.PREFIXES.values()}
        for actual, entry in self.recent.values():
            raw = self.before_live(actual)
            original = relocation_original_bytes(raw, entry)
            if actual.suffix in ('.py', '.json') and raw != original:
                raise ValueError('Scientific bytes changed: ' + str(actual))
            identical += raw == original
            rewritten += raw != original
            moved += entry['original'] != entry['destination']
            for new in expected:
                if entry['destination'].startswith(new + '/'):
                    expected[new].add(entry['destination'])
        for old, new in recent_runtime.PREFIXES.items():
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
            live_navigation_layer=dict(
                files=len(self.live),
                original_and_current_hashes_verified=True,
                frozen_recent_manifest_boundary_verified=True,
                layer_sha256=hashlib.sha256(LAYER.read_bytes()).hexdigest()),
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
    # Child Python processes need this same live layer. The target adjustment
    # is process-local and leaves both imported runtime source files intact.
    previous = recent_runtime.THIS
    recent_runtime.THIS = THIS
    try:
        with recent_runtime.execution_routes(store, runtime):
            yield
    finally:
        recent_runtime.THIS = previous


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--script')
    choice.add_argument('--verify-layout', action='store_true')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = LiveLayout()
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
