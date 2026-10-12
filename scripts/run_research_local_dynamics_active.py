"""Read-only replay through the active local-dynamics navigation layer.

Only four root navigation files are reversed before the frozen 1120--1158
merge. New-stage scripts retain their actual paths and current files; imported
historical scripts retain their own existing replay views.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import sys

import run_research_combined as previous
from research_layout import ROOT, ResearchLoader, no_archive_or_research_writes, relocation_original_bytes

THIS = Path(__file__).resolve()
STAGE = 'research_cognition_physics/archive_1159_'
LAYER = ROOT / 'research_cognition_physics/_migration/local_dynamics_active_1159/layer.json'
NAVIGATION = {'research_cognition_physics/' + name for name in
              ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')}


class LocalDynamicsLayout(previous.CombinedLayout):
    def __init__(self):
        raw = LAYER.read_bytes()
        self.dynamics_sha256 = previous.digest(raw)
        self.dynamics_layer = json.loads(raw.decode('utf-8-sig'))
        if self.dynamics_layer.get('schema') != 'local_dynamics_active_navigation_v1':
            raise ValueError('Unsupported local-dynamics navigation layer')
        entries = self.dynamics_layer.get('entries')
        if not isinstance(entries, list):
            raise ValueError('Navigation entries must be a list')
        self.dynamics_entries = {}
        for entry in entries:
            name = entry.get('original')
            if name not in NAVIGATION or entry.get('destination') != name or name in self.dynamics_entries:
                raise ValueError('Only the four in-place root navigation files may change')
            for field in ('original_sha256', 'current_sha256'):
                if not isinstance(entry.get(field), str) or not re.fullmatch(r'[0-9a-f]{64}', entry[field]):
                    raise ValueError('Invalid navigation hash: ' + name)
            for field in ('original_bytes', 'current_bytes'):
                if type(entry.get(field)) is not int or entry[field] < 0:
                    raise ValueError('Invalid navigation size: ' + name)
            edits = entry.get('edits')
            if not isinstance(edits, list):
                raise ValueError('Navigation edits must be a list: ' + name)
            for edit in edits:
                if (type(edit.get('start')) is not int or type(edit.get('end')) is not int
                        or edit['start'] < 0 or edit['end'] < edit['start']
                        or not isinstance(edit.get('before'), str) or not isinstance(edit.get('after'), str)):
                    raise ValueError('Invalid navigation edit: ' + name)
            self.dynamics_entries[name] = entry
        if set(self.dynamics_entries) != NAVIGATION:
            raise ValueError('Exactly four navigation files are required')
        super().__init__()
        for name, entry in self.dynamics_entries.items():
            prior = self.combined_entries[name]
            if (entry['original_sha256'] != prior['current_sha256']
                    or entry['original_bytes'] != prior['current_bytes']):
                raise ValueError('Navigation does not join the frozen combined boundary: ' + name)

    def before_dynamics(self, path):
        actual = Path(str(path))
        raw = actual.read_bytes()
        entry = self.dynamics_entries.get(previous.relative(actual))
        return relocation_original_bytes(raw, entry) if entry else raw

    def before_combined(self, path):
        old = self.combined_original_name(path)
        actual = self.combined_actual(path)
        if old and previous.moved_stage(old) and old not in self.combined_entries:
            raise FileNotFoundError('Historical file missing from combined manifest: ' + old)
        raw = self.before_dynamics(actual)
        entry = self.combined_entries.get(old)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def script_key(self, path):
        name = previous.relative(path)
        if name and name.startswith(STAGE + '/') and name.endswith('.py') and Path(str(path)).is_file():
            return name
        return super().script_key(path)

    def script_actual(self, key):
        return ROOT / key if key.startswith(STAGE + '/') else super().script_actual(key)

    def verify_dynamics(self):
        if previous.digest(LAYER.read_bytes()) != self.dynamics_sha256:
            raise ValueError('Navigation layer changed during verification')
        for name in self.dynamics_entries:
            self.before_dynamics(ROOT / name)
        support = self.dynamics_layer.get('unchanged_support_sha256')
        if not isinstance(support, dict) or not support:
            raise ValueError('Frozen support hashes are required')
        for name, expected in support.items():
            if previous.digest((ROOT / previous.safe_relative(name)).read_bytes()) != expected:
                raise ValueError('Frozen navigation support changed: ' + name)
        return dict(files=4, manifest_sha256=self.dynamics_sha256,
                    frozen_support_files=len(support), current_and_original_hashes_verified=True,
                    combined_boundary_verified=True, frozen_archives_modified=False,
                    scientific_assets_modified=False, directory_copy_used=False, zip_used=False)

    def verify_layout(self):
        return dict(local_dynamics_navigation=self.verify_dynamics(),
                    frozen_layers=super().verify_layout(), all_checks_passed=True)


class LocalDynamicsLoader(ResearchLoader):
    def exec_module(self, module):
        if not self.key.startswith(STAGE + '/'):
            return previous.CombinedLoader(self.runtime, self.key).exec_module(module)
        actual = self.runtime.store.script_actual(self.key)
        module.__file__ = str(actual)
        old_path = list(sys.path)
        sys.path.insert(0, str(actual.parent))
        try:
            exec(compile(actual.read_bytes(), str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = old_path


class LocalDynamicsRuntime(previous.CombinedRuntime):
    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, LocalDynamicsLoader(self, key),
                                              origin=str(self.store.script_actual(key)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--script')
    choice.add_argument('--verify-layout', action='store_true')
    choice.add_argument('--verify-navigation', action='store_true')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = LocalDynamicsLayout()
    runtime = LocalDynamicsRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    routes = previous.merged.previous.previous
    routes.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), routes.execution_routes(previous.ExecutionLayout(store), runtime):
        if args.verify_layout or args.verify_navigation:
            print(json.dumps(store.verify_layout() if args.verify_layout else store.verify_dynamics(),
                             ensure_ascii=False))
        else:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)


if __name__ == '__main__':
    main()
