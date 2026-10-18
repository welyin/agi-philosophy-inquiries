"""Read-only replay beneath round 1308 conditional-selection navigation.

New-stage scripts use their current files. For closed evidence the new root
navigation prefix is first reversed in memory, then the unchanged 1288--1307
closure and all inherited readers verify and recover their historical bytes.
"""
from __future__ import annotations

import argparse
import builtins
import importlib.util
import json
import os
from pathlib import Path
import sys

import run_research_structure_selection_closed as previous
from research_layout import ROOT, ResearchLoader, no_archive_or_research_writes
from research_layout import relocation_original_bytes

THIS = Path(__file__).resolve()
STAGE = 'research_cognition_physics/archive_1308_'
LAYER = ROOT / 'research_cognition_physics/_migration/conditional_selection_active_1308/layer.json'
VIEW = 'conditional_selection_active'
NAVIGATION = {'research_cognition_physics/' + name for name in
              ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')}
ENTRY_POINTS = [STAGE + '/README.md', STAGE + '/_shared/goal.md']
digest, relative, inside = previous.digest, previous.relative, previous.inside


class ConditionalSelectionLayout(previous.ClosedStructureSelectionLayout):
    def __init__(self):
        raw = LAYER.read_bytes()
        self.conditional_layer, self.conditional_sha256 = json.loads(raw), digest(raw)
        if (self.conditional_layer.get('schema') != 'conditional_selection_active_navigation_v1'
                or self.conditional_layer.get('stage_path') != STAGE
                or self.conditional_layer.get('stage_first_round') != 1308
                or self.conditional_layer.get('entry_points') != ENTRY_POINTS):
            raise ValueError('Unexpected conditional-selection navigation schema or stage')
        self.conditional_entries = {}
        for entry in self.conditional_layer['entries']:
            old, new = previous.matter.gr.validate_entry(entry)
            if old not in NAVIGATION or new != old or old in self.conditional_entries:
                raise ValueError('Only four distinct root navigation documents may change')
            if (len(entry['edits']) != 1 or entry['edits'][0]['before'] != ''
                    or entry['edits'][0]['start'] != entry['edits'][0]['end']):
                raise ValueError('Conditional navigation must be a single prefix insertion')
            self.conditional_entries[old] = entry
        if set(self.conditional_entries) != NAVIGATION:
            raise ValueError('Exactly four root navigation documents are required')
        super().__init__()
        for name, entry in self.conditional_entries.items():
            boundary = self.structure_closure_entries[name]
            if (entry['original_sha256'], entry['original_bytes']) != (
                    boundary['current_sha256'], boundary['current_bytes']):
                raise ValueError('Navigation does not join the structure-selection closure: ' + name)

    def before_conditional_selection_navigation(self, path):
        name = self.structure_closure_name(path)
        raw = self.structure_closure_actual(path).read_bytes()
        entry = self.conditional_entries.get(name)
        return relocation_original_bytes(raw, entry) if entry else raw

    def before_structure_closure(self, path):
        name = self.structure_closure_name(path)
        if previous.closed_stage(name) and name not in self.structure_closure_entries:
            raise FileNotFoundError('File absent from structure-selection closure: ' + str(name))
        raw = self.before_conditional_selection_navigation(path)
        entry = self.structure_closure_entries.get(name)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def script_key(self, path):
        name = relative(path)
        if inside(name, STAGE):
            return name if name.endswith('.py') and (ROOT / name).is_file() else None
        return super().script_key(path)

    def script_actual(self, key):
        return ROOT / key if inside(key, STAGE) else super().script_actual(key)

    def verify_conditional_navigation(self):
        if digest(LAYER.read_bytes()) != self.conditional_sha256:
            raise ValueError('Conditional-selection layer changed during replay')
        for name, entry in self.conditional_entries.items():
            old = self.before_conditional_selection_navigation(ROOT / name)
            enc = 'utf-8-sig' if old.startswith(b'\xef\xbb\xbf') else 'utf8'
            if entry['edits'][0]['start'] != old.decode(enc).index('\n') + 1:
                raise ValueError('Navigation insertion must follow the original title: ' + name)
        support = self.conditional_layer['unchanged_support_sha256']
        required = dict(self.structure_closure_manifest['unchanged_support_sha256'])
        required[previous.MANIFEST.relative_to(ROOT).as_posix()] = self.structure_closure_sha256
        required[previous.THIS.relative_to(ROOT).as_posix()] = digest(previous.THIS.read_bytes())
        if not isinstance(support, dict) or not set(required).issubset(support):
            raise ValueError('Missing inherited support or structure-closure bindings')
        for name, expected in support.items():
            if digest((ROOT / previous.safe_relative(name)).read_bytes()) != expected:
                raise ValueError('Frozen conditional-navigation support changed: ' + name)
            if name in required and expected != required[name]:
                raise ValueError('Inherited support binding changed: ' + name)
        for name in ENTRY_POINTS:
            if not (ROOT / name).is_file():
                raise FileNotFoundError('New-stage navigation target is missing: ' + name)
        return dict(status='PASS', files=4, manifest_sha256=self.conditional_sha256,
                    closed_boundary_verified=True, original_navigation_recoverable=True,
                    frozen_archives_modified=False, stage_first_round=1308,
                    entry_points=ENTRY_POINTS, frozen_support_files=len(support),
                    scientific_acceptance_checked=False)

    def verify_layout(self):
        return dict(conditional_selection_navigation=self.verify_conditional_navigation(),
                    frozen_layers=super().verify_layout(), all_checks_passed=True)


class ConditionalSelectionLoader(ResearchLoader):
    def exec_module(self, module):
        if not inside(self.key, STAGE):
            return previous.ClosedStructureSelectionLoader(self.runtime, self.key).exec_module(module)
        actual = self.runtime.store.script_actual(self.key)
        module.__file__, module.__builtins__ = str(actual), self.runtime.scoped_builtins[VIEW]
        old_path = list(sys.path)
        sys.path.insert(0, str(actual.parent))
        try:
            exec(compile(actual.read_bytes(), str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = old_path


class ConditionalSelectionRuntime(previous.ClosedStructureSelectionRuntime):
    def __init__(self, store):
        super().__init__(store)
        self.scoped_builtins[VIEW] = dict(vars(builtins))
        previous.matter.adapt_reader_imports(
            self, ConditionalSelectionLayout, ConditionalSelectionRuntime,
            ((previous, 'ClosedStructureSelectionLayout', 'ClosedStructureSelectionRuntime'),
             (previous.previous, 'StructureSelectionLayout', 'StructureSelectionRuntime'),
             (previous.matter, 'ClosedMatterGaugeLayout', 'ClosedMatterGaugeRuntime')))

    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, ConditionalSelectionLoader(self, key),
                                              origin=str(self.store.script_actual(key)))


HISTORICAL_FLAGS = ('--verify-structure-navigation', *previous.HISTORICAL_FLAGS)


def runner_modes():
    modes = previous.runner_modes()
    modes['scripts/run_research_structure_selection_active.py']['--verify-navigation'] = '--verify-structure-navigation'
    modes['scripts/run_research_structure_selection_closed.py'] = {
        '--verify-navigation': '--verify-structure-navigation', '--verify-closure': '--verify-closure',
        **{flag: flag for flag in previous.HISTORICAL_FLAGS}}
    return modes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--script')
    group.add_argument('--verify-navigation', action='store_true')
    group.add_argument('--verify-closure', action='store_true')
    group.add_argument('--verify-layout', action='store_true')
    for flag in HISTORICAL_FLAGS:
        group.add_argument(flag, action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = ConditionalSelectionLayout()
    runtime = ConditionalSelectionRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    previous.matter.routes.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), previous.matter.routes.execution_routes(previous.ExecutionLayout(store), runtime), previous.frozen_runner_routes(THIS, runner_modes()):
        if args.script:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)
        else:
            selected = next((flag for flag in previous.HISTORICAL_FLAGS if getattr(args, flag[2:].replace('-', '_'))), None)
            result = (previous.historical_check(store, selected) if selected else
                      store.verify_structure_navigation() if args.verify_structure_navigation else
                      store.verify_layout() if args.verify_layout else
                      store.verify_closure() if args.verify_closure else store.verify_conditional_navigation())
            print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
