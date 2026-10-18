"""Read-only historical replay beneath the round-1313 navigation layer."""
from __future__ import annotations

import argparse
import builtins
import importlib.util
import json
import os
from pathlib import Path
import sys

import run_research_conditional_selection_closed as previous
from research_layout import ROOT, ResearchLoader, no_archive_or_research_writes
from research_layout import relocation_original_bytes

THIS = Path(__file__).resolve()
STAGE = 'research_cognition_physics/archive_1313_'
LAYER = ROOT / 'research_cognition_physics/_migration/interaction_generation_active_1313/layer.json'
VIEW = 'interaction_generation_active'
NAVIGATION = {'research_cognition_physics/' + name for name in
              ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')}
ENTRY_POINTS = [STAGE + '/README.md', STAGE + '/_shared/goal.md']
digest, relative, inside = previous.digest, previous.relative, previous.inside


class InteractionGenerationLayout(previous.ClosedConditionalSelectionLayout):
    def __init__(self):
        raw = LAYER.read_bytes()
        self.interaction_layer, self.interaction_sha256 = json.loads(raw), digest(raw)
        if (self.interaction_layer.get('schema') != 'interaction_generation_active_navigation_v1'
                or self.interaction_layer.get('stage_path') != STAGE
                or self.interaction_layer.get('stage_first_round') != 1313
                or self.interaction_layer.get('entry_points') != ENTRY_POINTS):
            raise ValueError('Unexpected interaction-generation navigation schema or stage')
        self.interaction_entries = {}
        for entry in self.interaction_layer['entries']:
            old, new = previous.matter.gr.validate_entry(entry)
            if old not in NAVIGATION or new != old or old in self.interaction_entries:
                raise ValueError('Only four distinct root navigation documents may change')
            if (len(entry['edits']) != 1 or entry['edits'][0]['before'] != ''
                    or entry['edits'][0]['start'] != entry['edits'][0]['end']):
                raise ValueError('Active navigation must be a single prefix insertion')
            self.interaction_entries[old] = entry
        if set(self.interaction_entries) != NAVIGATION:
            raise ValueError('Exactly four root navigation documents are required')
        super().__init__()
        for name, entry in self.interaction_entries.items():
            boundary = self.conditional_closure_entries[name]
            if (entry['original_sha256'], entry['original_bytes']) != (
                    boundary['current_sha256'], boundary['current_bytes']):
                raise ValueError('Navigation does not join the 1308--1312 closure: ' + name)

    def before_interaction_navigation(self, path):
        name = self.conditional_closure_name(path)
        raw = self.conditional_closure_actual(path).read_bytes()
        entry = self.interaction_entries.get(name)
        return relocation_original_bytes(raw, entry) if entry else raw

    def before_conditional_closure(self, path):
        name = self.conditional_closure_name(path)
        if previous.closed_stage(name) and name not in self.conditional_closure_entries:
            raise FileNotFoundError('File absent from conditional-selection closure: ' + str(name))
        raw = self.before_interaction_navigation(path)
        entry = self.conditional_closure_entries.get(name)
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

    def verify_interaction_navigation(self):
        if digest(LAYER.read_bytes()) != self.interaction_sha256:
            raise ValueError('Interaction-generation layer changed during replay')
        for name, entry in self.interaction_entries.items():
            old = self.before_interaction_navigation(ROOT / name)
            enc = 'utf-8-sig' if old.startswith(b'\xef\xbb\xbf') else 'utf8'
            if entry['edits'][0]['start'] != old.decode(enc).index('\n') + 1:
                raise ValueError('Navigation insertion must follow the original title: ' + name)
        support = self.interaction_layer['unchanged_support_sha256']
        required = dict(self.conditional_closure_manifest['unchanged_support_sha256'])
        required[previous.MANIFEST.relative_to(ROOT).as_posix()] = self.conditional_closure_sha256
        required[previous.THIS.relative_to(ROOT).as_posix()] = digest(previous.THIS.read_bytes())
        if not isinstance(support, dict) or not set(required).issubset(support):
            raise ValueError('Missing inherited support or closure bindings')
        for name, expected in support.items():
            if digest((ROOT / previous.safe_relative(name)).read_bytes()) != expected:
                raise ValueError('Frozen navigation support changed: ' + name)
            if name in required and expected != required[name]:
                raise ValueError('Inherited support binding changed: ' + name)
        for name in ENTRY_POINTS:
            if not (ROOT / name).is_file():
                raise FileNotFoundError('New-stage navigation target is missing: ' + name)
        return dict(status='PASS', files=4, manifest_sha256=self.interaction_sha256,
                    closed_boundary_verified=True, original_navigation_recoverable=True,
                    frozen_archives_modified=False, stage_first_round=1313,
                    entry_points=ENTRY_POINTS, frozen_support_files=len(support),
                    scientific_acceptance_checked=False)

    def verify_layout(self):
        return dict(interaction_generation_navigation=self.verify_interaction_navigation(),
                    frozen_layers=super().verify_layout(), all_checks_passed=True)


class InteractionGenerationLoader(ResearchLoader):
    def exec_module(self, module):
        if not inside(self.key, STAGE):
            return previous.ClosedConditionalSelectionLoader(self.runtime, self.key).exec_module(module)
        actual = self.runtime.store.script_actual(self.key)
        module.__file__, module.__builtins__ = str(actual), self.runtime.scoped_builtins[VIEW]
        old_path = list(sys.path)
        sys.path.insert(0, str(actual.parent))
        try:
            exec(compile(actual.read_bytes(), str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = old_path


class InteractionGenerationRuntime(previous.ClosedConditionalSelectionRuntime):
    def __init__(self, store):
        super().__init__(store)
        self.scoped_builtins[VIEW] = dict(vars(builtins))
        previous.matter.adapt_reader_imports(
            self, InteractionGenerationLayout, InteractionGenerationRuntime,
            ((previous, 'ClosedConditionalSelectionLayout', 'ClosedConditionalSelectionRuntime'),
             (previous.previous, 'ConditionalSelectionLayout', 'ConditionalSelectionRuntime'),
             (previous.structure, 'ClosedStructureSelectionLayout', 'ClosedStructureSelectionRuntime'),
             (previous.matter, 'ClosedMatterGaugeLayout', 'ClosedMatterGaugeRuntime')))

    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, InteractionGenerationLoader(self, key),
                                              origin=str(self.store.script_actual(key)))


HISTORICAL_FLAGS = previous.HISTORICAL_FLAGS


def runner_modes():
    modes = previous.runner_modes()
    modes['scripts/run_research_conditional_selection_closed.py'] = {
        '--verify-navigation': '--verify-conditional-navigation',
        '--verify-closure': '--verify-closure',
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
    store = InteractionGenerationLayout()
    runtime = InteractionGenerationRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    previous.matter.routes.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), previous.matter.routes.execution_routes(previous.ExecutionLayout(store), runtime), previous.frozen_runner_routes(THIS, runner_modes()):
        if args.script:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)
        else:
            selected = next((flag for flag in HISTORICAL_FLAGS if getattr(args, flag[2:].replace('-', '_'))), None)
            result = (previous.historical_check(store, selected) if selected else
                      store.verify_layout() if args.verify_layout else
                      store.verify_closure() if args.verify_closure else store.verify_interaction_navigation())
            print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
