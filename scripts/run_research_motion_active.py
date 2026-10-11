"""Read-only historical replay after opening the actual-SR connection stage.

Four live root navigation files and the current SR paper share this reversible
outer layer. Frozen
SR closure, active, formula and older archive verifiers retain their exact code
and manifests. This wrapper does not implement new scientific permissions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

import run_research_sr_archive as previous
from research_layout import ROOT, relocation_original_bytes, no_archive_or_research_writes

THIS = Path(__file__).resolve()
LAYER = ROOT / 'research_cognition_physics/_migration/motion_active_1106/layer.json'
NAVIGATION = {'research_cognition_physics/' + name for name in
              ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')}
PAPER = 'research_cognition_physics/认知时空操作公理与洛伦兹运动学_阶段论文.md'
LIVE_MARKDOWN = NAVIGATION | {PAPER}


class BeforeMotionFile:
    def __init__(self, store, path):
        self.store, self.path = store, path

    def read_bytes(self):
        return self.store.before_motion(self.path)

    @property
    def suffix(self):
        return self.path.suffix


class BeforeMotionClosure:
    """Read-only adapter for the frozen verifier's listed changed Markdown readers."""
    def __init__(self, store):
        self.store = store
        self.sr_closure = {
            key: (BeforeMotionFile(store, path) if key in store.motion_entries else path, entry)
            for key, (path, entry) in store.sr_closure.items()
        }

    def __getattr__(self, name):
        return getattr(self.store, name)


class MotionLayout(previous.SRLayout):
    def __init__(self):
        raw = LAYER.read_bytes()
        self.motion_layer_sha256 = hashlib.sha256(raw).hexdigest()
        self.motion_layer = json.loads(raw.decode('utf-8-sig'))
        if self.motion_layer.get('schema') != 'motion_active_navigation_v1':
            raise ValueError('Unsupported motion navigation layer')
        self.motion_entries = {}
        for entry in self.motion_layer['entries']:
            name = entry['original']
            if name not in LIVE_MARKDOWN or name != entry['destination'] or name in self.motion_entries:
                raise ValueError('Only the four root navigation files and current SR paper may change')
            self.motion_entries[name] = entry
        if not NAVIGATION <= set(self.motion_entries) <= LIVE_MARKDOWN:
            raise ValueError('Incomplete motion navigation layer')
        super().__init__()

    def before_motion(self, path):
        actual = Path(str(path))
        raw = actual.read_bytes()
        try:
            name = actual.relative_to(ROOT).as_posix()
        except ValueError:
            return raw
        entry = self.motion_entries.get(name)
        return relocation_original_bytes(raw, entry) if entry else raw

    def before_closure(self, path):
        actual = self.actual(path, 'sr')
        raw = self.before_motion(actual)
        logical = self.sr_closure_reverse.get(previous.path_key(actual))
        if logical is None:
            return raw
        prior = relocation_original_bytes(raw, self.sr_closure[logical][1])
        if hasattr(self, 'read_paths'):
            self.read_paths.add(actual.relative_to(ROOT).as_posix())
            if raw != prior:
                self.restored_markdown.add(actual.relative_to(ROOT).as_posix())
        return prior

    def verify_closure_layer(self):
        return previous.SRLayout.verify_closure_layer(BeforeMotionClosure(self))

    def verify_motion_layer(self):
        if hashlib.sha256(LAYER.read_bytes()).hexdigest() != self.motion_layer_sha256:
            raise ValueError('Motion navigation manifest changed during verification')
        for name in self.motion_entries:
            self.before_motion(ROOT / name)
        for name, digest in self.motion_layer['unchanged_support_sha256'].items():
            if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
                raise ValueError('Previous frozen support changed: ' + name)
        return dict(files=len(self.motion_entries), manifest_sha256=self.motion_layer_sha256,
                    current_and_original_hashes_verified=True,
                    old_archive_modified=False, current_paper_revised=PAPER in self.motion_entries,
                    published_paper_recoverable=True, new_scientific_rounds=0,
                    directory_copy_used=False, zip_used=False)

    def verify_layout(self):
        current = self.verify_motion_layer()
        older = super().verify_layout()
        return dict(motion_navigation=current, frozen_layers=older, all_checks_passed=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--script')
    group.add_argument('--verify-layout', action='store_true')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = MotionLayout()
    runtime = previous.SRRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    # Child replays must keep this newest live-Markdown inverse layer.
    previous.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), previous.execution_routes(store, runtime):
        if args.verify_layout:
            print(json.dumps(store.verify_layout(), ensure_ascii=False))
        else:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)


if __name__ == '__main__':
    main()
