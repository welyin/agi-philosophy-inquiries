"""Read-only historical replay after opening the field-construction stage.

Only four live root navigation files belong to this reversible outer layer.
Their pre-update bytes join the frozen 1106--1119 merge exactly. Scientific
assets, earlier manifests, filenames and per-module historical views remain
under the existing runtime; no copied directory or ZIP is used.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import sys

import run_research_merged as previous
from research_layout import ROOT, no_archive_or_research_writes, relocation_original_bytes

THIS = Path(__file__).resolve()
LAYER = ROOT / 'research_cognition_physics/_migration/field_active_1120/layer.json'
NAVIGATION = {'research_cognition_physics/' + name for name in
              ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')}


class FieldLayout(previous.MergedLayout):
    def __init__(self):
        raw = LAYER.read_bytes()
        self.field_layer_sha256 = previous.digest(raw)
        self.field_layer = json.loads(raw.decode('utf-8-sig'))
        if self.field_layer.get('schema') != 'field_active_navigation_v1':
            raise ValueError('Unsupported field navigation layer')
        entries = self.field_layer.get('entries')
        if not isinstance(entries, list):
            raise ValueError('Field navigation entries must be a list')
        self.field_entries = {}
        for entry in entries:
            name = entry.get('original')
            if (name not in NAVIGATION or entry.get('destination') != name
                    or name in self.field_entries):
                raise ValueError('Only the four in-place root navigation files may change')
            for field in ('original_sha256', 'current_sha256'):
                if not isinstance(entry.get(field), str) or not re.fullmatch(r'[0-9a-f]{64}', entry[field]):
                    raise ValueError('Invalid field navigation hash: ' + name)
            for field in ('original_bytes', 'current_bytes'):
                if type(entry.get(field)) is not int or entry[field] < 0:
                    raise ValueError('Invalid field navigation size: ' + name)
            if not isinstance(entry.get('edits'), list):
                raise ValueError('Field navigation edits must be a list: ' + name)
            for edit in entry['edits']:
                if (type(edit.get('start')) is not int or type(edit.get('end')) is not int
                        or edit['start'] < 0 or edit['end'] < edit['start']
                        or not isinstance(edit.get('before'), str)
                        or not isinstance(edit.get('after'), str)):
                    raise ValueError('Invalid field navigation edit: ' + name)
            self.field_entries[name] = entry
        if set(self.field_entries) != NAVIGATION:
            raise ValueError('The field layer must cover exactly four navigation files')
        super().__init__()
        for name, entry in self.field_entries.items():
            prior = self.merge_entries[name]
            if (entry['original_sha256'] != prior['current_sha256']
                    or entry['original_bytes'] != prior['current_bytes']):
                raise ValueError('Field navigation does not join the frozen merge: ' + name)

    def before_field(self, path):
        actual = Path(str(path))
        raw = actual.read_bytes()
        entry = self.field_entries.get(previous.relative(actual))
        return relocation_original_bytes(raw, entry) if entry else raw

    def before_merge(self, path):
        # Replace only the parent's physical-read boundary. Its recorded
        # inverse still verifies the exact bytes before every older layer.
        old = self.original_name(path)
        actual = self.merge_actual(path)
        raw = self.before_field(actual)
        entry = self.merge_entries.get(old)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def verify_field_layer(self):
        if previous.digest(LAYER.read_bytes()) != self.field_layer_sha256:
            raise ValueError('Field navigation manifest changed during verification')
        for name in self.field_entries:
            self.before_field(ROOT / name)
        support = self.field_layer.get('unchanged_support_sha256', {})
        if not isinstance(support, dict):
            raise ValueError('Field frozen-support hashes must be a mapping')
        for name, expected in support.items():
            if previous.digest((ROOT / previous.safe_relative(name)).read_bytes()) != expected:
                raise ValueError('Previous frozen support changed: ' + name)
        return dict(files=len(self.field_entries), manifest_sha256=self.field_layer_sha256,
                    frozen_support_files=len(support), current_and_original_hashes_verified=True,
                    frozen_merge_boundary_verified=True, scientific_python_json_modified=False,
                    old_archives_modified=False, directory_copy_used=False, zip_used=False)

    def verify_layout(self):
        current = self.verify_field_layer()
        older = super().verify_layout()
        return dict(field_navigation=current, frozen_layers=older, all_checks_passed=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--script')
    group.add_argument('--verify-layout', action='store_true')
    group.add_argument('--verify-field', action='store_true')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = FieldLayout()
    runtime = previous.MergedRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    routes = previous.previous.previous
    # Child checks must retain this outer inverse and physical cwd mapping.
    routes.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), routes.execution_routes(previous.ExecutionLayout(store), runtime):
        if args.verify_layout or args.verify_field:
            print(json.dumps(store.verify_layout() if args.verify_layout else store.verify_field_layer(),
                             ensure_ascii=False))
        else:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)


if __name__ == '__main__':
    main()
