"""Read-only replay through the current active-document revision layer.

The outer layer verifies current Markdown and recovers its previous bytes before
the frozen formula, SR-navigation and archive layers run. No old manifest or
scientific asset is rewritten. The active layer contains reversible text edits,
not a second collection of document snapshots.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import sys
import types

import run_research_formatted as formatted
from research_layout import (
    ROOT,
    Layout,
    no_archive_or_research_writes,
    relocation_original_bytes,
)

THIS = Path(__file__).resolve()
ACTIVE_LAYER = ROOT / 'research_cognition_physics/_migration/active_documents/layer.json'


def _path_key(path):
    return os.path.normcase(os.path.abspath(path))


class _BeforeActiveFile:
    """The unchanged formula verifier only needs ROOT / name and read_bytes."""

    def __init__(self, store, path):
        self.store = store
        self.path = path

    def __truediv__(self, part):
        return _BeforeActiveFile(self.store, self.path / part)

    def read_bytes(self):
        return self.store.before_active(self.path)


class ActiveLayout(formatted.FormattedLayout):
    def __init__(self):
        layer_bytes = ACTIVE_LAYER.read_bytes()
        self.active_layer_sha256 = hashlib.sha256(layer_bytes).hexdigest()
        self.active_layer = json.loads(layer_bytes.decode('utf-8-sig'))
        if self.active_layer.get('schema') != 'research_active_documents_v1':
            raise ValueError('Unsupported active-document layer')
        entries = self.active_layer.get('entries')
        if not isinstance(entries, list):
            raise ValueError('Active-document entries must be a list')
        self.active = {}
        resolved_paths = set()
        for entry in entries:
            name = entry.get('destination')
            if not isinstance(name, str) or entry.get('original') != name:
                raise ValueError('Active revisions must be in-place Markdown edits')
            relative = PurePosixPath(name)
            if (not name or '\\' in name or ':' in name or relative.is_absolute()
                    or '..' in relative.parts or relative.as_posix() != name
                    or relative.suffix != '.md'
                    or set(relative.parts) & {'.git', '.codex', '.agents', '.research_runtime'}):
                raise ValueError('Invalid active-document path: ' + name)
            path = ROOT / name
            resolved = path.resolve()
            if not resolved.is_relative_to(ROOT):
                raise ValueError('Active-document path escapes the project: ' + name)
            key = _path_key(path)
            if key in self.active or resolved in resolved_paths:
                raise ValueError('Duplicate active-document path: ' + name)
            for field in ('current_sha256', 'original_sha256'):
                if not isinstance(entry.get(field), str) or not re.fullmatch(r'[0-9a-f]{64}', entry[field]):
                    raise ValueError('Invalid active-document hash: ' + name)
            for field in ('current_bytes', 'original_bytes'):
                if type(entry.get(field)) is not int or entry[field] < 0:
                    raise ValueError('Invalid active-document size: ' + name)
            if not isinstance(entry.get('edits'), list):
                raise ValueError('Active-document edits must be a list: ' + name)
            for edit in entry['edits']:
                if (type(edit.get('start')) is not int or type(edit.get('end')) is not int
                        or edit['start'] < 0 or edit['end'] < edit['start']
                        or not isinstance(edit.get('before'), str)
                        or not isinstance(edit.get('after'), str)):
                    raise ValueError('Invalid active-document edit: ' + name)
            self.active[key] = (path, entry)
            resolved_paths.add(resolved)
        super().__init__()

    def before_active(self, path):
        """Read real current bytes; verify both hashes before returning prior text."""
        actual = Path(str(path))
        raw = actual.read_bytes()
        item = self.active.get(_path_key(actual))
        return relocation_original_bytes(raw, item[1]) if item else raw

    def physical_bytes(self, path):
        actual = Path(str(path))
        raw = self.before_active(actual)
        key = _path_key(actual)
        if key in self.formats:
            raw = formatted.formatting.restore_bytes(raw, self.formats[key])
        if key in self.sr:
            raw = relocation_original_bytes(raw, self.sr[key])
        return raw

    def verify_active_layer(self):
        # Recheck the manifest too: a concurrent replacement must not be hidden
        # by the parsed in-memory version.
        if hashlib.sha256(ACTIVE_LAYER.read_bytes()).hexdigest() != self.active_layer_sha256:
            raise ValueError('Active-document manifest changed during verification')
        changed = 0
        for path, entry in self.active.values():
            original = self.before_active(path)
            changed += entry['current_sha256'] != hashlib.sha256(original).hexdigest()
        return {
            'schema': self.active_layer['schema'],
            'manifest_sha256': self.active_layer_sha256,
            'files': len(self.active),
            'reversible_markdown_files': changed,
            'current_and_original_hashes_verified': True,
            'scope': 'listed in-place active Markdown revisions only',
            'scientific_python_json_modified': False,
            'snapshot_copies_used': False,
        }

    def verify_formula_layer(self):
        # Run the original verifier's exact code with only its ROOT reader
        # replaced. It sees post-format / pre-active bytes, then validates its
        # own frozen inventory, hashes, formula bodies and exact inverse. No
        # global pathlib patch and no replacement digest is used.
        original_verify = formatted.formatting.verify
        namespace = dict(original_verify.__globals__)
        namespace['ROOT'] = _BeforeActiveFile(self, ROOT)
        verify = types.FunctionType(
            original_verify.__code__, namespace, original_verify.__name__,
            original_verify.__defaults__, original_verify.__closure__,
        )
        return verify(self.format_layer)

    def verify_layout(self):
        active = self.verify_active_layer()
        formula = self.verify_formula_layer()
        prior = formatted.spatial.SpatialLayout.verify_layout(self)
        old = Layout.verify(self)
        return {
            'active_document_layer': active,
            'format_layer': formula,
            'sr_navigation': {'files': len(self.sr), 'inverse_hashes_verified': True},
            'spatial_and_prior_layout': prior,
            'older_source_layout': old,
            'all_checks_passed': True,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--script')
    group.add_argument('--verify-layout', action='store_true')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = ActiveLayout()
    runtime = formatted.spatial.SpatialRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    # Child replays must retain this outer layer instead of dropping back to
    # the formula-only wrapper.
    formatted.spatial.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), formatted.spatial.execution_routes(store, runtime):
        if args.verify_layout:
            print(json.dumps(store.verify_layout(), ensure_ascii=False))
        else:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)


if __name__ == '__main__':
    main()
