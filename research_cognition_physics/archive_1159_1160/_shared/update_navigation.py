"""Maintain one reversible 1159-onward navigation notice.

Edit navigation_notice.md beside this script, then run this file with --check
(or no option) to preview, and --apply to publish. Links in the notice are
relative to research_cognition_physics, where the four navigation files live.
Repeated application replaces the one current notice; it never stacks copies.
The frozen combined manifest and every older layer remain unchanged. Do not
run older navigation maintainers or edit the root navigation files directly.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research_cognition_physics'
HERE = Path(__file__).resolve().parent
NOTICE = HERE / 'navigation_notice.md'
LAYER = BASE / '_migration/local_dynamics_active_1159/layer.json'
FROZEN = BASE / '_migration/merge_1120_1158/manifest.json'
NAMES = ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')
sys.path.insert(0, str(ROOT / 'scripts'))
from research_layout import relocation_original_bytes


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group()
    choice.add_argument('--apply', action='store_true')
    choice.add_argument('--check', action='store_true')
    args = parser.parse_args()
    notice_raw = NOTICE.read_bytes()
    notice = notice_raw.decode('utf-8-sig').replace('\r\n', '\n').strip()
    if not notice.startswith('## ') or notice.count('$$') % 2:
        raise ValueError('Notice must begin with a level-two heading and have balanced display math')
    frozen_bytes = FROZEN.read_bytes()
    frozen = json.loads(frozen_bytes)
    if frozen.get('schema') != 'research_merge_1120_1158_v1':
        raise ValueError('Unexpected frozen combined manifest')
    prior = {e['destination']: e for e in frozen['entries']}
    previous_bytes = LAYER.read_bytes() if LAYER.exists() else None
    previous = json.loads(previous_bytes) if previous_bytes is not None else None
    if previous and previous.get('schema') != 'local_dynamics_active_navigation_v1':
        raise ValueError('Unexpected current navigation layer')
    existing = {e['original']: e for e in previous['entries']} if previous else {}
    names = {BASE.joinpath(name).relative_to(ROOT).as_posix() for name in NAMES}
    if existing and (set(existing) != names or len(previous['entries']) != 4):
        raise ValueError('The existing layer must cover exactly four navigation files')
    entries, changes = [], []
    for name in NAMES:
        path = BASE / name
        relative = path.relative_to(ROOT).as_posix()
        current = path.read_bytes()
        original = relocation_original_bytes(current, existing[relative]) if existing else current
        if sha(original) != prior[relative]['current_sha256'] or len(original) != prior[relative]['current_bytes']:
            raise ValueError('Navigation does not match the frozen combined boundary: ' + name)
        encoding = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf-8'
        before = original.decode(encoding)
        if not before.startswith('# ') or '\n' not in before:
            raise ValueError('Navigation must start with its existing document heading: ' + name)
        newline = '\r\n' if '\r\n' in before else '\n'
        split = before.index('\n') + 1
        inserted = newline + notice.replace('\n', newline) + newline + newline
        after = before[:split] + inserted + before[split:]
        new = after.encode(encoding)
        entry = dict(original=relative, destination=relative,
                     original_sha256=sha(original), current_sha256=sha(new),
                     original_bytes=len(original), current_bytes=len(new),
                     edits=[dict(start=split, end=split, before='', after=inserted,
                                 kind='active_stage_navigation_notice')])
        if relocation_original_bytes(new, entry) != original:
            raise ValueError('Navigation inverse failed: ' + name)
        entries.append(entry)
        changes.append((path, current, new))
    support_names = ['scripts/run_research_combined.py', FROZEN.relative_to(ROOT).as_posix(),
                     'research_cognition_physics/_migration/field_active_1120/layer.json']
    support = {name: sha((ROOT / name).read_bytes()) for name in support_names}
    if previous and previous.get('unchanged_support_sha256') != support:
        raise ValueError('Frozen support changed; do not refresh its hashes silently')
    layer = dict(schema='local_dynamics_active_navigation_v1', date=datetime.now().date().isoformat(),
                 notice_path=NOTICE.relative_to(ROOT).as_posix(), notice_sha256=sha(notice_raw),
                 entries=entries, unchanged_support_sha256=support)
    if args.apply:
        if FROZEN.read_bytes() != frozen_bytes or NOTICE.read_bytes() != notice_raw:
            raise ValueError('Inputs changed while preparing navigation')
        if (LAYER.read_bytes() if LAYER.exists() else None) != previous_bytes:
            raise ValueError('Navigation layer changed concurrently')
        for path, current, new in changes:
            if path.read_bytes() != current:
                raise ValueError('Navigation changed concurrently: ' + str(path))
        LAYER.parent.mkdir(parents=True, exist_ok=True)
        written = []
        try:
            for path, current, new in changes:
                if current != new:
                    path.write_bytes(new)
                    written.append((path, current, new))
            LAYER.write_text(json.dumps(layer, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        except Exception:
            for path, current, new in reversed(written):
                if path.read_bytes() == new:
                    path.write_bytes(current)
            if previous_bytes is not None:
                LAYER.write_bytes(previous_bytes)
            raise
    print(json.dumps(dict(status='APPLIED' if args.apply else 'CHECKED', files=4,
                          changed_files=sum(current != new for _, current, new in changes),
                          frozen_combined_unchanged=True, original_navigation_recoverable=True,
                          scientific_files_modified=False), ensure_ascii=False))


if __name__ == '__main__':
    main()
