"""Physically close rounds 1063--1085 with byte-preserving evidence mapping.

Only Markdown directory/runtime references may change, through exact inverse
patches. Existing Python/JSON, scientific receipts and older runtimes stay intact.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

from research_layout import relocation_original_bytes
from organize_research_1044_1062 import text_edits

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'research_cognition_physics'
SOURCE = BASE / 'archive_1063_'
TARGET = BASE / 'archive_1063_1085'
MIG = BASE / '_migration/closure_1063_1085_20261010'
MANIFEST = MIG / 'manifest.json'
RUNTIMES = [f'scripts/{name}.py' for name in
            ('research_layout', 'run_research_current', 'run_research_closed',
             'run_research_recent', 'run_research_live', 'test_research_layout')]
SKIP = {'.git', '.agents', '.codex', 'node_modules', '__pycache__'}


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def relative(path):
    return path.relative_to(ROOT).as_posix()


def guard_targets():
    workspace = ROOT.resolve()
    if SOURCE.parent.resolve() != TARGET.parent.resolve() or SOURCE.parent.resolve() != BASE.resolve():
        raise ValueError('Source and target must share the research parent')
    for path in (SOURCE, TARGET, MIG):
        if not path.resolve().is_relative_to(workspace) or path.resolve() in (workspace, BASE.resolve()):
            raise ValueError('Unsafe migration target: ' + str(path))
    if SOURCE.name != 'archive_1063_' or TARGET.name != 'archive_1063_1085':
        raise ValueError('Unexpected stage names')


def changed_markdown(path, raw):
    if path.suffix != '.md':
        return raw, []
    encoding = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    before = raw.decode(encoding)
    after = re.sub(r'archive_1063_(?!\d)', 'archive_1063_1085', before)
    after = after.replace('scripts/run_research_live.py', 'scripts/run_research_spatial_archive.py')
    edits = text_edits(before, after)
    for edit in edits:
        edit['kind'] = 'stage_directory_or_replay_entry_reference'
    return after.encode(encoding), edits


def plan():
    guard_targets()
    if not SOURCE.is_dir() or TARGET.exists() or MANIFEST.exists():
        raise ValueError('Source missing, destination occupied, or manifest already exists')
    all_stage = sorted(p for p in SOURCE.rglob('*') if p.is_file())
    if any(p.is_symlink() or getattr(p, 'is_junction', lambda: False)() for p in [SOURCE, *SOURCE.rglob('*')]):
        raise ValueError('Reparse points are not allowed in the stage tree')
    if any(p.suffix in ('.zip', '.pyc', '.pyo') for p in all_stage):
        raise ValueError('Unexpected archive or bytecode in source evidence')
    entries = []
    external = []
    for path in ROOT.rglob('*.md'):
        if path.is_relative_to(SOURCE) or any(part in SKIP for part in path.parts):
            continue
        raw = path.read_bytes()
        changed, _ = changed_markdown(path, raw)
        if changed != raw:
            external.append(path)
    for path in sorted(all_stage + external):
        raw = path.read_bytes()
        updated, edits = changed_markdown(path, raw)
        destination = TARGET / path.relative_to(SOURCE) if path.is_relative_to(SOURCE) else path
        entry = dict(original=relative(path), destination=relative(destination),
                     original_sha256=digest(raw), current_sha256=digest(updated),
                     original_bytes=len(raw), current_bytes=len(updated), edits=edits)
        if relocation_original_bytes(updated, entry) != raw:
            raise ValueError('Inverse patch mismatch: ' + str(path))
        entries.append(entry)
    old_evidence = [BASE / '_migration/layout_1044_1062_20261009' / name
                    for name in ('manifest.json', 'checks.json', 'replay_checks.json')]
    return dict(schema='spatial_stage_closure_v1', date='2026-10-10',
                source=relative(SOURCE), destination=relative(TARGET), entries=entries,
                unchanged_runtime_sha256={p: digest((ROOT / p).read_bytes()) for p in RUNTIMES},
                unchanged_evidence_sha256={relative(p): digest(p.read_bytes()) for p in old_evidence},
                source_file_count=len(all_stage), source_suffix_counts=dict(Counter(p.suffix for p in all_stage)),
                external_markdown_files=len(external), new_scientific_groups=0,
                cumulative_scientific_calibrations=3860, stage_scientific_calibrations=22,
                scientific_python_json_unchanged=True, old_directory_shell=False,
                directory_copy_used=False, zip_created=False)


def write_exclusive(path, value):
    with path.open('x', encoding='utf8', newline='\n') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2)
        stream.write('\n')


def apply():
    data = plan()
    if MIG.exists():
        raise ValueError('Maintenance destination already exists')
    # Check the complete original inventory immediately before any mutation.
    for entry in data['entries']:
        if digest((ROOT / entry['original']).read_bytes()) != entry['original_sha256']:
            raise ValueError('Source changed after planning: ' + entry['original'])
    MIG.mkdir(parents=True)
    write_exclusive(MANIFEST, data)
    SOURCE.rename(TARGET)
    for entry in data['entries']:
        if not entry['edits']:
            continue
        path = ROOT / entry['destination']
        raw = path.read_bytes()
        if digest(raw) != entry['original_sha256']:
            raise ValueError('File changed during migration: ' + str(path))
        updated, edits = changed_markdown(path, raw)
        if edits != entry['edits'] or digest(updated) != entry['current_sha256']:
            raise ValueError('Planned patch differs: ' + str(path))
        path.write_bytes(updated)
    print(json.dumps(summary(data) | {'renamed': True}, ensure_ascii=False))


def summary(data):
    return {k: v for k, v in data.items() if k not in ('entries', 'unchanged_runtime_sha256', 'unchanged_evidence_sha256')} | {
        'manifest_entries': len(data['entries']),
        'reversible_markdown_files': sum(bool(e['edits']) for e in data['entries'])}


def verify():
    guard_targets()
    data = json.loads(MANIFEST.read_text(encoding='utf8'))
    for entry in data['entries']:
        raw = (ROOT / entry['destination']).read_bytes()
        original = relocation_original_bytes(raw, entry)
        if Path(entry['destination']).suffix in ('.py', '.json') and raw != original:
            raise ValueError('Scientific bytes changed')
    for field in ('unchanged_runtime_sha256', 'unchanged_evidence_sha256'):
        for name, expected in data[field].items():
            if digest((ROOT / name).read_bytes()) != expected:
                raise ValueError('Previously frozen support changed: ' + name)
    inventory = {relative(p) for p in TARGET.rglob('*') if p.is_file()}
    expected = {e['destination'] for e in data['entries'] if e['original'].startswith(relative(SOURCE) + '/')}
    if SOURCE.exists() or inventory != expected:
        raise ValueError('Final directory inventory mismatch')
    print(json.dumps(summary(data) | {'all_byte_and_inventory_checks_passed': True}, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--plan', action='store_true')
    actions.add_argument('--apply', action='store_true')
    actions.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.apply:
        apply()
    elif args.verify:
        verify()
    else:
        print(json.dumps(summary(plan()), ensure_ascii=False))


if __name__ == '__main__':
    main()
