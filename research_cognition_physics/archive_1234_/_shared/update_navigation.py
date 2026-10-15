"""Publish a reversible 1234 notice above the closed 1219--1233 boundary."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research_cognition_physics'
HERE = Path(__file__).resolve().parent
NOTICE = HERE / 'navigation_notice.md'
LAYER = BASE / '_migration/gr_response_active_1234/layer.json'
FROZEN = BASE / '_migration/closure_1219_1233/manifest.json'
NAMES = ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')
sys.path.insert(0, str(ROOT / 'scripts'))
from research_layout import relocation_original_bytes


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--apply', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    frozen_raw = FROZEN.read_bytes()
    frozen = json.loads(frozen_raw)
    require(frozen['schema'] == 'research_closure_1219_1233_v1', 'Unexpected closure schema')
    boundary = {entry['destination']: entry for entry in frozen['entries']}
    notice_raw = NOTICE.read_bytes()
    notice = notice_raw.decode('utf-8-sig').replace('\r\n', '\n').strip()
    require(notice.startswith('## ') and notice.count('$$') % 2 == 0, 'Invalid active notice')
    old_raw = LAYER.read_bytes() if LAYER.exists() else None
    old = json.loads(old_raw) if old_raw else None
    if old:
        require(old['schema'] == 'gr_response_active_navigation_v1', 'Unexpected navigation schema')
    existing = {entry['original']: entry for entry in old['entries']} if old else {}
    keys = {(BASE / name).relative_to(ROOT).as_posix() for name in NAMES}
    if old:
        require(set(existing) == keys and len(old['entries']) == 4, 'Incomplete existing navigation')
    entries, changes = [], []
    for name in NAMES:
        path = BASE / name
        key = path.relative_to(ROOT).as_posix()
        current = path.read_bytes()
        original = relocation_original_bytes(current, existing[key]) if old else current
        require((sha(original), len(original)) == (boundary[key]['current_sha256'], boundary[key]['current_bytes']),
                'Current navigation does not restore frozen closure: ' + name)
        enc = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf-8'
        text = original.decode(enc)
        require(text.startswith('# ') and '\n' in text, 'Invalid document heading: ' + name)
        newline = '\r\n' if '\r\n' in text else '\n'
        split = text.index('\n') + 1
        insertion = newline + notice.replace('\n', newline) + newline + newline
        after = (text[:split] + insertion + text[split:]).encode(enc)
        entry = dict(original=key, destination=key, original_sha256=sha(original), current_sha256=sha(after),
                     original_bytes=len(original), current_bytes=len(after),
                     edits=[dict(start=split, end=split, before='', after=insertion,
                                 kind='active_stage_navigation_notice')])
        require(relocation_original_bytes(after, entry) == original, 'Navigation cannot be reversed')
        entries.append(entry)
        changes.append((path, current, after))
    supports = ['scripts/run_research_thermal_geometry_closed.py', FROZEN.relative_to(ROOT).as_posix(),
                'research_cognition_physics/_migration/thermal_geometry_active_1219/layer.json',
                'research_cognition_physics/_migration/closure_1219_1233/runtime_checks.json']
    support = {name: sha((ROOT / name).read_bytes()) for name in supports}
    if old:
        require(old['unchanged_support_sha256'] == support, 'Frozen navigation support changed')
    layer = dict(schema='gr_response_active_navigation_v1', date='2026-10-15',
                 notice_path=NOTICE.relative_to(ROOT).as_posix(), notice_sha256=sha(notice_raw),
                 entries=entries, unchanged_support_sha256=support)
    if args.apply:
        require(FROZEN.read_bytes() == frozen_raw and NOTICE.read_bytes() == notice_raw, 'Concurrent source change')
        require((LAYER.read_bytes() if LAYER.exists() else None) == old_raw, 'Concurrent layer change')
        require(all(path.read_bytes() == current for path, current, after in changes), 'Concurrent navigation change')
        LAYER.parent.mkdir(parents=True, exist_ok=True)
        written = []
        try:
            for path, current, after in changes:
                if current != after:
                    path.write_bytes(after)
                    written.append((path, current, after))
            LAYER.write_text(json.dumps(layer, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        except Exception:
            for path, current, after in reversed(written):
                if path.read_bytes() == after:
                    path.write_bytes(current)
            if old_raw is not None:
                LAYER.write_bytes(old_raw)
            raise
    print(json.dumps(dict(status='APPLIED' if args.apply else 'CHECKED', files=4,
                          changed_files=sum(current != after for path, current, after in changes),
                          closed_boundary_unchanged=True, original_navigation_recoverable=True,
                          scientific_files_modified=False), ensure_ascii=False))


if __name__ == '__main__':
    main()
