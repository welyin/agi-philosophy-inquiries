"""One-time, byte-preserving closure of rounds 1086--1105.

This maintenance tool is separate from historical scientific programs. It moves
the directory once, adds the authorized closure heading to six navigation files,
and records reversible Markdown path edits. Existing Python/JSON stay unchanged.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
from organize_research_1044_1062 import text_edits
from research_layout import relocation_original_bytes

BASE = ROOT / 'research_cognition_physics'
OLD, NEW = 'archive_1086_', 'archive_1086_1105'
SOURCE, TARGET = BASE / OLD, BASE / NEW
MAINTENANCE = Path(__file__).resolve().parent
SNAPSHOT = MAINTENANCE / 'preparation_snapshot.json'
MANIFEST = MAINTENANCE / 'manifest.json'
PUBLICATION_RECEIPT = 'research_cognition_physics/' + NEW + '/_shared/stage_paper_closure/publication_checks.json'
NAVIGATION = {BASE / n for n in ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')}
NAVIGATION |= {SOURCE / 'README.md', SOURCE / '文件索引.md'}
SKIP = {'.git', '.agents', '.codex', 'node_modules', '__pycache__', '_history', '_migration', '.research_runtime'}
PREFIX = '## 1086—1105阶段结项'


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def relative(path):
    return path.relative_to(ROOT).as_posix()


def guard():
    workspace = ROOT.resolve()
    for path in (SOURCE, TARGET, MAINTENANCE):
        if not path.resolve().is_relative_to(workspace) or path.resolve() in (workspace, BASE.resolve()):
            raise ValueError('Unsafe maintenance path: ' + str(path))
    if SOURCE.parent.resolve() != BASE.resolve() or TARGET.parent.resolve() != BASE.resolve():
        raise ValueError('Archive targets must be direct research children')
    if SOURCE.name != OLD or TARGET.name != NEW:
        raise ValueError('Unexpected archive names')


def verify_snapshot():
    snapshot = json.loads(SNAPSHOT.read_text('utf8'))
    for entry in snapshot['entries']:
        actual = ROOT / entry['path']
        if not actual.is_file() or sha(actual.read_bytes()) != entry['sha256']:
            raise ValueError('Pre-closure asset changed concurrently: ' + entry['path'])
    return snapshot


def decode(raw):
    encoding = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    return raw.decode(encoding), encoding


def changed_markdown(path, raw, template, paper_name, closure_name):
    if path.suffix.lower() != '.md':
        return raw, []
    before, encoding = decode(raw)
    after = re.sub(r'archive_1086_(?!\d)', NEW, before)
    after = after.replace('scripts/run_research_active.py', 'scripts/run_research_sr_archive.py')
    if path in NAVIGATION:
        if PREFIX in before:
            raise ValueError('Closure heading already exists: ' + str(path))
        local = path.parent == SOURCE
        links = {
            '{{PAPER}}': ('../' if local else '') + paper_name,
            '{{STAGE}}': 'README.md' if local else NEW + '/README.md',
            '{{CLOSURE}}': ('_shared/stage_paper_closure/' if local else NEW + '/_shared/stage_paper_closure/') + closure_name,
        }
        section = template.strip()
        for token, value in links.items():
            section = section.replace(token, value)
        if re.search(r'\{\{[^}]+\}\}', section):
            raise ValueError('Unresolved navigation token')
        if not section.startswith(PREFIX):
            raise ValueError('Unexpected closure navigation title')
        newline = '\r\n' if '\r\n' in before else '\n'
        section = section.replace('\r\n', '\n').replace('\n', newline)
        front = section + newline * 2 + '---' + newline * 2 + '## 历史导航（截至1105，保留原研究时点）' + newline * 2
        # Keep the original document title first, then the new closure notice.
        match = re.match(r'(#[^\r\n]*(?:\r\n|\n))', after)
        if match:
            after = match.group(1) + newline + front + after[match.end():].lstrip('\r\n')
        else:
            after = front + after
    edits = text_edits(before, after)
    for edit in edits:
        edit['kind'] = 'authorized_closure_navigation' if path in NAVIGATION else 'stage_directory_or_replay_path_only'
    return after.encode(encoding), edits


def plan(paper_name, closure_name):
    guard()
    if not SOURCE.is_dir() or TARGET.exists() or MANIFEST.exists():
        raise ValueError('Source missing, target occupied, or migration already applied')
    snapshot = verify_snapshot()
    if Path(paper_name).name != paper_name or not (BASE / paper_name).is_file():
        raise ValueError('Expected the completed stage paper in the research root')
    if Path(closure_name).name != closure_name:
        raise ValueError('Closure target must be a file name')
    closure = SOURCE / '_shared/stage_paper_closure'
    if not (closure / closure_name).is_file():
        raise ValueError('Closure anchor is missing')
    template_path = closure / 'navigation_section.md'
    template = template_path.read_text('utf-8-sig')
    stage = sorted(p for p in SOURCE.rglob('*') if p.is_file())
    if any(p.is_symlink() or getattr(p, 'is_junction', lambda: False)() for p in [SOURCE, *SOURCE.rglob('*')]):
        raise ValueError('Reparse points are not allowed in the stage tree')
    if any(p.suffix.lower() in ('.zip', '.pyc', '.pyo') for p in stage):
        raise ValueError('Unexpected archive or bytecode in stage')
    external = []
    for path in ROOT.rglob('*.md'):
        if path.is_relative_to(SOURCE) or any(part in SKIP for part in path.parts):
            continue
        raw = path.read_bytes()
        updated, _ = changed_markdown(path, raw, template, paper_name, closure_name)
        if updated != raw:
            external.append(path)
    entries = []
    for path in sorted(stage + external):
        raw = path.read_bytes()
        updated, edits = changed_markdown(path, raw, template, paper_name, closure_name)
        destination = TARGET / path.relative_to(SOURCE) if path.is_relative_to(SOURCE) else path
        entry = dict(original=relative(path), destination=relative(destination),
                     original_sha256=sha(raw), current_sha256=sha(updated),
                     original_bytes=len(raw), current_bytes=len(updated), edits=edits)
        if relocation_original_bytes(updated, entry) != raw:
            raise ValueError('Inverse edit failed: ' + relative(path))
        entries.append(entry)
    support = {e['path']: e['sha256'] for e in snapshot['entries']
               if e['path'].startswith('scripts/') or e['path'].startswith('research_cognition_physics/_migration/active_documents/')}
    return dict(schema='sr_stage_closure_v1', source=relative(SOURCE), destination=relative(TARGET),
                date='2026-10-10', entries=entries, unchanged_support_sha256=support,
                preparation_snapshot_sha256=sha(SNAPSHOT.read_bytes()),
                source_file_count=len(stage), source_suffix_counts=dict(Counter(p.suffix for p in stage)),
                external_markdown_files=len(external), stage_paper=relative(BASE / paper_name),
                closure_anchor=relative(TARGET / '_shared/stage_paper_closure' / closure_name),
                stage_scope='axiomatized conditional kinematics; complete physical SR and roadmap remain unfinished',
                actual_goal_status='paused', new_research_rounds=0, scientific_calibration_increment=0,
                post_publication_files=[PUBLICATION_RECEIPT],
                scientific_python_json_unchanged=True, previous_active_layer_unchanged=True,
                old_directory_shell=False, directory_copy_used=False, zip_created=False)


def summary(data):
    return {k: v for k, v in data.items() if k not in ('entries', 'unchanged_support_sha256')} | {
        'manifest_entries': len(data['entries']),
        'reversible_markdown_files': sum(bool(e['edits']) for e in data['entries'])}


def apply_edits(raw, edits):
    text, encoding = decode(raw)
    for edit in reversed(edits):
        if text[edit['start']:edit['end']] != edit['before']:
            raise ValueError('Forward patch mismatch')
        text = text[:edit['start']] + edit['after'] + text[edit['end']:]
    return text.encode(encoding)


def apply(paper_name, closure_name):
    data = plan(paper_name, closure_name)
    # Validate all final paths and bytes immediately before the sole directory move.
    guard()
    for entry in data['entries']:
        if sha((ROOT / entry['original']).read_bytes()) != entry['original_sha256']:
            raise ValueError('Source changed after planning: ' + entry['original'])
    with MANIFEST.open('x', encoding='utf8', newline='\n') as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
    SOURCE.rename(TARGET)
    for entry in data['entries']:
        if not entry['edits']:
            continue
        path = ROOT / entry['destination']
        raw = path.read_bytes()
        if sha(raw) != entry['original_sha256']:
            raise ValueError('Concurrent edit during move: ' + entry['destination'])
        updated = apply_edits(raw, entry['edits'])
        if sha(updated) != entry['current_sha256']:
            raise ValueError('Forward patch hash mismatch')
        path.write_bytes(updated)
    print(json.dumps(summary(data) | {'renamed': True}, ensure_ascii=False))


def verify():
    guard()
    data = json.loads(MANIFEST.read_text('utf8'))
    expected = set()
    for entry in data['entries']:
        raw = (ROOT / entry['destination']).read_bytes()
        prior = relocation_original_bytes(raw, entry)
        if Path(entry['destination']).suffix.lower() in ('.py', '.json') and raw != prior:
            raise ValueError('Scientific bytes changed')
        if entry['original'].startswith(relative(SOURCE) + '/'):
            expected.add(entry['destination'])
    for name, digest in data['unchanged_support_sha256'].items():
        if sha((ROOT / name).read_bytes()) != digest:
            raise ValueError('Frozen support changed: ' + name)
    allowed = data.get('post_publication_files', [])
    if allowed != [PUBLICATION_RECEIPT] or PUBLICATION_RECEIPT in expected:
        raise ValueError('Unexpected post-publication allowance')
    inventory = {relative(p) for p in TARGET.rglob('*') if p.is_file()}
    if SOURCE.exists() or not expected <= inventory or inventory - expected - set(allowed):
        raise ValueError('Stage inventory mismatch')
    additions = {name: sha((ROOT / name).read_bytes()) for name in allowed if name in inventory}
    print(json.dumps(summary(data) | {'all_byte_and_inventory_checks_passed': True,
                                     'post_publication_receipt_sha256': additions}, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--plan', action='store_true')
    action.add_argument('--apply', action='store_true')
    action.add_argument('--verify', action='store_true')
    parser.add_argument('--paper')
    parser.add_argument('--closure', default='README.md')
    args = parser.parse_args()
    if args.verify:
        verify()
    elif not args.paper:
        parser.error('--paper is required for plan/apply')
    elif args.apply:
        apply(args.paper, args.closure)
    else:
        print(json.dumps(summary(plan(args.paper, args.closure)), ensure_ascii=False))


if __name__ == '__main__':
    main()
