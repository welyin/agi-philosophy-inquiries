"""Read-only replay after reclassifying the 1323--1326 paper as interim work.

Current documents stay on disk.  Only manifest-listed Markdown is restored in
memory for historical verification; no old paper, ZIP, or backup is created.
"""
from __future__ import annotations

import argparse
import builtins
from contextlib import contextmanager
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import sys
from urllib.parse import unquote

from research_layout import ROOT, relocation_original_bytes

MANIFEST = ROOT / 'research_cognition_physics/_migration/reclassify_1326_interim/manifest.json'
NAVIGATION_START = ROOT / 'research_cognition_physics/_migration/long_horizon_start_1327/manifest.json'
CURRENT_PROGRESS = NAVIGATION_START.with_name('current_progress_manifest.json')
PROGRESS_BEGIN = '<!-- CURRENT_RESEARCH_PROGRESS_BEGIN -->'
PROGRESS_END = '<!-- CURRENT_RESEARCH_PROGRESS_END -->'
NAVIGATION_NAMES = {'research_cognition_physics/' + name for name in (
    'README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md', 'archive_1323_/README.md')}
OLD_PAPER = 'research_cognition_physics/16_未知参考交接的作用选择与动态来源_阶段论文.md'
NEW_PAPER = 'research_cognition_physics/archive_1323_/未知参考交接的作用选择与动态来源_中间成果.md'


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def safe_name(value):
    if not isinstance(value, str) or '\\' in value:
        raise ValueError('Expected a project-relative POSIX path')
    path = PurePosixPath(value)
    if path.is_absolute() or '..' in path.parts or path.as_posix() != value or ':' in value:
        raise ValueError('Unsafe project-relative path: ' + value)
    if not (ROOT / value).resolve().is_relative_to(ROOT):
        raise ValueError('Path leaves the workspace: ' + value)
    return value


def path_key(path):
    if not isinstance(path, (str, bytes, os.PathLike)):
        return None
    return os.path.normcase(os.path.abspath(os.fsdecode(path)))


class CurrentProgress:
    """One replaceable progress block; no growing history of navigation copies."""
    def __init__(self, navigation_raw):
        self.raw = CURRENT_PROGRESS.read_bytes() if CURRENT_PROGRESS.exists() else None
        self.entries, self.evidence = {}, {}
        if self.raw is None:
            return
        data = json.loads(self.raw)
        if data.get('schema') != 'current_research_progress_v1' or navigation_raw is None:
            raise ValueError('Unexpected current-progress schema or missing frozen navigation baseline')
        if data['navigation_start_manifest_sha256'] != digest(navigation_raw):
            raise ValueError('Progress manifest is not bound to the frozen navigation-start layer')
        for entry in data['entries']:
            name = safe_name(entry['path'])
            if name in self.entries:
                raise ValueError('Duplicate current-progress navigation path')
            self.entries[name] = entry
        if set(self.entries) != NAVIGATION_NAMES:
            raise ValueError('Progress may only change the five declared navigation blocks')
        if not isinstance(data['evidence_sha256'], dict) or not data['evidence_sha256']:
            raise ValueError('Progress must bind actual evidence files')
        for name, expected in data['evidence_sha256'].items():
            safe_name(name)
            if not name.startswith('research_cognition_physics/archive_1323_/'):
                raise ValueError('Progress evidence must belong to the active research archive')
            if not re.fullmatch(r'[0-9a-f]{64}', expected):
                raise ValueError('Invalid progress evidence digest: ' + name)
            raw = (ROOT / name).read_bytes()
            if digest(raw) != expected:
                raise ValueError('Current-progress evidence changed: ' + name)
            self.evidence[name] = raw

    def strip(self, name, raw):
        if name not in self.entries:
            return raw
        entry = self.entries[name]
        if digest(raw) != entry['current_sha256'] or len(raw) != entry['current_bytes']:
            raise ValueError('Current progress navigation changed: ' + name)
        enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
        text = raw.decode(enc)
        newline = '\r\n' if '\r\n' in text else '\n'
        cut = text.index('\n') + 1
        if text.count(PROGRESS_BEGIN) != 1 or text.count(PROGRESS_END) != 1:
            raise ValueError('Progress markers must be unique: ' + name)
        opening, ending = newline + PROGRESS_BEGIN + newline, PROGRESS_END + newline
        if not text[cut:].startswith(opening):
            raise ValueError('Progress must immediately follow the first title line: ' + name)
        stop = text.index(ending, cut) + len(ending)
        block = text[cut:stop].encode('utf8')
        if digest(block) != entry['block_sha256'] or len(block) != entry['block_bytes']:
            raise ValueError('Progress block changed: ' + name)
        previous = (text[:cut] + text[stop:]).encode(enc)
        if digest(previous) != entry['baseline_sha256'] or len(previous) != entry['baseline_bytes']:
            raise ValueError('Content outside progress block changed: ' + name)
        return previous

    def verify_unchanged(self):
        if self.raw is not None and CURRENT_PROGRESS.read_bytes() != self.raw:
            raise ValueError('Current-progress manifest changed during verification')
        for name, raw in self.evidence.items():
            if (ROOT / name).read_bytes() != raw:
                raise ValueError('Current-progress evidence changed during verification: ' + name)


class InterimReclassification:
    def __init__(self):
        self.manifest_raw = MANIFEST.read_bytes()
        self.data = json.loads(self.manifest_raw)
        if self.data.get('schema') != 'interim_reclassification_v1':
            raise ValueError('Unexpected interim reclassification schema')
        self.navigation_raw = NAVIGATION_START.read_bytes() if NAVIGATION_START.exists() else None
        self.navigation = json.loads(self.navigation_raw) if self.navigation_raw is not None else None
        self.navigation_entries = {}
        if self.navigation is not None:
            if self.navigation.get('schema') != 'long_horizon_navigation_start_v1':
                raise ValueError('Unexpected navigation-start schema')
            if self.navigation['previous_manifest_sha256'] != digest(self.manifest_raw):
                raise ValueError('Navigation start is not bound to the unchanged earlier migration')
            for entry in self.navigation['entries']:
                name = safe_name(entry['destination'])
                if name != safe_name(entry['original']) or name in self.navigation_entries:
                    raise ValueError('Navigation start may only add one layer to each existing path')
                self.navigation_entries[name] = entry
            if set(self.navigation_entries) != NAVIGATION_NAMES:
                raise ValueError('Exactly the five current navigation files must be bound')
            receipt_path = ROOT / safe_name(self.navigation['goal_receipt'])
            receipt_raw = receipt_path.read_bytes()
            if digest(receipt_raw) != self.navigation['goal_receipt_sha256']:
                raise ValueError('Goal-start receipt changed')
            if json.loads(receipt_raw)['tool_result']['goal']['status'] != 'active':
                raise ValueError('Navigation start lacks an active goal-start receipt')
        self.progress = CurrentProgress(self.navigation_raw)
        self.cache, self.current, self.entries = {}, {}, {}
        destinations = set()
        for entry in self.data['entries']:
            old, new = safe_name(entry['original']), safe_name(entry['destination'])
            if not old.endswith('.md') or not new.endswith('.md'):
                raise ValueError('Only Markdown can be overlaid')
            if old in self.entries or new in destinations:
                raise ValueError('Duplicate original or destination')
            if new != old and (old, new) != (OLD_PAPER, NEW_PAPER):
                raise ValueError('Unexpected moved document: ' + old)
            raw = (ROOT / new).read_bytes()
            baseline = self.progress.strip(new, raw)
            previous = relocation_original_bytes(baseline, self.navigation_entries[new]) if new in self.navigation_entries else baseline
            restored = relocation_original_bytes(previous, entry)
            self.entries[old] = entry
            destinations.add(new)
            self.current[new] = raw
            self.cache[path_key(ROOT / old)] = restored
        if OLD_PAPER not in self.entries or self.entries[OLD_PAPER]['destination'] != NEW_PAPER:
            raise ValueError('The withdrawn paper must have its declared interim destination')
        if (ROOT / OLD_PAPER).exists():
            raise ValueError('The withdrawn numbered paper still exists on disk')
        self.protected = set(self.cache) | {path_key(ROOT / p) for p in destinations}
        self.scientific = self.data['scientific_sha256']
        if not isinstance(self.scientific, dict) or len(self.scientific) != 40:
            raise ValueError('Exactly forty frozen scientific files must be bound')
        for name, expected in self.scientific.items():
            safe_name(name)
            if not re.fullmatch(r'[0-9a-f]{64}', expected):
                raise ValueError('Invalid scientific digest: ' + name)
            if name in self.entries or name in destinations:
                raise ValueError('Scientific evidence cannot be rewritten by the overlay: ' + name)
            if digest((ROOT / name).read_bytes()) != expected:
                raise ValueError('Frozen scientific evidence changed: ' + name)

    def verify(self):
        if MANIFEST.read_bytes() != self.manifest_raw:
            raise ValueError('Reclassification manifest changed during verification')
        if self.navigation_raw is not None and NAVIGATION_START.read_bytes() != self.navigation_raw:
            raise ValueError('Navigation-start manifest changed during verification')
        self.progress.verify_unchanged()
        links = 0
        for name, raw in self.current.items():
            path = ROOT / name
            if path.read_bytes() != raw:
                raise ValueError('Current Markdown changed during verification: ' + name)
            text = raw.decode('utf-8-sig')
            for match in re.finditer(r'!?\[[^\]\n]*\]\((<[^>]+>|[^)\n]+)\)', text):
                target = match.group(1).strip()
                target = target[1:-1] if target.startswith('<') and target.endswith('>') else re.sub(r'\s+"[^"\n]*"$', '', target)
                if re.match(r'^(https?|mailto|app|codex|plugin):', target, re.I):
                    continue
                target = unquote(re.sub(r':\d+$', '', target.partition('#')[0]))
                if not target:
                    continue
                if re.match(r'^/[A-Za-z]:[/\\]', target):
                    target = target[1:]
                linked = Path(target)
                if not linked.is_absolute():
                    linked = path.parent / linked
                if not linked.exists():
                    raise FileNotFoundError('Broken current Markdown link: ' + name + ' -> ' + target)
                links += 1
        return dict(status='PASS', schema=self.data['schema'], manifest_sha256=digest(self.manifest_raw),
                    markdown_files=len(self.entries), local_links=links, frozen_scientific_files=len(self.scientific),
                    current_and_original_hashes_verified=True, withdrawn_paper_absent=True,
                    original_markdown_recovered_in_memory=True, files_written=0,
                    navigation_start_files=len(self.navigation_entries),
                    navigation_start_manifest_sha256=digest(self.navigation_raw) if self.navigation_raw is not None else None,
                    current_progress_files=len(self.progress.entries),
                    current_progress_evidence_files=len(self.progress.evidence),
                    current_progress_manifest_sha256=digest(self.progress.raw) if self.progress.raw is not None else None,
                    scientific_claims_revalidated=False)

    @contextmanager
    def historical_view(self):
        """Patch only the listed paths; prefetching avoids recursive I/O."""
        saved = {name: getattr(Path, name) for name in ('read_bytes', 'read_text', 'exists', 'is_file', 'open')}
        saved_open = builtins.open
        cache, protected = self.cache, self.protected

        def stream_for(file, mode='r', buffering=-1, encoding=None, errors=None, newline=None, closefd=True, opener=None):
            key = path_key(file)
            if key in protected and any(flag in mode for flag in 'wax+'):
                raise PermissionError('Reclassification replay is read-only: ' + str(file))
            if key not in cache:
                return saved_open(file, mode, buffering, encoding, errors, newline, closefd, opener)
            if mode not in ('r', 'rb', 'rt') or opener is not None or not closefd:
                raise ValueError('Unsupported historical read mode or custom opener')
            raw = io.BytesIO(cache[key])
            return raw if 'b' in mode else io.TextIOWrapper(raw, encoding=encoding or 'utf8', errors=errors, newline=newline)

        def read_bytes(path):
            key = path_key(path)
            return cache[key] if key in cache else saved['read_bytes'](path)

        def read_text(path, encoding=None, errors=None, newline=None):
            if path_key(path) not in cache:
                return saved['read_text'](path, encoding=encoding, errors=errors, **({'newline': newline} if newline is not None else {}))
            with stream_for(path, 'r', encoding=encoding, errors=errors, newline=newline) as stream:
                return stream.read()

        def exists(path, **kwargs):
            return True if path_key(path) in cache else saved['exists'](path, **kwargs)

        def is_file(path, **kwargs):
            return True if path_key(path) in cache else saved['is_file'](path, **kwargs)

        def path_open(path, mode='r', buffering=-1, encoding=None, errors=None, newline=None):
            if path_key(path) in protected:
                return stream_for(path, mode, buffering, encoding, errors, newline)
            return saved['open'](path, mode, buffering, encoding, errors, newline)

        try:
            Path.read_bytes, Path.read_text = read_bytes, read_text
            Path.exists, Path.is_file, Path.open = exists, is_file, path_open
            builtins.open = stream_for
            yield
        finally:
            for name, method in saved.items():
                setattr(Path, name, method)
            builtins.open = saved_open


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--verify-reclassification', action='store_true')
    mode.add_argument('--verify-navigation', action='store_true')
    mode.add_argument('--script')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
    if arguments not in ([], ['--check']):
        parser.error('Only default read-only verification or --check is supported; publishing and write-new are forbidden')
    if not args.script and arguments:
        parser.error('Script arguments require --script')
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    view = InterimReclassification()
    result = view.verify()
    if args.verify_reclassification:
        print(json.dumps(result, ensure_ascii=False))
        return
    if args.script:
        script = Path(args.script)
        if not script.is_absolute():
            script = ROOT / script
        script = script.resolve()
        if not script.is_relative_to(ROOT) or script.suffix != '.py' or not script.is_file():
            parser.error('Script must be an existing Python file inside this workspace')
        forwarded = ['--script', script.relative_to(ROOT).as_posix(), '--', *arguments]
    else:
        forwarded = ['--verify-navigation']
    import run_research_interaction_rules_active as previous
    old_argv = sys.argv
    try:
        sys.argv = [str(ROOT / 'scripts/run_research_interaction_rules_active.py'), *forwarded]
        with view.historical_view():
            previous.main()
    finally:
        sys.argv = old_argv


if __name__ == '__main__':
    main()
