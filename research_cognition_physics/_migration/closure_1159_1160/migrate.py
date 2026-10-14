"""Prepare/apply link-only metadata for the separately authorized directory rename.

The physical Rename-Item is deliberately performed in PowerShell after explicit
absolute-path checks. This script neither moves directories nor copies evidence.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research_cognition_physics'
OLD = 'research_cognition_physics/archive_1159_'
NEW = 'research_cognition_physics/archive_1159_1160'
HERE = Path(__file__).resolve().parent
MANIFEST = HERE / 'manifest.json'
sys.path.insert(0, str(ROOT / 'scripts'))
from research_layout import relocation_original_bytes


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def prepare():
    if MANIFEST.exists() or not (ROOT / OLD).is_dir() or (ROOT / NEW).exists():
        raise ValueError('Expected original directory and no destination/manifest')
    completion_path = ROOT / OLD / '_shared/completion_receipt.json'
    completion_raw = completion_path.read_bytes()
    completion = json.loads(completion_raw)
    assets_raw = (ROOT / OLD / '_shared/stage_assets_check.json').read_bytes()
    if completion['goal']['status'] != 'complete' or completion['stage_assets_sha256'] != sha(assets_raw):
        raise ValueError('Completed baseline receipt is not bound to its assets')
    for name, expected in json.loads(assets_raw)['sha256'].items():
        if sha((ROOT / name).read_bytes()) != expected:
            raise ValueError('Completed baseline asset changed: ' + name)
    stage_files = sorted(p for p in (ROOT / OLD).rglob('*') if p.is_file())
    external = [BASE / name for name in ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')]
    external += [BASE / '07_局域相对论过程的区域收支与量子记录_阶段论文.md']
    current_links = set(external + [ROOT / OLD / '_shared/navigation_notice.md'])
    entries = []
    for path in [*stage_files, *external]:
        old = path.relative_to(ROOT).as_posix()
        new = NEW + old[len(OLD):] if old.startswith(OLD + '/') else old
        raw = path.read_bytes()
        edits = []
        encoding = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf-8'
        if path in current_links:
            text = raw.decode(encoding)
            edits = [dict(start=m.start(), end=m.end(), before=m.group(),
                          after='archive_1159_1160', kind='closed_stage_link')
                     for m in re.finditer(r'archive_1159_(?=[/\\])', text)]
        elif old in (OLD + '/README.md', OLD + '/文件索引.md'):
            text = raw.decode(encoding)
            replacements = {'# 1159起：': '# 1159—1160：', '# 1159起阶段文件索引': '# 1159—1160阶段文件索引'}
            if path.name == 'README.md':
                before = '根部导航使用本阶段活动可逆层维护，正文源为[_shared/navigation_notice.md](_shared/navigation_notice.md)，通过[_shared/update_navigation.py](_shared/update_navigation.py)的`--check`与`--apply`维护。历史复算入口为 scripts/run_research_local_dynamics_active.py；它先恢复此前导航视图，再调用已有只读兼容层。'
                after = '本目录已按完成轮次归档为archive_1159_1160。旧导航正文与维护器保留为历史证据，不再直接运行旧维护器。当前只读复算入口为 scripts/run_research_gr_active.py，支持本目录的新路径与历史archive_1159_路径，并恢复各脚本原来的路径和文档视图；迁移说明见[可逆归档记录](../_migration/closure_1159_1160/README.md)。'
                if before not in text:
                    raise ValueError('Expected completed-stage README instruction not found')
                replacements[before] = after
            else:
                replacements['- [可逆导航维护脚本](_shared/update_navigation.py)：预览用`--check`，发布用`--apply`；根导航不手工重复改写。'] = '- [旧可逆导航维护脚本](_shared/update_navigation.py)：仅保留历史证据，不再直接发布。当前复算与导航接入见[可逆归档记录](../_migration/closure_1159_1160/README.md)。'
            for before, after in replacements.items():
                for match in re.finditer(re.escape(before), text):
                    edits.append(dict(start=match.start(), end=match.end(), before=before,
                                      after=after, kind='closed_stage_navigation'))
        new_raw = raw
        if edits:
            text = raw.decode(encoding)
            edits.sort(key=lambda e: e['start'])
            for edit in reversed(edits):
                text = text[:edit['start']] + edit['after'] + text[edit['end']:]
            new_raw = text.encode(encoding)
        entry = dict(original=old, destination=new, original_sha256=sha(raw),
                     current_sha256=sha(new_raw), original_bytes=len(raw),
                     current_bytes=len(new_raw), edits=edits)
        if relocation_original_bytes(new_raw, entry) != raw:
            raise ValueError('Inverse verification failed: ' + old)
        entries.append(entry)
    support_paths = set(ROOT.joinpath('scripts').glob('run_research*.py'))
    support_paths.add(ROOT / 'scripts/research_layout.py')
    support_paths.update(p for p in (BASE / '_migration').rglob('*.json')
                         if not p.is_relative_to(HERE) and 'gr_active_1161' not in p.parts)
    support_paths.discard(ROOT / 'scripts/run_research_gr_active.py')
    support = {p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in sorted(support_paths)}
    data = dict(schema='research_closure_1159_1160_v1', date=datetime.now().isoformat(timespec='seconds'),
                source=OLD, destination=NEW, entries=entries,
                baseline_completion_sha256=sha(completion_raw), baseline_assets_sha256=sha(assets_raw),
                baseline_goal_status='complete', scientific_rounds=[1159, 1160], scientific_check_groups=8,
                unchanged_support_sha256=support, created_stage_files=[],
                historical_embedded_paths='Frozen PY/JSON, goal text and review commands retain their recorded paths; replay maps these paths without refreshing evidence.')
    with MANIFEST.open('x', encoding='utf-8') as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
    print(json.dumps(dict(status='PREPARED', entries=len(entries), moved=len(stage_files),
                          reversible_markdown_files=sum(bool(e['edits']) for e in entries),
                          frozen_support_files=len(support)), ensure_ascii=False))


def apply():
    data = json.loads(MANIFEST.read_text(encoding='utf-8'))
    if (ROOT / OLD).exists() or not (ROOT / NEW).is_dir():
        raise ValueError('Physical rename must be completed first')
    pending = []
    for entry in data['entries']:
        path = ROOT / entry['destination']
        raw = path.read_bytes()
        if sha(raw) != entry['original_sha256'] or len(raw) != entry['original_bytes']:
            raise ValueError('Unexpected post-rename bytes: ' + entry['destination'])
        if entry['edits']:
            encoding = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf-8'
            text = raw.decode(encoding)
            for edit in reversed(entry['edits']):
                if text[edit['start']:edit['end']] != edit['before']:
                    raise ValueError('Original edit mismatch')
                text = text[:edit['start']] + edit['after'] + text[edit['end']:]
            new = text.encode(encoding)
            if relocation_original_bytes(new, entry) != raw:
                raise ValueError('Inverse mismatch')
            pending.append((path, raw, new))
    for name, expected in data['unchanged_support_sha256'].items():
        if sha((ROOT / name).read_bytes()) != expected:
            raise ValueError('Frozen support changed: ' + name)
    for path, raw, new in pending:
        if path.read_bytes() != raw:
            raise ValueError('Concurrent edit: ' + str(path))
        path.write_bytes(new)
    print(json.dumps(dict(status='APPLIED', changed_markdown=len(pending), scientific_bytes_modified=False)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'apply'))
    prepare() if parser.parse_args().action == 'prepare' else apply()
