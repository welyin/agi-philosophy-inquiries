"""Archive the completed thermal-geometry stage with reversible Markdown edits."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import subprocess

from organize_research_1044_1062 import text_edits
from research_layout import relocation_original_bytes

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'research_cognition_physics'
SOURCE = BASE / 'archive_1219_'
TARGET = BASE / 'archive_1219_1233'
MIG = BASE / '_migration/closure_1219_1233'
MANIFEST = MIG / 'manifest.json'
RUNNER = 'scripts/run_research_thermal_geometry_closed.py'
NAV = {'README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md'}
LINK = re.compile(r'(?P<lead>\[[^\]\n]*\]\()(?P<target>[^)\n]+)(?P<tail>\))')


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def rel(path):
    return path.relative_to(ROOT).as_posix()


def mapped(path):
    return TARGET / path.relative_to(SOURCE) if path.is_relative_to(SOURCE) else path


def rename_text(text):
    return re.sub(r'archive_1219_(?!\d)', 'archive_1219_1233', text)


def notice(prefix):
    return ('## 2026-10-15：1219—1233归档与广义相对论阶段准备\n\n'
            '热力学—几何对应阶段已结项，完整目录为`archive_1219_1233`；'
            '15轮报告、代码、结果和原始验收保留，当前复算使用归档入口。'
            '当前没有进行中的科学轮次；下一轮编号为1234，本次未开启新的应用目标。\n\n'
            '**下一阶段方向：广义相对论的认知操作重建——从完整热流与几何熵到共同几何响应。** '
            '优先研究几何响应的选择条件，在同一有效域中连接完整来源、局部平衡和独立几何熵，'
            '再调用成熟重建定理；判清Newton常数、宇宙学项和高阶修正中的剩余自由。'
            '复用既有热学和空间成果，不把制造某台通用加速测温设备作为整个阶段的先决门槛。'
            '本次没有新增公理或宣称Einstein动力学已经推出。\n\n'
            f'[归档核验与复算]({prefix}_migration/closure_1219_1233/README.md) · '
            f'[下一阶段准备]({prefix}_migration/closure_1219_1233/next_stage_plan.md)。'
            '下方旧“当前／下一步”属于历史时点。\n\n')


def rewrite(path, raw):
    if path.suffix.lower() != '.md':
        return raw, []
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf-8'
    before = raw.decode(enc)
    after = LINK.sub(lambda m: m['lead'] + rename_text(m['target']) + m['tail'], before)
    if path.parent == BASE and path.name in NAV:
        head, rest = after.split('\n', 1)
        after = head + '\n\n' + notice('') + rest.lstrip('\r\n')
    if path.parent in (SOURCE, TARGET) and path.name == 'README.md':
        _, rest = after.split('\n', 1)
        after = '# 1219—1233：热力学—几何对应（研究档案）\n\n' + notice('../') + rest.lstrip('\r\n')
    if path.parent in (SOURCE, TARGET) and path.name == '文件索引.md':
        head, rest = after.split('\n', 1)
        after = head + '\n\n归档后的复算统一使用[只读入口](../_migration/closure_1219_1233/README.md)；下方旧命令属于原时点。\n\n' + rest.lstrip('\r\n')
    edits = text_edits(before, after)
    for edit in edits:
        edit['kind'] = 'archive_link_or_current_closure_navigation'
    return after.encode(enc), edits


def guards():
    for p in (SOURCE, TARGET, MIG):
        if not p.resolve().is_relative_to(ROOT.resolve()) or p.resolve() in (ROOT.resolve(), BASE.resolve()):
            raise ValueError('Unsafe target: ' + str(p))
    if SOURCE.parent.resolve() != BASE.resolve() or TARGET.parent.resolve() != BASE.resolve():
        raise ValueError('Only the named sibling archive may move')


def plan():
    guards()
    if MANIFEST.exists() or TARGET.exists() or not SOURCE.is_dir():
        raise ValueError('Unexpected source, target or existing manifest')
    receipt = json.loads((SOURCE/'_shared/goal_completion_receipt.json').read_text(encoding='utf-8'))
    if receipt['goal']['status'] != 'complete':
        raise ValueError('Stage goal has not completed')
    moved = []
    for p in [SOURCE, *SOURCE.rglob('*')]:
        if p.is_symlink() or getattr(p, 'is_junction', lambda: False)():
            raise ValueError('Reparse point: ' + str(p))
        if p.is_file():
            moved.append(p)
    if any(p.suffix.lower() in {'.zip', '.pyc', '.pyo'} for p in moved):
        raise ValueError('Unexpected archive or bytecode')
    found = subprocess.run(['rg', '-l', '--hidden', '-g', '*.md', '-g', '!**/.git/**',
        '-g', '!**/.research_runtime/**', '-g', '!**/closure_1219_1233/**',
        'archive_1219_', '.'], cwd=ROOT, capture_output=True, text=True, encoding='utf-8', check=True)
    candidates = {ROOT/name.removeprefix('.\\').removeprefix('./') for name in found.stdout.splitlines()}
    candidates.update(BASE/name for name in NAV)
    external = [p for p in candidates if p not in moved and rewrite(p, p.read_bytes())[0] != p.read_bytes()]
    entries = []
    for p in sorted(set(moved + external)):
        raw = p.read_bytes()
        current, edits = rewrite(p, raw)
        entry = dict(original=rel(p), destination=rel(mapped(p)), original_sha256=sha(raw),
            current_sha256=sha(current), original_bytes=len(raw), current_bytes=len(current), edits=edits)
        if relocation_original_bytes(current, entry) != raw:
            raise ValueError('Inverse edit failed: ' + str(p))
        entries.append(entry)
    support = [p for p in (ROOT/'scripts').glob('*.py') if p.name not in {Path(__file__).name, Path(RUNNER).name}]
    support += [p for p in (BASE/'_migration').rglob('*') if p.is_file() and not p.is_relative_to(MIG)
                and p.suffix.lower() in {'.py', '.json'}]
    return dict(schema='research_closure_1219_1233_v1', date='2026-10-15',
        stages=[dict(source=rel(SOURCE), destination=rel(TARGET), first_round=1219, last_round=1233)],
        entries=entries, unchanged_support_sha256={rel(p):sha(p.read_bytes()) for p in sorted(set(support))},
        source_file_counts={'archive_1219_':len(moved)}, source_suffix_counts=dict(Counter(p.suffix for p in moved)),
        external_markdown_files=len(external), scientific_python_json_unchanged=True,
        new_scientific_rounds=0, GR_goal_completed=False, thermal_geometry_goal_status='complete',
        next_round=1234, new_goal_started=False, directory_copy_used=False, old_directory_shell=False,
        zip_created=False, replay_entry=RUNNER)


def summary(data):
    return {k:v for k,v in data.items() if k not in {'entries','unchanged_support_sha256'}} | dict(
        manifest_entries=len(data['entries']), reversible_markdown_files=sum(bool(e['edits']) for e in data['entries']),
        frozen_support_files=len(data['unchanged_support_sha256']))


def apply():
    data = plan()
    for e in data['entries']:
        if sha((ROOT/e['original']).read_bytes()) != e['original_sha256']:
            raise ValueError('Concurrent modification: '+e['original'])
    MIG.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
    renamed, patched = False, []
    try:
        SOURCE.rename(TARGET)
        renamed = True
        for e in data['entries']:
            if not e['edits']:
                continue
            p = ROOT/e['destination']
            raw = p.read_bytes()
            if sha(raw) != e['original_sha256']:
                raise ValueError('Concurrent write: '+str(p))
            current, edits = rewrite(p, raw)
            if edits != e['edits'] or sha(current) != e['current_sha256']:
                raise ValueError('Plan differs: '+str(p))
            p.write_bytes(current)
            patched.append((p,raw,current))
    except Exception:
        for p,raw,current in reversed(patched):
            if p.read_bytes() != current:
                raise RuntimeError('Concurrent change blocks rollback: '+str(p))
            p.write_bytes(raw)
        if renamed:
            TARGET.rename(SOURCE)
        MANIFEST.unlink()
        raise
    print(json.dumps(summary(data)|{'status':'MOVED'}, ensure_ascii=False))


def verify():
    guards()
    data = json.loads(MANIFEST.read_text(encoding='utf-8'))
    for e in data['entries']:
        raw = (ROOT/e['destination']).read_bytes()
        old = relocation_original_bytes(raw,e)
        if Path(e['destination']).suffix.lower() in {'.py','.json'} and raw != old:
            raise ValueError('Scientific bytes changed')
    for name,expected in data['unchanged_support_sha256'].items():
        if sha((ROOT/name).read_bytes()) != expected:
            raise ValueError('Frozen support changed: '+name)
    expected = {e['destination'] for e in data['entries'] if e['original'].startswith(rel(SOURCE)+'/')}
    actual = {rel(p) for p in TARGET.rglob('*') if p.is_file()}
    if SOURCE.exists() or actual != expected:
        raise ValueError('Archive inventory mismatch')
    print(json.dumps(summary(data)|{'status':'PASS','old_directories_absent':True}, ensure_ascii=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group()
    mode.add_argument('--apply', action='store_true')
    mode.add_argument('--verify', action='store_true')
    args = p.parse_args()
    if args.apply:
        apply()
    elif args.verify:
        verify()
    else:
        print(json.dumps(summary(plan()), ensure_ascii=False))
