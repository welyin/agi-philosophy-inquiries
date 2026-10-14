"""Close two stage directories without rewriting frozen Python/JSON evidence.

Markdown navigation edits are reversible. The manifest is a path/hash map, not
a second archive. Only the exact two user-requested sibling directories move.
"""
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
MIG = BASE / '_migration/closure_1161_1218'
MANIFEST = MIG / 'manifest.json'
PAIRS = {'archive_1161_': 'archive_1161_1200', 'archive_1201_': 'archive_1201_1218'}
RUNNER = 'scripts/run_research_thermodynamics_closed.py'
NAV = {'README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md'}
PAPER = '08_认知热学角色与热力学四律_阶段论文.md'
LINK = re.compile(r'(?P<lead>\[[^\]\n]*\]\()(?P<target>[^)\n]+)(?P<tail>\))')


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def rel(path):
    return path.relative_to(ROOT).as_posix()


def mapped(path):
    for old, new in PAIRS.items():
        if path.is_relative_to(BASE / old):
            return BASE / new / path.relative_to(BASE / old)
    return path


def replace_names(text):
    for old, new in PAIRS.items():
        text = re.sub(re.escape(old) + r'(?!\d)', new, text)
    return text


def closure_paragraph(prefix):
    return ('## 2026-10-15：目录归档与后继交接\n\n'
            '1161—1200和1201—1218已分别归档为`archive_1161_1200`、`archive_1201_1218`。'
            '热力学四律阶段已在声明范围内结项；1161—1200的GR目标保持原有未完成状态。'
            '下一研究轮从1219起，方向为同一实际过程的热力学—几何连接；本次仅归档，尚未启动新目标或新阶段目录。\n\n'
            f'[归档说明与复算入口]({prefix}_migration/closure_1161_1218/README.md)。'
            '历史代码、结果、回执及科学范围保留；旧路径通过只读入口解析到当前唯一文件，文档引用修改可逆核验。\n\n')


def rewrite(path, raw):
    if path.suffix.lower() != '.md':
        return raw, []
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf-8'
    before = raw.decode(enc)
    # Preserve historical prose and commands; update the destinations of links.
    after = LINK.sub(lambda m: m['lead'] + replace_names(m['target']) + m['tail'], before)
    if path.parent == BASE and path.name in NAV:
        head, rest = after.split('\n', 1)
        after = head + '\n\n' + closure_paragraph('') + rest.lstrip('\r\n')
    if path.name == 'README.md' and path.parent.name in set(PAIRS) | set(PAIRS.values()):
        head, rest = after.split('\n', 1)
        if path.parent.name in ('archive_1161_', 'archive_1161_1200'):
            head = '# 1161—1200：广义相对论的认知操作重建（研究档案）'
        after = head + '\n\n' + closure_paragraph('../') + rest.lstrip('\r\n')
    if path.name == '文件索引.md' and path.parent.name in set(PAIRS) | set(PAIRS.values()):
        head, rest = after.split('\n', 1)
        after = head + '\n\n归档后的代码复算与历史哈希核对统一使用[只读入口](../_migration/closure_1161_1218/README.md)。\n\n' + rest.lstrip('\r\n')
    if path == BASE / PAPER:
        old_command = "-B -X utf8 'research_cognition_physics/archive_1201_/_shared/after1218_common_prediction_check.py'"
        new_command = "-B -X utf8 'scripts/run_research_thermodynamics_closed.py' --script 'research_cognition_physics/archive_1201_1218/_shared/after1218_common_prediction_check.py'"
        if old_command not in after:
            raise ValueError('Expected paper replay command absent')
        after = after.replace(old_command, new_command)
    edits = text_edits(before, after)
    for edit in edits:
        edit['kind'] = 'archive_link_or_current_closure_navigation'
    return after.encode(enc), edits


def guards():
    for old, new in PAIRS.items():
        source, target = BASE / old, BASE / new
        if source.parent.resolve() != BASE.resolve() or target.parent.resolve() != BASE.resolve():
            raise ValueError('Unexpected directory parent')
        for path in (source, target, MIG):
            if not path.resolve().is_relative_to(ROOT.resolve()) or path.resolve() in (ROOT.resolve(), BASE.resolve()):
                raise ValueError('Unsafe target: ' + str(path))
    return True


def plan():
    guards()
    if MANIFEST.exists() or any((BASE / new).exists() for new in PAIRS.values()):
        raise ValueError('Manifest or destination already exists')
    moved = []
    for old in PAIRS:
        source = BASE / old
        if not source.is_dir():
            raise FileNotFoundError(source)
        for p in [source, *source.rglob('*')]:
            if p.is_symlink() or getattr(p, 'is_junction', lambda: False)():
                raise ValueError('Reparse point in source tree: ' + str(p))
            if p.is_file():
                moved.append(p)
    if any(p.suffix.lower() in {'.zip', '.pyc', '.pyo'} for p in moved):
        raise ValueError('Unexpected ZIP or bytecode in source evidence')
    found = subprocess.run(['rg', '-l', '--hidden', '-g', '*.md', '-g', '!**/.git/**',
        '-g', '!**/.research_runtime/**', '-g', '!**/closure_1161_1218/**',
        'archive_1161_|archive_1201_', '.'], cwd=ROOT, capture_output=True, text=True, encoding='utf-8', check=True)
    candidates = {ROOT / name.removeprefix('.\\').removeprefix('./') for name in found.stdout.splitlines()}
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
    support = [p for p in (ROOT / 'scripts').glob('*.py') if p.name not in {Path(__file__).name, Path(RUNNER).name}]
    support += [p for p in (BASE / '_migration').rglob('*')
                if p.is_file() and not p.is_relative_to(MIG) and p.suffix.lower() in {'.py', '.json'}]
    return dict(schema='research_closure_1161_1218_v1', date='2026-10-15',
        stages=[dict(source=rel(BASE/old), destination=rel(BASE/new),
                     first_round=1161 if old=='archive_1161_' else 1201,
                     last_round=1200 if old=='archive_1161_' else 1218) for old,new in PAIRS.items()],
        entries=entries, unchanged_support_sha256={rel(p):sha(p.read_bytes()) for p in sorted(set(support))},
        source_file_counts={old:sum(p.is_relative_to(BASE/old) for p in moved) for old in PAIRS},
        source_suffix_counts=dict(Counter(p.suffix for p in moved)), external_markdown_files=len(external),
        scientific_python_json_unchanged=True, new_scientific_rounds=0,
        GR_goal_completed=False, four_laws_goal_status='complete', next_round=1219,
        new_goal_started=False, directory_copy_used=False, old_directory_shell=False, zip_created=False,
        replay_entry=RUNNER)


def summary(data):
    return {k:v for k,v in data.items() if k not in {'entries','unchanged_support_sha256'}} | {
        'manifest_entries':len(data['entries']), 'reversible_markdown_files':sum(bool(e['edits']) for e in data['entries']),
        'frozen_support_files':len(data['unchanged_support_sha256'])}


def apply():
    data = plan()
    for e in data['entries']:
        if sha((ROOT/e['original']).read_bytes()) != e['original_sha256']:
            raise ValueError('Concurrent modification: ' + e['original'])
    MIG.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
    renamed, patched = [], []
    try:
        for old,new in PAIRS.items():
            source,target=BASE/old,BASE/new
            source.rename(target)
            renamed.append((source,target))
        for e in data['entries']:
            if not e['edits']:
                continue
            p=ROOT/e['destination']
            raw=p.read_bytes()
            if sha(raw) != e['original_sha256']:
                raise ValueError('Concurrent write during relocation: '+str(p))
            current,edits=rewrite(p,raw)
            if edits != e['edits'] or sha(current) != e['current_sha256']:
                raise ValueError('Relocation plan differs: '+str(p))
            p.write_bytes(current)
            patched.append((p,raw,current))
    except Exception:
        for p,raw,current in reversed(patched):
            if p.read_bytes() != current:
                raise RuntimeError('Concurrent change blocks safe rollback: '+str(p))
            p.write_bytes(raw)
        for source,target in reversed(renamed):
            target.rename(source)
        MANIFEST.unlink()
        raise
    print(json.dumps(summary(data)|{'status':'MOVED'},ensure_ascii=False))


def verify():
    guards()
    data=json.loads(MANIFEST.read_text(encoding='utf-8'))
    for e in data['entries']:
        raw=(ROOT/e['destination']).read_bytes()
        original=relocation_original_bytes(raw,e)
        if Path(e['destination']).suffix.lower() in {'.py','.json'} and raw != original:
            raise ValueError('Scientific bytes changed')
    for name,expected in data['unchanged_support_sha256'].items():
        if sha((ROOT/name).read_bytes()) != expected:
            raise ValueError('Frozen support changed: '+name)
    for old,new in PAIRS.items():
        expected={e['destination'] for e in data['entries'] if e['original'].startswith(rel(BASE/old)+'/')}
        actual={rel(p) for p in (BASE/new).rglob('*') if p.is_file()}
        if (BASE/old).exists() or actual != expected:
            raise ValueError('Wrong destination inventory: '+new)
    print(json.dumps(summary(data)|{'status':'PASS','old_directories_absent':True},ensure_ascii=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--apply',action='store_true')
    modes.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    if args.apply: apply()
    elif args.verify: verify()
    else: print(json.dumps(summary(plan()),ensure_ascii=False))
