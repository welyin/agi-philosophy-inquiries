"""One-time byte-preserving rename and reversible navigation edits.

No extraction, directory copy, ZIP, junction, or old-directory shell is used.
--apply exclusively records the closure manifest and renames the real folder.
--verify checks original bytes, links and the closed-stage replay interfaces.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

from organize_research_231_775 import links
from research_layout import relocation_original_bytes

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'research_cognition_physics'
SOURCE = BASE/'archive_1009_'
TARGET = BASE/'archive_1009_1043'
MIG = BASE/'_migration/closure_1009_1043_20261008'
MANIFEST = MIG/'manifest.json'
CHECKS = MIG/'checks.json'
NAV = {'README.md', '文件索引.md'}
RUNTIME = ['scripts/research_layout.py', 'scripts/run_research_current.py',
           'scripts/test_research_layout.py']


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def relative(path):
    return path.relative_to(ROOT).as_posix()


def assert_targets():
    project = ROOT.resolve()
    assert SOURCE.parent.resolve() == TARGET.parent.resolve() == BASE.resolve()
    assert SOURCE.resolve().is_relative_to(project)
    assert TARGET.resolve().is_relative_to(project)
    assert SOURCE.name == 'archive_1009_' and TARGET.name == 'archive_1009_1043'
    assert SOURCE != TARGET


def navigation_change(raw, name):
    encoding = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    text = raw.decode(encoding)
    first = text.splitlines(keepends=True)[0]
    newline = '\r\n' if first.endswith('\r\n') else '\n'
    heading = '# 1009—1043轮：条件生成链与剩余自由' if name == 'README.md' else '# 1009—1043轮：文件索引'
    old_heading = first.rstrip('\r\n')
    addition = (newline+'> 本阶段已归档至1043轮。下文的历史轮次与命令保留原时点语义；'
                '冻结脚本统一通过 `scripts/run_research_closed.py --script '
                'archive_1009_1043/1043/verify_round1043.py` 只读复算。'
                '[阶段论文](../条件生成链与物理选择的剩余自由_阶段论文.md) · '
                '[下一阶段准备](../archive_1044_/README.md) · '
                '[目录迁移及核验](../_migration/closure_1009_1043_20261008/README.md)。'
                '此顶部说明优先于下文旧轮次的“当前下一项”。'+newline)
    edits = [dict(start=0, end=len(old_heading), before=old_heading, after=heading,
                  kind='closed_stage_heading'),
             dict(start=len(first), end=len(first), before='', after=addition,
                  kind='closed_stage_replay_notice')]
    pieces, previous = [], 0
    for edit in edits:
        pieces.extend((text[previous:edit['start']], edit['after']))
        previous = edit['end']
    pieces.append(text[previous:])
    return ''.join(pieces).encode(encoding), edits


def plan():
    assert_targets()
    assert SOURCE.is_dir() and not SOURCE.is_symlink()
    assert not TARGET.exists(), TARGET
    entries = []
    for path in sorted(p for p in SOURCE.rglob('*') if p.is_file()):
        assert not path.is_symlink()
        assert path.suffix.lower() not in ('.zip', '.pyc', '.pyo'), path
        rel = path.relative_to(SOURCE)
        raw = path.read_bytes()
        new, edits = navigation_change(raw, rel.as_posix()) if rel.as_posix() in NAV else (raw, [])
        entries.append(dict(original=relative(path), destination=relative(TARGET/rel),
                            original_sha256=sha(raw), current_sha256=sha(new),
                            original_bytes=len(raw), current_bytes=len(new), edits=edits))
    return dict(schema='byte_preserving_stage_closure_v1', date='2026-10-08',
                source=relative(SOURCE), destination=relative(TARGET), entries=entries,
                unchanged_runtime_sha256={p:sha((ROOT/p).read_bytes()) for p in RUNTIME},
                frozen_1043_receipt_sha256=sha((SOURCE/'1043/research_round_1043_checks.json').read_bytes()),
                scientific_python_json_unchanged=True, new_scientific_groups=0,
                cumulative_scientific_groups_unchanged=3819,
                copied_directory=False, old_directory_shell=False, zip_created=False,
                root_navigation_owned_by_mainline=True)


def write_exclusive(path, value):
    with path.open('x', encoding='utf8') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2)
        stream.write('\n')


def apply():
    data = plan()
    assert not MANIFEST.exists() and not MIG.exists(), MIG
    MIG.mkdir(parents=True)
    write_exclusive(MANIFEST, data)
    readme = '''# 1009—1043阶段关闭与目录迁移

2026-10-08。仅将实际目录`archive_1009_`改名为`archive_1009_1043`。新增科学组0，累计3819不变。全部科学Python、JSON及冻结正文原字节保留；仅阶段README和索引加入归档标题／复算说明，其改动以可逆字符补丁保存，不保存整份导航副本。

[逐文件映射与补丁](manifest.json)登记源、目标、前后哈希及原字节数；核验入口是项目根目录的`python -B -X utf8 scripts/close_research_1009_1043.py --verify`。首次执行核验后本目录的`checks.json`保存执行证据。

冻结的`research_layout.py`、`run_research_current.py`与其测试原字节保持。新入口复用这些旧工具，仅加入此关闭阶段的路径、runpy及已知Python子进程路由。不会生成旧目录壳、junction、ZIP或第二份研究目录。历史逻辑路径只存在于只读内存视图中。

## 复算

从项目根目录运行：

```powershell
python -B -X utf8 scripts/run_research_closed.py --verify-closure
python -B -X utf8 scripts/run_research_closed.py --script archive_1009_1043/1043/verify_round1043.py
python -B -X utf8 scripts/run_research_closed.py --script archive_1009_1043/1043/legacy_1009_scope_check.py
```

旧逻辑路径亦可作为`--script`值。直接执行历史脚本仍可能按旧目录查依赖，应使用上述入口。1009原默认全目录比较失败继续保留，范围复核通过不改称原入口通过。目录迁移不补充物理证明或改变阶段结论。

只有两份导航的反向补丁用于验证原字节；它们以后若再编辑，须记录新的维护层，不能绕过旧哈希检查。新增阶段资料应保存在后继目录中。
'''
    with (MIG/'README.md').open('x', encoding='utf8') as stream:
        stream.write(readme)
    # Absolute, same-parent targets have been checked above. One filesystem
    # rename moves the actual folder; there is no recursive shell operation.
    SOURCE.rename(TARGET)
    for entry in data['entries']:
        if entry['edits']:
            path = ROOT/entry['destination']
            old = path.read_bytes()
            assert sha(old) == entry['original_sha256']
            new, edits = navigation_change(old, path.name)
            assert edits == entry['edits'] and sha(new) == entry['current_sha256']
            path.write_bytes(new)
    print(json.dumps(dict(renamed=True, source=data['source'], destination=data['destination'],
                          files=len(data['entries']), modified_markdown_files=2,
                          scientific_python_json_unchanged=True), ensure_ascii=False))


def verify(record):
    assert_targets()
    assert not SOURCE.exists() and TARGET.is_dir()
    data = json.loads(MANIFEST.read_text('utf8'))
    counts = dict(files=0, byte_identical=0, reversible_markdown=0, python=0, json=0)
    for entry in data['entries']:
        path = ROOT/entry['destination']
        raw = path.read_bytes()
        original = relocation_original_bytes(raw, entry)
        counts['files'] += 1
        counts['byte_identical'] += raw == original
        counts['reversible_markdown'] += raw != original
        if path.suffix in ('.py', '.json'):
            assert raw == original
            counts['python' if path.suffix == '.py' else 'json'] += 1
    assert {relative(p) for p in TARGET.rglob('*') if p.is_file()} == {e['destination'] for e in data['entries']}
    for path, digest in data['unchanged_runtime_sha256'].items():
        assert sha((ROOT/path).read_bytes()) == digest
    assert sha((TARGET/'1043/research_round_1043_checks.json').read_bytes()) == data['frozen_1043_receipt_sha256']
    missing, number = [], 0
    for path in sorted(TARGET.rglob('*.md')):
        for _, _, target, local in links(path.read_text('utf-8-sig')):
            number += 1
            actual = (path.parent/local.replace('\\', '/')).resolve()
            if not actual.exists():
                missing.append(dict(source=relative(path), target=target))
    assert not missing, missing
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1', PYTHONUTF8='1')
    commands = [
        [sys.executable, '-B', '-X', 'utf8', str(ROOT/'scripts/run_research_closed.py'), '--verify-closure'],
        [sys.executable, '-B', '-X', 'utf8', str(ROOT/'scripts/run_research_closed.py'), '--script', 'archive_1009_1043/1043/verify_round1043.py'],
        [sys.executable, '-B', '-X', 'utf8', str(ROOT/'scripts/run_research_closed.py'), '--script', 'archive_1009_1043/1043/legacy_1009_scope_check.py'],
    ]
    runs = []
    for command in commands:
        run = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True, encoding='utf8', timeout=180)
        runs.append(dict(command=command, exit_code=run.returncode, stdout=run.stdout, stderr=run.stderr))
    result = dict(schema='closed_stage_migration_checks_v1', date='2026-10-08',
                  all_checks_passed=all(r['exit_code'] == 0 for r in runs), counts=counts,
                  local_markdown_links_checked=number, missing_local_links=missing,
                  manifest_sha256=sha(MANIFEST.read_bytes()),
                  frozen_1043_receipt_sha256=data['frozen_1043_receipt_sha256'],
                  unchanged_runtime_sha256=data['unchanged_runtime_sha256'],
                  maintenance_scripts_sha256={relative(p):sha(p.read_bytes()) for p in
                                             (Path(__file__), ROOT/'scripts/run_research_closed.py')},
                  default_1009_failure_preserved=True, original_33_successes_not_reclassified=True,
                  new_scientific_groups=0, cumulative_unchanged=3819,
                  old_directory_absent=True, zip_created=False, duplicate_directory_created=False,
                  runs=runs)
    if record:
        assert result['all_checks_passed'], result
        write_exclusive(CHECKS, result)
    elif CHECKS.exists():
        old = json.loads(CHECKS.read_text('utf8'))
        for key in ('counts', 'manifest_sha256', 'frozen_1043_receipt_sha256',
                    'unchanged_runtime_sha256', 'maintenance_scripts_sha256'):
            assert result[key] == old[key], key
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if not result['all_checks_passed']:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--plan', action='store_true')
    modes.add_argument('--apply', action='store_true')
    modes.add_argument('--verify', action='store_true')
    parser.add_argument('--record', action='store_true', help='Create only the maintenance check receipt')
    args = parser.parse_args()
    if args.apply:
        apply()
    elif args.verify:
        verify(args.record)
    else:
        data = plan()
        print(json.dumps({k:v for k,v in data.items() if k != 'entries'} | {'files':len(data['entries'])}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
