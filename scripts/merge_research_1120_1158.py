"""Move two closed research stages into one archive; retain frozen evidence.

Only Markdown references and current navigation change, with reversible edits.
Scientific Python/JSON and prior runtime/manifests retain their original bytes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

from organize_research_231_775 import links
from organize_research_1044_1062 import text_edits
from research_layout import relocation_original_bytes

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'research_cognition_physics'
NEW = BASE / 'archive_1120_1158'
MIG = BASE / '_migration/merge_1120_1158'
OLD = {'archive_1120_': '1120_1156', 'archive_1157_1158': '1157_1158'}
NAV = {'README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def rel(path):
    return path.relative_to(ROOT).as_posix()


def safe(path):
    path = path.resolve()
    assert path.is_relative_to(ROOT) and path != ROOT, ('Outside workspace', path)
    assert not set(path.relative_to(ROOT).parts) & {'.git', '.codex', '.agents', '.research_runtime'}
    return path


def mapped(path):
    path = Path(os.path.abspath(path))
    for old, group in OLD.items():
        source = BASE / old
        if path == source:
            return NEW
        if path.is_relative_to(source):
            tail = path.relative_to(source)
            if tail.parts[0] == '_shared':
                return NEW / '_shared' / group / Path(*tail.parts[1:])
            if tail.as_posix() in ('README.md', '文件索引.md'):
                return NEW / '_shared' / group / tail
            return NEW / tail
    return path


def plain_paths(text):
    # Keep bare historical stage names/dates intact. Only explicit path prefixes
    # in prose/code are changed; Markdown destinations are handled separately.
    for old, group in OLD.items():
        for sep in ('/', '\\'):
            prefix = old + sep
            target = 'archive_1120_1158' + sep
            for name in ('README.md', '文件索引.md'):
                text = text.replace(prefix + name, target + '_shared' + sep + group + sep + name)
            text = text.replace(prefix + '_shared' + sep, target + '_shared' + sep + group + sep)
            text = text.replace(prefix, target)
    return text


def notice():
    return '''## 1120—1158合并归档

按用户要求，1120—1156与1157—1158现合并为[archive_1120_1158](archive_1120_1158/README.md)：39份编号报告、39轮正式验收、152组构造复算。先记录旧公理下的组合缺口，再接续CC组合公理的充分性与共同材料模型；两段原目标和结项范围分别保留。本次仅整理目录，不新增研究轮次、公理或经验校准。

[统一轮次索引](archive_1120_1158/文件索引.md) · [06 认知组合公理与局域量子过程](06_认知组合公理与局域量子过程_阶段论文.md) · [CC正式采用](archive_1120_1158/_shared/1157_1158/adoption_and_dependencies.md) · [迁移及历史复算](_migration/merge_1120_1158/README.md)。旧目录已实际移除，历史代码通过统一只读入口解析到现存文件。

以下分段栏目保留各自研究时点；目录合并不把1120—1156的限定不足与新增CC后的正向结论混为同一组前提。

---

'''


def rewrite(path, raw):
    if path.suffix.lower() != '.md':
        return raw, []
    encoding = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf-8'
    before = raw.decode(encoding)
    after = before
    destination = mapped(path)
    changes = []
    for start, end, target, local in links(before):
        old_target = Path(os.path.abspath(path.parent / local))
        target_new = mapped(old_target)
        if target_new != old_target or destination.parent != path.parent:
            updated = os.path.relpath(target_new, destination.parent).replace('\\', '/')
            if '#' in target:
                updated += '#' + target.split('#', 1)[1]
            changes.append((start, end, updated))
    for start, end, value in reversed(changes):
        after = after[:start] + value + after[end:]
    after = plain_paths(after)
    if path.parent == BASE and path.name in NAV:
        newline = '\r\n' if '\r\n' in after else '\n'
        split = after.index('\n') + 1
        after = after[:split] + newline + notice().replace('\n', newline) + after[split:]
    return after.encode(encoding), text_edits(before, after)


def inventory():
    paths = set()
    for name in OLD:
        source = BASE / name
        assert source.is_dir(), source
        paths.update(p for p in source.rglob('*') if p.is_file())
    command = ['rg', '--no-ignore', '--files', '-g', '*.md', '-g', '!**/.git/**',
               '-g', '!**/.research_runtime/**', '-g', '!**/node_modules/**',
               '-g', '!**/.codex/**', '-g', '!**/.agents/**']
    for name in subprocess.run(command, cwd=ROOT, check=True, capture_output=True,
                               encoding='utf-8').stdout.splitlines():
        path = ROOT / name
        if path in paths:
            continue
        raw = path.read_bytes()
        if any(old.encode() in raw for old in OLD):
            paths.add(path)
    return sorted(paths)


def readme():
    return '''# 1120—1158：局域实现的组合缺口与认知组合公理

2026-10-12。按用户要求合并1120—1156与1157—1158两段已完成研究。39份编号报告在本目录根部，证明、代码、结果与审阅保留在1120至1158各编号目录中。

## 研究顺序与结论

1. **1120—1156：共同来源、材料接续与局域实现的缺口。** 37轮、144组构造复算。最终给出现行旧公理范围中的限定不足结论：完整量子资源、实际原钟、端点及四P，并不单独保证原事件和期限的因果局域实现。[该段结项](_shared/1120_1156/stage_closure/stage_completion.md)。具体候选失败、某条路线失败与明确公理不足保持区分。
2. **1157—1158：认知组合公理与局域量子过程。** 2轮、8组构造复算。正式采用CC1—CC3，证明其在完整有限母类型和实际原语资格下保证全部许可有限组合的局域实现；共同模型保留纠缠、原钟、完整资源、四P与许可未来。[正式采用](_shared/1157_1158/adoption_and_dependencies.md) · [该段结项](_shared/1157_1158/stage_completion.md)。

合计39轮正式验收、152组构造复算；新增公理一组CC（三条），科学校准增量0、历史累计3860。两个目标均已完成。目录整理不新增轮次或更改科学结论，也不自动启动下一阶段。

## 阅读入口

- [逐轮索引](文件索引.md)
- [06 认知组合公理与局域量子过程阶段论文](../06_认知组合公理与局域量子过程_阶段论文.md)
- [1120—1156原段说明](_shared/1120_1156/README.md)与[原索引](_shared/1120_1156/文件索引.md)
- [1157—1158原段说明](_shared/1157_1158/README.md)与[原索引](_shared/1157_1158/文件索引.md)
- [迁移和只读复算说明](../_migration/merge_1120_1158/README.md)

原共享目标、准入、方向审计、依赖账与结项收据分别位于`_shared/1120_1156/`和`_shared/1157_1158/`，同名文件不覆盖。历史文稿中的旧目录名及“当前／尚未采用”保留对应时点；实际可点击引用已经迁移。

## 复算

从项目根目录使用统一入口，例如：

```powershell
python -B -X utf8 scripts/run_research_combined.py --script research_cognition_physics/archive_1120_1158/1158/check.py
python -B -X utf8 scripts/run_research_combined.py --verify-layout
```

统一入口按每个历史模块的原逻辑路径读取唯一现存文件；对Markdown只在内存中恢复路径修改前的字节，用于原哈希校验。科学代码和JSON保持原字节，不依赖ZIP、复制树或旧目录壳。旧冻结校验程序请通过此入口运行。

结论限各自公开的有限有效范围。具体材料、场论和引力等后续自由仍按原文保留。
'''


def index():
    lines = ['# 1120—1158轮次索引', '', '[阶段说明](README.md) · [研究总目录](../README.md)', '',
             '|轮次|报告|状态|材料|', '|---|---|---|---|']
    for n in range(1120, 1159):
        source = BASE / ('archive_1120_' if n <= 1156 else 'archive_1157_1158')
        title = next(s[2:] for s in (source / f'research_note_{n}.md').read_text(
            encoding='utf-8-sig').splitlines() if s.startswith('# ')).replace('|', '／')
        acceptance = json.loads((source / str(n) / 'acceptance.json').read_text(encoding='utf-8'))
        assert acceptance['status'] == 'PASS', n
        lines.append(f'|{n}|[{title}](research_note_{n}.md)|已验收|[{n}/]({n}/)|')
    lines += ['', '## 两段原始管理记录', '',
              '- [1120—1156依赖与验收账](_shared/1120_1156/dependency_ledger.json)',
              '- [1157—1158依赖与验收账](_shared/1157_1158/dependency_ledger.json)',
              '- [CC采用记录](_shared/1157_1158/adoption_record.json)及[阶段论文](../06_认知组合公理与局域量子过程_阶段论文.md)',
              '- [完整迁移清单](../_migration/merge_1120_1158/manifest.json)', '',
              '39轮、152组复算；本次合并仅维护文件位置与引用，原审阅和科学范围不变。', '']
    return '\n'.join(lines)


def build():
    assert not NEW.exists(), 'Merged archive already exists'
    assert not (MIG / 'manifest.json').exists(), 'Migration already recorded'
    changes, entries, destinations = {}, [], set()
    preserved_links = 0
    for path in inventory():
        source = safe(path)
        destination = safe(mapped(source))
        raw = source.read_bytes()
        after, edits = rewrite(source, raw)
        if source == destination and raw == after:
            continue
        assert destination not in destinations, ('Collision', destination)
        destinations.add(destination)
        assert source == destination or not destination.exists(), ('Would overwrite', destination)
        entry = dict(original=rel(source), destination=rel(destination), original_sha256=sha(raw),
                     current_sha256=sha(after), original_bytes=len(raw), current_bytes=len(after), edits=edits)
        if source.suffix.lower() in ('.pyc', '.pyo') and '__pycache__' in source.parts:
            # Preserve the existing local file, but never make a Git-ignored
            # bytecode cache a prerequisite for a clean checkout's replay.
            entry['optional_local_cache'] = True
        assert relocation_original_bytes(after, entry) == raw
        if source.suffix.lower() != '.md':
            assert raw == after
        if source.suffix.lower() == '.md':
            before_links, after_links = list(links(raw.decode('utf-8-sig'))), list(links(after.decode('utf-8-sig')))
            if source.parent == BASE and source.name in NAV:
                after_links = after_links[len(list(links(notice()))):]
            assert len(before_links) == len(after_links), source
            for a, b in zip(before_links, after_links):
                target = Path(os.path.abspath(source.parent / a[3]))
                if target.exists():
                    actual = Path(os.path.abspath(destination.parent / b[3]))
                    assert actual == mapped(target), (source, a[2], b[2], mapped(target))
                    preserved_links += 1
        entries.append(entry)
        changes[source] = (destination, raw, after)
    created = {NEW / 'README.md': readme().encode('utf-8'), NEW / '文件索引.md': index().encode('utf-8')}
    support = {rel(p): sha(p.read_bytes()) for p in (ROOT / 'scripts').glob('*.py')
               if p.name not in {Path(__file__).name, 'run_research_combined.py'}}
    support[rel(BASE / '_migration/field_active_1120/layer.json')] = sha((BASE / '_migration/field_active_1120/layer.json').read_bytes())
    support[rel(BASE / '_migration/merge_1106_1119/manifest.json')] = sha((BASE / '_migration/merge_1106_1119/manifest.json').read_bytes())
    manifest = dict(schema='research_merge_1120_1158_v1', date='2026-10-12', entries=entries,
                    created_files=[dict(path=rel(p), sha256=sha(b), bytes=len(b)) for p, b in created.items()],
                    unchanged_support_sha256=support, old_stages=[rel(BASE / x) for x in OLD],
                    merged_stage=rel(NEW), formal_rounds=list(range(1120, 1159)), numerical_groups=152,
                    scientific_python_json_byte_identical=True, directory_copy_used=False, zip_used=False,
                    preserved_local_link_targets=preserved_links,
                    ancillary_cache_files=[e['destination'] for e in entries if e['destination'].endswith('.pyc')])
    return manifest, changes, created


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    manifest, changes, created = build()
    print(json.dumps(dict(entries=len(changes), moved_files=sum(p != d for p, (d, _, _) in changes.items()),
                         changed_markdown=sum(old != new for _, old, new in changes.values()),
                         preserved_local_link_targets=manifest['preserved_local_link_targets'],
                         frozen_python_json=sum(p.suffix in ('.py', '.json') for p in changes),
                         ancillary_cache_files=len(manifest['ancillary_cache_files']),
                         formal_rounds=39, numerical_groups=152), ensure_ascii=False))
    if not args.apply:
        return
    # All absolute sources and destinations checked before moving any file.
    for source, (destination, raw, _) in changes.items():
        assert source.read_bytes() == raw, ('Concurrent edit', source)
        safe(source)
        safe(destination)
        assert source == destination or not destination.exists()
    for path in created:
        safe(path)
        assert not path.exists()
    MIG.mkdir(parents=True, exist_ok=True)
    with (MIG / 'manifest.json').open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(manifest, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
    for source, (destination, raw, after) in changes.items():
        destination.parent.mkdir(parents=True, exist_ok=True)
        if source != destination:
            source.rename(destination)
        if raw != after:
            destination.write_bytes(after)
    for path, raw in created.items():
        path.write_bytes(raw)
    # Remove only verified empty original directories; never recursively delete.
    for old in OLD:
        directory = safe(BASE / old)
        for path in sorted(directory.rglob('*'), key=lambda p: len(p.parts), reverse=True):
            assert path.is_dir(), ('Unmoved file', path)
            safe(path).rmdir()
        directory.rmdir()
    print('APPLIED: actual files moved, original directories absent, no copies or ZIPs.')


if __name__ == '__main__':
    main()
