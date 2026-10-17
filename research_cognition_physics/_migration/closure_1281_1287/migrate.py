"""Prepare/verify reversible Markdown edits for the 1281--1287 closure.

The directory rename is performed by native PowerShell after explicit path
validation. This utility never moves/deletes files or rewrites scientific data.
"""
import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research_cognition_physics'
HERE = Path(__file__).resolve().parent
OLD = 'research_cognition_physics/archive_1281_'
NEW = 'research_cognition_physics/archive_1281_1287'
NAV = [BASE / n for n in ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')]
PAPER = BASE / '11_认知比较条件与局域物质场的共同来源_阶段论文.md'
sys.path.insert(0, str(ROOT / 'scripts'))
from research_layout import relocation_original_bytes


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def dump_new(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(value, f, ensure_ascii=False, indent=2)
        f.write('\n')


def edit_bytes(raw, edits):
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    text = raw.decode(enc)
    for e in reversed(edits):
        assert text[e['start']:e['end']] == e['before']
        text = text[:e['start']] + e['after'] + text[e['end']:]
    return text.encode(enc)


def prepare():
    old, new = ROOT / OLD, ROOT / NEW
    assert old.is_dir() and not new.exists()
    before = json.loads((HERE / 'independent_before.json').read_text(encoding='utf-8'))
    files = sorted(p for p in old.rglob('*') if p.is_file())
    frozen = {x['path']: x for x in before['files']}
    assert {p.relative_to(ROOT).as_posix() for p in files} == set(frozen)
    for p in files:
        item = frozen[p.relative_to(ROOT).as_posix()]
        assert sha(p.read_bytes()) == item['sha256'], p
    external = NAV + [PAPER]
    entries = []
    for p in sorted(files + external):
        raw = p.read_bytes()
        name = p.relative_to(ROOT).as_posix()
        dest = NEW + name[len(OLD):] if name.startswith(OLD + '/') else name
        edits = []
        if p.suffix == '.md':
            text = raw.decode('utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8')
            for m in re.finditer(r'archive_1281_(?![0-9])', text):
                edits.append(dict(start=m.start(), end=m.end(), before=m.group(), after='archive_1281_1287'))
            line_end = text.index('\n') + 1
            if p in NAV:
                addition = ('\n## 1281—1287阶段结项与归档\n\n'
                    '已归档为[archive_1281_1287](archive_1281_1287/README.md)，共7轮、7组作者与7组独立核验。'
                    '[论文11](11_认知比较条件与局域物质场的共同来源_阶段论文.md)与'
                    '[最终验收](archive_1281_1287/_shared/stage_acceptance_after1287.json)已完成。'
                    '成果限条件性结构约束及具名共同实例；公理审计保留强输入和剩余自由，不宣称全部物理已由一般认知唯一生成。\n\n'
                    '[目录迁移与复算](_migration/closure_1281_1287/README.md)。下方启动状态属于历史；后继安排以最新阶段导航为准。\n')
                edits.append(dict(start=line_end, end=line_end, before='', after=addition))
            if name == OLD + '/README.md':
                heading = text[:line_end].rstrip('\r\n')
                edits.append(dict(start=0, end=len(heading), before=heading,
                                  after='# 1281—1287：物质场与规范相互作用的认知来源'))
                addition = ('\n## 正式归档\n\n'
                    '2026-10-18按用户要求收尾为 `archive_1281_1287`。7轮的报告、代码、结果和审阅全部保留；'
                    '科学Python／JSON字节不变，Markdown只更新路径与归档导航。'
                    '[迁移证据与当前复算入口](../_migration/closure_1281_1287/README.md)。\n\n'
                    '后继从[archive_1288_](../archive_1288_/README.md)启动，'
                    '以共同操作要求约束规范群、物质表示及相互作用的选择。历史结项时“未开启下一目标”不覆盖这一后继安排。\n')
                edits.append(dict(start=line_end, end=line_end, before='', after=addition))
        edits.sort(key=lambda e: (e['start'], e['end']))
        current = edit_bytes(raw, edits)
        e = dict(original=name, destination=dest, original_sha256=sha(raw), current_sha256=sha(current),
                 original_bytes=len(raw), current_bytes=len(current), edits=edits)
        assert relocation_original_bytes(current, e) == raw
        entries.append(e)
    active_path = BASE / '_migration/matter_gauge_active_1281/layer.json'
    support = dict(json.loads(active_path.read_text(encoding='utf-8'))['unchanged_support_sha256'])
    for p in (active_path, ROOT / 'scripts/run_research_matter_gauge_active.py'):
        support[p.relative_to(ROOT).as_posix()] = sha(p.read_bytes())
    for name, expected in support.items():
        assert sha((ROOT / name).read_bytes()) == expected, name
    plan = dict(schema='research_closure_1281_1287_v1', date='2026-10-18',
                stages=[dict(source=OLD, destination=NEW, first_round=1281, last_round=1287)],
                entries=entries, unchanged_support_sha256=support,
                moved_files=len(files), external_markdown_files=len(external),
                directory_copy_used=False, zip_used=False, scientific_results_rewritten=False)
    dump_new(HERE / 'plan.json', plan)
    print(json.dumps(dict(status='PLAN_READY', moved_files=len(files), entries=len(entries),
                         markdown_edits=sum(bool(e['edits']) for e in entries)), ensure_ascii=False))


def apply():
    assert not (ROOT / OLD).exists() and (ROOT / NEW).is_dir()
    plan = json.loads((HERE / 'plan.json').read_text(encoding='utf-8'))
    for entry in plan['entries']:
        p = ROOT / entry['destination']
        raw = p.read_bytes()
        # Preflight the entire batch before writing, retaining safe re-entry.
        assert sha(raw) in (entry['original_sha256'], entry['current_sha256']), p
    for entry in plan['entries']:
        p = ROOT / entry['destination']
        raw = p.read_bytes()
        if sha(raw) != entry['current_sha256']:
            current = edit_bytes(raw, entry['edits'])
            assert sha(current) == entry['current_sha256']
            p.write_bytes(current)
        relocation_original_bytes(p.read_bytes(), entry)
    dump_new(HERE / 'manifest.json', plan)
    # Retain the final manifest only, not a byte-identical planning duplicate.
    assert (HERE / 'plan.json').read_bytes() == (HERE / 'manifest.json').read_bytes()
    (HERE / 'plan.json').unlink()
    print(json.dumps(dict(status='MIGRATED_AND_REVERSED', manifest_sha256=sha((HERE / 'manifest.json').read_bytes()),
                         moved_files=plan['moved_files'], entries=len(plan['entries'])), ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['prepare', 'apply'])
    args = parser.parse_args()
    prepare() if args.mode == 'prepare' else apply()
