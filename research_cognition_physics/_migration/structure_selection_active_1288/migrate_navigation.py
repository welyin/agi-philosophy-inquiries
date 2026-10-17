"""Record one reversible navigation prefix after the 1281--1287 closure."""
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
BASE = ROOT / 'research_cognition_physics'
MANIFEST = BASE / '_migration/closure_1281_1287/manifest.json'
CLOSED_RUNNER = ROOT / 'scripts/run_research_matter_gauge_closed.py'
NAV = ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')
STAGE = 'research_cognition_physics/archive_1288_'
sys.path.insert(0, str(ROOT / 'scripts'))
from research_layout import relocation_original_bytes


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    if (HERE / 'layer.json').exists():
        raise FileExistsError('Navigation layer already recorded')
    closure = json.loads(MANIFEST.read_text(encoding='utf-8'))
    boundary = {e['destination']: e for e in closure['entries']}
    addition = ('\n## 当前目标：共同操作要求下的规范结构与物质选择\n\n'
        '2026-10-18接续已归档的1281—1287；新阶段从1288开始，'
        '见[阶段导航](archive_1288_/README.md)、[正式目标](archive_1288_/_shared/goal.md)和'
        '[复用与选择计划](archive_1288_/_shared/reuse_and_selection_plan.md)。\n\n'
        '研究局域量子接续、动态内部参考、完整记录与资源／几何来源的共同要求，'
        '怎样约束规范群、物质表示／粒子种类和相互作用。先固定最小共同任务、复用成熟结构结论，'
        '仅对影响选择的分支投入共同实现；不预押标准模型群唯一，不把单实例失败扩大为群的排除。'
        '整体公理压缩后置，防循环、相容性和输入登记同步进行。\n\n'
        '**启动状态：1288进入查重与命题准入，尚无本阶段新增科学验收。** '
        '本次归档与目标登记不计科学轮次，不自动采用新公理。下方历次启动状态属于历史记录，'
        '当前安排以本节及新阶段导航为准。\n')
    changes = []
    for filename in NAV:
        p = BASE / filename
        raw = p.read_bytes()
        name = p.relative_to(ROOT).as_posix()
        assert sha(raw) == boundary[name]['current_sha256'], name
        encoding = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
        text = raw.decode(encoding)
        at = text.index('\n') + 1
        current = (text[:at] + addition + text[at:]).encode(encoding)
        entry = dict(original=name, destination=name, original_sha256=sha(raw),
                     current_sha256=sha(current), original_bytes=len(raw), current_bytes=len(current),
                     edits=[dict(start=at, end=at, before='', after=addition)])
        assert relocation_original_bytes(current, entry) == raw
        changes.append((p, raw, current, entry))
    support = dict(closure['unchanged_support_sha256'])
    for p in (MANIFEST, CLOSED_RUNNER):
        support[p.relative_to(ROOT).as_posix()] = sha(p.read_bytes())
    for name, expected in support.items():
        assert sha((ROOT / name).read_bytes()) == expected, name
    entry_points = [STAGE + '/README.md', STAGE + '/_shared/goal.md']
    assert all((ROOT / p).is_file() for p in entry_points)
    layer = dict(schema='structure_selection_active_navigation_v1', date='2026-10-18',
                 stage_path=STAGE, stage_first_round=1288, entry_points=entry_points,
                 startup_state='1288 admission and reuse; no new scientific acceptance',
                 entries=[change[3] for change in changes], unchanged_support_sha256=support)
    # Preflight all sources before the first mutation.
    assert all(p.read_bytes() == raw for p, raw, _, _ in changes)
    for p, _, current, _ in changes:
        p.write_bytes(current)
    with (HERE / 'layer.json').open('x', encoding='utf-8', newline='\n') as f:
        json.dump(layer, f, ensure_ascii=False, indent=2)
        f.write('\n')
    print(json.dumps(dict(status='NAVIGATION_LAYER_CREATED', files=4,
                         layer_sha256=sha((HERE / 'layer.json').read_bytes()),
                         scientific_files_modified=False)))


if __name__ == '__main__':
    main()
