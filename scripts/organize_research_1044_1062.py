"""Archive the 19 existing rounds; preserve evidence and reversible navigation edits.

This is maintenance, not scientific completion. No ZIPs, copied trees, symlinks,
or historical directory shells are created. --apply refuses existing targets.
"""
from __future__ import annotations

import argparse
from collections import Counter
from difflib import SequenceMatcher
import hashlib
import json
from pathlib import Path
import re

from organize_research_231_775 import links
from research_layout import relocation_original_bytes

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'research_cognition_physics'
MIG = BASE / '_migration/layout_1044_1062_20261009'
PAIRS = {'archive_1044_': 'archive_1044_1045', 'archive_1046_': 'archive_1046_1062'}
RUNTIME = ('scripts/research_layout.py', 'scripts/run_research_current.py',
           'scripts/run_research_closed.py', 'scripts/test_research_layout.py')
# Two source expressions use mathematical evaluation syntax, not Markdown links.
MATH_NOT_LINKS = {('1044/finite_window_proof.md', 'x−X'), ('1051/proof.md', 'σ')}
MARKERS = {
    'README.md': '## 当前入口：1044—1045阶段已结项',
    'research_direction.md': '以下保存上一阶段结项及历史方向，不覆盖当前目标。',
    'RESEARCH_STATE.md': '以下为1044—1045结项快照和更早历史，不覆盖本栏。',
}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def rel(path):
    return path.relative_to(ROOT).as_posix()


def mapped(path):
    for old, new in PAIRS.items():
        if path.is_relative_to(BASE / old):
            return BASE / new / path.relative_to(BASE / old)
    return path


def paths_rewritten(text):
    for old, new in PAIRS.items():
        text = re.sub(re.escape(old) + r'(?!\d)', new, text)
    return text


def text_edits(before, after):
    if before == after:
        return []
    # Line matching avoids quadratic character matching in repeated long logs.
    old_lines, new_lines = before.splitlines(keepends=True), after.splitlines(keepends=True)
    old_offsets, new_offsets = [0], [0]
    for line in old_lines:
        old_offsets.append(old_offsets[-1] + len(line))
    for line in new_lines:
        new_offsets.append(new_offsets[-1] + len(line))
    return [dict(start=old_offsets[a], end=old_offsets[b],
                 before=before[old_offsets[a]:old_offsets[b]],
                 after=after[new_offsets[c]:new_offsets[d]], kind='navigation_or_archive_reference')
            for tag, a, b, c, d in SequenceMatcher(None, old_lines, new_lines, autojunk=False).get_opcodes()
            if tag != 'equal']


SUMMARIES = {
    1044: '实际散射、目标反冲与未知后态共同约束下，仍保可辨耦合自由。',
    1045: '实际热写后态接续同源任务，精确分类必须保留的记录资料。',
    1046: '同一接收器由实际光子波包形成完整方向效果。',
    1047: '束缚源连接Newton与TT出口；遗漏势场应力不能由单一源标定修复。',
    1048: '微观偶极与有限记录材料的匹配条件及误差。',
    1049: '固定校准下，重top阈值同时运输散射与Higgs来源。',
    1050: '真实Coulomb过程连接自然后态与完整初始电荷仪器。',
    1051: '有限预算动作商提供有效Lie接口，保微观非光滑的可能。',
    1052: '同一弱流的两个能标接入公开资料与共同预测检验。',
    1053: '1050实际Coulomb后态进一步接入电子Newton来源字典。',
    1054: '同一光—原子作用共同给方向与局部位置读数。',
    1055: '真实相对受力、有限质量参考反冲接入同一光学读数。',
    1056: '同一弱衰变中，电子粗记录不足以恢复总角能流来源任务。',
    1057: '有限并行归一条件约束近似因果次序。',
    1058: '共享准备的四份一致性进一步约束因果可分近似。',
    1059: '固定公开α标定合同出现张力，不能只修改磁矩一行修复。',
    1060: '指定全族完整混合反常消去仍不足以选择标准超荷。',
    1061: '核与介子弱流的CKM保守外包非空，尚非完整联合拟合。',
    1062: '自由初片下同应力可保绝对质量差，实验室有限读口仍可辨。',
}


def round_table(start, end):
    out = ['|轮次|主结果与范围|材料|', '|---|---|---|']
    for n in range(start, end + 1):
        out.append(f'|[{n}](research_note_{n}.md)|{SUMMARIES[n]}|[{n}/]({n}/)|')
    return '\n'.join(out)


def stage_readme(old):
    early = old == 'archive_1044_'
    start, end = (1044, 1045) if early else (1046, 1062)
    title = '共同物理过程中的耦合自由与记录资料' if early else '共同物理接口、公开检验与整体证据综合'
    status = ('原有限阶段目标已结项；两轮科学校准，累计至3821。'
              if early else '整体路线的进展归档，未结项；17轮科学校准，累计至3838。研究目标保持暂停。')
    body = f'''# {start}—{end}轮：{title}

2026-10-09整理。{status} 本次仅整理19轮既有内容，不增加轮次、科学组、经验组或认知公理。目录范围闭合不代表完整ROADMAP完成。

[研究总目录](../README.md) · [完整文件索引](文件索引.md) · [原路线图](../ROADMAP.md) · [迁移及只读复算](../_migration/layout_1044_1062_20261009/README.md)

## 正式轮次

{round_table(start, end)}

## 成果与验收入口

'''
    if early:
        body += '''- [阶段综合与输入账](_shared/阶段综合与输入账.md)、[独立完成审阅](_shared/independent_completion_review.md)、[实际结项回执](_shared/goal_completion.json)。
- 1044保未知输入、实际指针、反冲及双方后态；有限共同合同仍有可辨耦合自由，未反证全部认知原则。
- 1045分类其声明后续任务所需资料；静态来源接续不等于动态几何或所有未来操作均闭合。
- 两轮属于不同父对象，不能合称一个完整宇宙模型。[后继进展](../archive_1046_1062/README.md)另行归档。
'''
    else:
        body += '''- [最新整体综合](_shared/overall_evidence_after1062.md)：14类输入、16组现象的已证、采用、限定自由及未决项。
- [综合核验](_shared/overall_evidence_after1062_checks.json)：本批综合计0；不是整条路线图完成证书。
- [1062主线验收](1062/mainline_acceptance.json)：最新正式轮次与累计依据。
- [1044—1045前阶段](../archive_1044_1045/README.md)：已有有限目标结项，不与当前整体目标混同。

## 建议阅读次序

1. **1046—1051：实际载体与有效接口。** 从光子接收、束缚源、材料匹配到自然后态及有限动作商。
2. **1052—1056：共同读数、来源与运动。** 公开弱流、电子来源、光学位置与真实反冲相互接续。
3. **1057—1062：组合条件、经验检验与剩余自由。** 因果组合、α／CKM共同标定、超荷与绝对质量的限定结论。

这是阅读分组，不是三份互相独立的物理实现。1050→1053、1046→1054→1055、1057→1058须沿原依赖复用。

## 跨轮补充与待审内容

成熟文献采用、准入审计、范围修正及整体综合均不另编轮次。跨轮材料仅维护一份，分别保存在`_shared/`与`_admission/`，详见[文件索引](文件索引.md)。

**氢谱共同匹配尚未终签：** [作者稿](_shared/hydrogen_spectral_matching_adoption.md)、[暂停审阅检查点](_shared/hydrogen_spectral_matching_adoption_review.md)。作者轻量复算已完成，缪氢提取来源等审查仍待完成；不得列为已验收采用或1063轮。恢复研究时先处理这一待审项。

已接受的局部共同描述、条件性425接口及有限资源闭合均保留；三维认知生成、全部物理输入的选择与完整统一目标仍未完成。1059的限定不相容、1061的外包未排除、完整联合拟合通过是不同结论。
'''
    return body + '''
## 文件与复算约定

正式报告置于本目录；代码、结果、证明及审阅置于同号目录。历史报告中的“待独审”保留其写作时语义，当前是否验收以主线验收及本索引为准。代码、JSON结果和冻结收据未改字节；Markdown仅更正路径，导航另作可逆记录。

冻结入口使用项目根目录下的`python -B -X utf8 scripts/run_research_recent.py --script <当前脚本路径>`；具体示例见[迁移说明](../_migration/layout_1044_1062_20261009/README.md)。旧逻辑路径只在内存中解析，不恢复旧目录，不依赖ZIP。
'''


def title(path):
    text = path.read_text(encoding='utf-8-sig')
    for line in text.splitlines():
        if line.startswith('# '):
            return line[2:].replace('|', '／').replace('[', '（').replace(']', '）')
    return path.stem


def file_index(old):
    phase = BASE / old
    early = old == 'archive_1044_'
    start, end = (1044, 1045) if early else (1046, 1062)
    out = [f'# {start}—{end}轮文件索引', '', '[阶段说明](README.md) · [研究总目录](../README.md)', '',
           '本索引登记现存文件；编号报告与补充审计分开。冻结收据中的旧路径保留当时语义，复算由迁移入口定位当前文件。', '',
           '## 正式报告及材料', '', round_table(start, end), '', '### 各轮复算与审阅', '']
    for n in range(start, end + 1):
        files = sorted(f for f in (phase/str(n)).iterdir() if f.is_file())
        chosen = [f for f in files if f.suffix == '.py' or 'review' in f.stem or 'acceptance' in f.stem or f.name == 'NEXT.md']
        out.append(f'- **{n}：** ' + ' · '.join(f'[{f.name}]({n}/{f.name})' for f in chosen))
    for folder, label in [('_shared', '跨轮综合、成熟采用及维护'), ('_admission', '准入、去重与未入轮审计')]:
        out += ['', f'## {label}', '', '以下条目不额外计为正式研究轮次。', '']
        if folder == '_shared' and not early:
            out += ['**待终审：** [氢谱作者稿](_shared/hydrogen_spectral_matching_adoption.md)及[暂停检查点](_shared/hydrogen_spectral_matching_adoption_review.md)尚未验收。其余条目仍以各自明确签收范围为准。', '']
        docs = sorted((phase/folder).rglob('*.md'))
        for p in docs:
            relative = p.relative_to(phase).as_posix()
            flag = '【待审】' if 'hydrogen_spectral_matching' in relative else ''
            out.append(f'- {flag}[{title(p)}]({relative})')
    extras = [p for p in phase.glob('*.md') if p.name not in ('README.md','文件索引.md') and not p.name.startswith('research_note_')]
    if extras:
        out += ['', '## 阶段准备', '']
        out += [f'- [{title(p)}]({p.name})' for p in extras]
    out += ['', '## 完整文件与版本记录', '',
            '[迁移清单](../_migration/layout_1044_1062_20261009/manifest.json)逐文件登记位置、大小及前后哈希，覆盖全部代码、JSON、证明与补充材料。正文及来源页链接到各自代码／结果；无需复制一份材料到每个引用它的轮次。', '']
    return '\n'.join(out)


def root_notice():
    return '''## 1044—1062轮归档（2026-10-09）

|阶段目录|范围与状态|
|---|---|
|[archive_1044_1045](archive_1044_1045/README.md)|1044—1045：共同物理过程的耦合自由与记录资料；原有限目标已结项|
|[archive_1046_1062](archive_1046_1062/README.md)|1046—1062：共同接口、公开检验及整体综合；进展归档，整体目标未结项|

最新正式1062，累计科学校准3838；本次整理新增0轮。正式报告在阶段根目录，材料在各轮编号目录，跨轮材料仍仅维护一份。

[整体证据综合](archive_1046_1062/_shared/overall_evidence_after1062.md)已核14类输入与16组现象；[氢谱审阅检查点](archive_1046_1062/_shared/hydrogen_spectral_matching_adoption_review.md)仍未终签。研究目标保持暂停，未自动启动1063，也未修改原ROADMAP目标。恢复前先核待审内容。

[迁移、引用与历史复算说明](_migration/layout_1044_1062_20261009/README.md)。历史代码／结果保持原字节，旧路径经只读解析器定位；历史导航改动只记一份可逆补丁，不保存重复目录或ZIP。下文历史状态不覆盖本栏。

'''


def proposed(path, raw):
    if path.suffix != '.md':
        return raw, []
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    text = raw.decode(enc)
    new = text
    for old in PAIRS:
        if path == BASE/old/'README.md':
            new = stage_readme(old)
        elif path == BASE/old/'文件索引.md':
            new = file_index(old)
    if path.parent == BASE and path.name in MARKERS:
        cut = text.index(MARKERS[path.name])
        new = text.splitlines()[0] + '\n\n' + root_notice() + text[cut:]
    elif path == BASE/'ROADMAP.md':
        cut = text.index('\n') + 1
        new = text[:cut] + '\n' + root_notice() + text[cut:]
    elif path == BASE/'_migration/README.md':
        cut = text.index('\n') + 1
        new = text[:cut] + '\n\n## 2026-10-09：1044—1062归档\n\n[19轮归档、路径核验与只读复算](layout_1044_1062_20261009/README.md)。新入口`scripts/run_research_recent.py`承接1009关闭层及本次两阶段路径；不恢复旧目录、不复制科学材料。\n' + text[cut:]
    new = paths_rewritten(new)
    return new.encode(enc), text_edits(text, new)


def plan():
    assert not MIG.exists(), 'This migration has already been started.'
    for old, new in PAIRS.items():
        src, dst = BASE/old, BASE/new
        assert src.resolve().parent == dst.resolve().parent == BASE.resolve()
        assert src.resolve().is_relative_to(ROOT.resolve()) and dst.resolve().is_relative_to(ROOT.resolve())
        assert src.is_dir() and not src.is_symlink() and not dst.exists()
    moved = {p for old in PAIRS for p in (BASE/old).rglob('*') if p.is_file()}
    candidates = moved | set(BASE.rglob('*.md')) | set((ROOT/'猜想').rglob('*.md'))
    entries = []
    for p in sorted(candidates):
        assert not p.is_symlink()
        raw = p.read_bytes()
        new, edits = proposed(p, raw)
        if p not in moved and raw == new:
            continue
        assert p.suffix not in ('.zip','.pyc','.pyo'), p
        if p.suffix in ('.py','.json'):
            assert raw == new
        e = dict(original=rel(p), destination=rel(mapped(p)), original_sha256=sha(raw),
                 current_sha256=sha(new), original_bytes=len(raw), current_bytes=len(new), edits=edits)
        assert relocation_original_bytes(new, e) == raw
        entries.append(e)
    return dict(schema='research_round_archive_layout_v1', date='2026-10-09',
                pairs={rel(BASE/a): rel(BASE/b) for a,b in PAIRS.items()}, entries=entries,
                moved_file_count=len(moved), latest_formal_round=1062, formal_rounds=19,
                cumulative_scientific_calibrations=3838, new_scientific_groups=0,
                research_goal_status='paused', whole_roadmap_completed=False,
                hydrogen_adoption_accepted=False,
                unchanged_runtime_sha256={p:sha((ROOT/p).read_bytes()) for p in RUNTIME},
                project_readme_sha256=sha((ROOT/'README.md').read_bytes()),
                copied_trees=False, old_directory_shells=False, zip_dependency=False)


def apply():
    data = plan()
    # Validate every original before the first move; refuse concurrent changes.
    for e in data['entries']:
        assert sha((ROOT/e['original']).read_bytes()) == e['original_sha256']
    MIG.mkdir(parents=True)
    with (MIG/'manifest.json').open('x', encoding='utf8', newline='\n') as out:
        json.dump(data, out, ensure_ascii=False, indent=2)
        out.write('\n')
    # Same-parent absolute targets are checked above; move actual directories.
    for old,new in PAIRS.items():
        (BASE/old).rename(BASE/new)
    finish_edits(data)


def finish_edits(data):
    # Resume only an already recorded, fully moved inventory; no new planning.
    for old,new in PAIRS.items():
        assert not (BASE/old).exists() and (BASE/new).is_dir()
    for e in data['entries']:
        assert sha((ROOT/e['destination']).read_bytes()) in (e['original_sha256'],e['current_sha256'])
    for e in data['entries']:
        p = ROOT/e['destination']
        raw = p.read_bytes()
        if sha(raw) == e['current_sha256']:
            continue
        assert sha(raw) == e['original_sha256']
        if not e['edits']:
            continue
        enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
        text = raw.decode(enc)
        for edit in reversed(e['edits']):
            assert text[edit['start']:edit['end']] == edit['before']
            text = text[:edit['start']] + edit['after'] + text[edit['end']:]
        new = text.encode(enc)
        assert sha(new) == e['current_sha256']
        p.write_bytes(new)
    print(json.dumps(dict(moved_files=data['moved_file_count'], manifest_entries=len(data['entries']),
                          changed_markdown=sum(bool(e['edits']) for e in data['entries']), new_rounds=0)))


def verify():
    data = json.loads((MIG/'manifest.json').read_text(encoding='utf8'))
    counts = Counter()
    for e in data['entries']:
        path = ROOT/e['destination']
        raw = path.read_bytes()
        original = relocation_original_bytes(raw,e)
        counts['files'] += 1
        counts['unchanged'] += raw == original
        counts['changed_markdown'] += raw != original
        if path.suffix in ('.py','.json'):
            assert raw == original
            counts['science_code_results'] += 1
    for old,new in PAIRS.items():
        assert not (BASE/old).exists()
        actual = {rel(p) for p in (BASE/new).rglob('*') if p.is_file()}
        expected = {e['destination'] for e in data['entries'] if e['destination'].startswith(rel(BASE/new)+'/')}
        assert actual == expected, (actual-expected,expected-actual)
    for p,h in data['unchanged_runtime_sha256'].items():
        assert sha((ROOT/p).read_bytes()) == h
    assert sha((ROOT/'README.md').read_bytes()) == data['project_readme_sha256']
    missing=[]
    checked=0
    for e in data['entries']:
        path=ROOT/e['destination']
        if path.suffix != '.md':
            continue
        for _,_,target,local in links(path.read_text(encoding='utf-8-sig')):
            if any(path.as_posix().endswith('/'+source) and local==token
                   for source,token in MATH_NOT_LINKS):
                continue
            checked += 1
            if not (path.parent/local.replace('\\','/')).resolve().exists():
                missing.append(dict(source=rel(path),target=target))
    assert not missing,missing
    result=dict(all_layout_checks_passed=True,counts=dict(counts),local_links_checked=checked,
                missing_links=missing,manifest_sha256=sha((MIG/'manifest.json').read_bytes()),
                latest_formal_round=1062,new_rounds=0,goal_status='paused',whole_roadmap_completed=False)
    print(json.dumps(result,ensure_ascii=False))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--plan',action='store_true')
    group.add_argument('--apply',action='store_true')
    group.add_argument('--verify',action='store_true')
    group.add_argument('--finish-recorded-edits',action='store_true')
    args=parser.parse_args()
    if args.apply:
        apply()
    elif args.verify:
        verify()
    elif args.finish_recorded_edits:
        data=json.loads((MIG/'manifest.json').read_text(encoding='utf8'))
        assert data['schema']=='research_round_archive_layout_v1'
        finish_edits(data)
    else:
        data=plan()
        print(json.dumps({k:v for k,v in data.items() if k!='entries'} | {'manifest_entries':len(data['entries']),
                          'changed_markdown':sum(bool(e['edits']) for e in data['entries'])},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
