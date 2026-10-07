"""Physically split completed rounds; retain byte-identical science and reversible references.

Run --plan first, then --apply. No ZIP, mirror tree, or research round is created.
"""
from pathlib import Path
import argparse
from collections import Counter
import hashlib
import json
import os
import re
from urllib.parse import unquote

from organize_research_231_775 import links

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'research_cognition_physics'
OLD = BASE / 'archive_764_'
MIG = BASE / '_migration/layout_764_1008_20261008'
SPEC = ROOT / 'scripts/research_phases_764_1008.json'
PLAN = ROOT / '.research_runtime/split_764_1008_plan.json'
TOP_MOVES = {
    'README.md': '_shared/notes/research_progress_764_1008.md',
    '文件索引.md': '_history/indexes/archive_764_file_index_through1008.md',
    '阶段成果总览.md': '_shared/notes/阶段成果总览.md',
    '跨阶段主题索引.md': '_shared/notes/跨阶段主题索引.md',
}
EXCLUDE = {'.git', '.codex', '.agents', '.research_runtime', '.tmp', '__pycache__', 'node_modules'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def rel(path, base=ROOT):
    return os.path.relpath(path, base).replace('\\', '/')


def safe(path, bound=BASE):
    value = path.resolve()
    assert value.is_relative_to(bound.resolve()) and value != bound.resolve(), value
    return value


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    assert not path.exists(),path
    atomic(path,(json.dumps(value,ensure_ascii=False,indent=2)+'\n').encode('utf8'))


def atomic(path, raw):
    safe(path,ROOT)
    temporary=path.with_name(path.name+'.migration_tmp')
    if temporary.exists():
        assert temporary.read_bytes()==raw,('Unrecognized temporary file',temporary)
    else:
        with temporary.open('xb') as out:
            out.write(raw);out.flush();os.fsync(out.fileno())
    assert sha(temporary.read_bytes())==sha(raw)
    os.replace(temporary,path)


def apply_edits(raw, entry):
    assert sha(raw)==entry['original_sha256'],entry['original']
    if not entry['edits']:return raw
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    text=raw.decode(enc)
    for edit in reversed(entry['edits']):
        assert text[edit['start']:edit['end']]==edit['before']
        text=text[:edit['start']]+edit['after']+text[edit['end']:]
    revised=text.encode(enc)
    assert sha(revised)==entry['current_sha256']
    return revised


def documents():
    for parent, dirs, files in os.walk(ROOT):
        dirs[:] = [name for name in dirs if name not in EXCLUDE]
        for name in files:
            if name.lower().endswith('.md'):
                yield Path(parent) / name


def prepare():
    phases = json.loads(SPEC.read_text('utf8'))['phases']
    assert [n for p in phases for n in range(p['start'], p['end']+1)] == list(range(764,1009))
    assert len(phases) == 10 and OLD.is_dir() and not (MIG/'manifest.json').exists()
    phase_of = {n: f"archive_{p['start']}_{p['end']}" for p in phases for n in range(p['start'],p['end']+1)}
    for phase in set(phase_of.values()):
        assert not (BASE/phase).exists(), phase
    mapping = {}
    for source in sorted(OLD.rglob('*')):
        if not source.is_file():
            continue
        parts = source.relative_to(OLD).parts
        if parts[0].isdigit():
            target = BASE/phase_of[int(parts[0])]/Path(*parts)
        elif re.match(r'^research_note_(\d+)', parts[0]):
            n = int(re.match(r'^research_note_(\d+)',parts[0])[1])
            target = BASE/phase_of[n]/Path(*parts)
        elif parts[0] in ('_shared','_migration','_history'):
            target = BASE/Path(*parts)
        else:
            assert len(parts)==1 and parts[0] in TOP_MOVES, parts
            target = BASE/TOP_MOVES[parts[0]]
        safe(source, OLD); safe(target)
        assert not target.exists(), target
        mapping[source.resolve()] = target.resolve()
    assert len({str(p).casefold() for p in mapping.values()}) == len(mapping)
    directory_map = {OLD.resolve(): BASE.resolve()}
    for name in ('_shared','_migration','_history'):
        directory_map[(OLD/name).resolve()] = (BASE/name).resolve()
    for n,phase in phase_of.items():
        directory_map[(OLD/str(n)).resolve()] = (BASE/phase/str(n)).resolve()

    def mapped(path):
        path = path.resolve()
        if path in mapping:
            return mapping[path]
        for olddir in sorted(directory_map, key=lambda x:len(x.parts),reverse=True):
            if path.is_relative_to(olddir):
                # Keep nonexistent future/historical file references identifiable.
                return directory_map[olddir]/path.relative_to(olddir)
        return path

    raws = {p:p.read_bytes() for p in mapping}
    for path in documents():
        raws.setdefault(path.resolve(),path.read_bytes())
    contents, entries, inherited = {}, [], []
    for src, raw in raws.items():
        dst = mapping.get(src,src)
        edits = []
        changed = raw
        if src.suffix.lower()=='.md':
            enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
            before = raw.decode(enc)
            for a,b,target,local in links(before):
                original_target = (src.parent/local.replace('\\','/')).resolve()
                new_target = mapped(original_target)
                suffix = '#'+target.split('#',1)[1] if '#' in target else ''
                after = rel(new_target,dst.parent)+suffix
                if not original_target.exists():
                    inherited.append([rel(dst),after])
                if after != target:
                    edits.append(dict(start=a,end=b,before=target,after=after,kind='markdown_target'))
            # Concrete project paths in examples keep pointing to current assets.
            occupied = [(e['start'],e['end']) for e in edits]
            pattern = re.compile(r'archive_764_([\\/])(?:(research_note_)(\d+)|(\d+)(?=[\\/])|(_shared|_migration|_history)(?=[\\/]))')
            for match in pattern.finditer(before):
                if any(a<=match.start()<b or a<match.end()<=b for a,b in occupied):
                    continue
                slash=match[1]
                if match[5]:
                    after=match[5]
                else:
                    n=int(match[3] or match[4])
                    if n not in phase_of:
                        continue
                    after=phase_of[n]+slash+(('research_note_'+match[3]) if match[2] else match[4])
                edits.append(dict(start=match.start(),end=match.end(),before=match[0],after=after,kind='literal_path'))
            edits.sort(key=lambda e:e['start'])
            assert all(a['end']<=b['start'] for a,b in zip(edits,edits[1:])),src
            revised = before
            for e in reversed(edits):
                assert revised[e['start']:e['end']]==e['before']
                revised=revised[:e['start']]+e['after']+revised[e['end']:]
            changed=revised.encode(enc)
        if src!=dst or raw!=changed:
            entries.append(dict(original=rel(src),destination=rel(dst),original_sha256=sha(raw),
                current_sha256=sha(changed),original_bytes=len(raw),current_bytes=len(changed),edits=edits))
            contents[src]=(dst,raw,changed)
    return phases, mapping, contents, entries, inherited


def add_edit(entry, original, before_raw, after_raw):
    """Navigation-only changes recorded as exact reversible text patches."""
    import difflib
    enc='utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
    before,after=original.decode(enc),after_raw.decode(enc)
    # Line anchors prevent quadratic character matching in long repetitive navigation.
    # Exact replacement spans still refer to character offsets in the original file.
    left,right=before.splitlines(keepends=True),after.splitlines(keepends=True)
    lo=[0];ro=[0]
    for line in left:lo.append(lo[-1]+len(line))
    for line in right:ro.append(ro[-1]+len(line))
    seq=difflib.SequenceMatcher(None,left,right,autojunk=False)
    entry['edits']=[dict(start=lo[a],end=lo[b],before=before[lo[a]:lo[b]],after=after[ro[c]:ro[d]],kind='navigation_or_reference')
        for tag,a,b,c,d in seq.get_opcodes() if tag!='equal']
    entry['current_sha256']=sha(after_raw); entry['current_bytes']=len(after_raw)


def update_navigation(phases, contents, entries):
    by_original={e['original']:e for e in entries}
    def change(path, fn):
        assert path in contents, path
        dst,original,current=contents[path]
        enc='utf-8-sig' if current.startswith(b'\xef\xbb\xbf') else 'utf8'
        newline='\r\n' if b'\r\n' in current else '\n'
        revised=fn(current.decode(enc).replace('\r\n','\n')).replace('\n',newline).encode(enc)
        add_edit(by_original[rel(path)],original,current,revised)
        contents[path]=(dst,original,revised)
    root_block=('## 764—1008轮已按主题归档（2026-10-08）\n\n'
        '原`archive_764_`已实际拆分为10个连续阶段。正式报告仍放各阶段根目录，代码、结果和草稿放同号轮次目录；'
        '跨阶段材料集中在`_shared`，迁移和复算说明集中在`_migration`。'
        '[本次分期与核验](_migration/layout_764_1008_20261008/README.md)。'
        '科学编号仍止于1008、累计3786；没有新增研究轮次或启动新的目标。\n\n')
    def research_readme(s):
        first,rest=s.split('\n',1)
        s=first+'\n\n'+root_block+rest.lstrip('\n')
        s=s.replace('231轮以后的15个阶段保留既有划分，末阶段 `archive_764_` 为持续研究目录。',
            '231—763原有14个阶段保持；764—1008新增10个主题分期，231—1008合计24个阶段。')
        row=re.search(r'^\|15\|764起[^\n]+$',s,re.M)
        assert row
        rows='\n'.join(f"|{p['id']}|{p['start']}—{p['end']}|[{p['slug']}](archive_{p['start']}_{p['end']}/README.md)|" for p in phases)
        return s[:row.start()]+rows+s[row.end():]
    change(BASE/'README.md',research_readme)
    for name in ('research_direction.md','RESEARCH_STATE.md'):
        def fn(s):
            first,rest=s.split('\n',1)
            s=first+'\n\n'+root_block+rest.lstrip('\n')
            if name=='RESEARCH_STATE.md':
                s=s.replace('## 目标修订：整体机制与决定性检验（2026-10-07）','## 历史目标同步：整体机制与决定性检验（2026-10-07）',1)
                s=s.replace('research_direction.md#当前目标正文从认知操作压缩物理的独立假设',
                    'research_direction.md#当前目标正文完善整体认知操作解释')
            return s
        change(BASE/name,fn)
    def old_progress(s):
        first,rest=s.split('\n',1)
        return ('# 764—1008轮累计进展记录\n\n原持续目录的进展正文保留于此；原文件已按主题实际分期。'
            '[当前阶段目录](../../README.md) · [迁移说明](../../_migration/layout_764_1008_20261008/README.md)。'
            '下面的“当前”“下一轮”等措辞属于记录当时，不表示重新启动已完成工作。\n\n'+rest.lstrip('\n'))
    change(OLD/'README.md',old_progress)
    def shared(s):
        return ('# 跨阶段公共材料\n\n没有唯一轮次归属的专题、代码和结果只保留一份。'
            '[研究总目录](../README.md) · [阶段成果总览](notes/阶段成果总览.md) · '
            '[跨阶段主题](notes/跨阶段主题索引.md) · [累计进展](notes/research_progress_764_1008.md) · '
            '[迁移与复算](../_migration/README.md)。\n')
    change(OLD/'_shared/README.md',shared)
    def migration(s):
        first,rest=s.split('\n',1)
        return (first+'\n\n## 2026-10-08：764—1008再次分期\n\n'
            '旧迁移元数据已集中到本目录，科学文件由当前物理位置读取。'
            '[本次10阶段划分、完整核验及复算说明](layout_764_1008_20261008/README.md)。'
            '新增一层精确路径与文本变更映射，保留旧清单和旧收据，不创建旧目录或ZIP。'
            '776—1008也通过统一入口运行，例如：\n\n'
            '```powershell\npython -B -X utf8 scripts/run_research_current.py --script archive_990_1008/1008/overall_completion_audit.py\n```\n\n'+rest.lstrip('\n'))
    change(OLD/'_migration/README.md',migration)


def generated_documents(phases, entries, mapping):
    generated={}
    for p in phases:
        folder=BASE/f"archive_{p['start']}_{p['end']}"
        result=f"# {p['start']}—{p['end']}轮：{p['slug']}\n\n"
        result+='[研究总目录](../README.md) · [文件索引](文件索引.md) · [复算说明](../_migration/README.md)\n\n'
        result+='正式报告直接放本目录；其余材料按轮次放在编号目录中。原文件只归属一次，跨阶段复用通过引用连接。\n\n'
        result+=f"## 研究问题\n\n{p['question']}\n\n## 阶段成果\n\n"+'\n'.join('- '+x for x in p['results'])+'\n\n'
        result+=f"## 适用范围与未完成项\n\n{p['limits']}\n\n本目录是研究主题分期，不能据此认定其中全部开放问题已解决。原报告的认知动机、额外输入、解析证明、数值校准和物理解释保持区分。\n\n"
        result+='## 关键报告\n\n'+'、'.join(f'[{n}](research_note_{n}.md)' for n in p['keys'])+'\n\n'
        result+=f"## 与后续研究的连接\n\n{p['transition']}\n\n## 全部轮次\n\n|轮次|报告|\n|---|---|\n"
        for n in range(p['start'],p['end']+1):
            name=f'research_note_{n}.md'
            source=OLD/name
            if not source.exists():source=mapping[source.resolve()]
            title=source.read_text('utf-8-sig').splitlines()[0].lstrip('# ').replace('|','／')
            result+=f'|{n}|[{title}]({name})|\n'
        generated[folder/'README.md']=result.encode('utf8')
        index='# 文件索引\n\n[阶段README](README.md) · [统一复算入口](../_migration/README.md)。科学代码与结果保持原字节；报告引用已迁移。\n\n|轮次|文件|\n|---|---|\n'
        for e in entries:
            target=ROOT/e['destination']
            if target.is_relative_to(folder) and '__pycache__' not in target.parts:
                rr=target.relative_to(folder)
                owner=rr.parts[0] if rr.parts[0].isdigit() else re.search(r'research_note_(\d+)',rr.name)[1]
                index+=f'|{owner}|[{rr.as_posix()}]({rr.as_posix()})|\n'
        generated[folder/'文件索引.md']=index.encode('utf8')
    text='# 764—1008轮分期与迁移\n\n2026-10-08。按研究对象和论证转向分为10个连续阶段，245份正式报告全部保留。'
    text+='目录整理不新增科学轮次、不改目标、不运行图像检查；正式仍1008，累计3786。\n\n'
    text+='|轮次|主题|\n|---|---|\n'
    for p in phases:
        text+=f"|{p['start']}—{p['end']}|[{p['slug']}](../../archive_{p['start']}_{p['end']}/README.md)|\n"
    text+='\n## 文件与证据规则\n\n科学Python、JSON结果及旧迁移清单逐字节保留。Markdown只迁移引用，现行导航另更新归档状态；全部修改都有精确可逆记录。'
    text+='公共材料集中在研究根目录的`_shared`、`_history`、`_migration`，原`archive_764_`不保留壳目录。没有ZIP或整套副本。\n\n'
    text+='[本次映射](manifest.json) · [迁移检查](checks.json) · [代表性复算](replay_checks.json) · [全部复算说明](../README.md)\n\n'
    text+='## 冻结脚本的复算\n\n从项目根目录运行统一入口；它将旧逻辑路径映射到真实当前文件，逐层逆转已登记的引用改动后核旧哈希。'
    text+='运行只读，不依赖旧目录。直接运行冻结脚本仍可能按旧布局找跨轮次材料，应使用：\n\n'
    text+='```powershell\npython -B -X utf8 scripts/run_research_current.py --script archive_990_1008/1008/overall_completion_audit.py\npython -B -X utf8 scripts/run_research_current.py --script archive_764_796/775/joint_wick_source_anomaly.py\n```\n\n'
    text+='旧发布核验中“当时导航必须逐字相同”的要求具有历史范围。复算报告分别记录科学复算、历史证据和现行布局，不以迁移成功代替物理证明，也不声称已重跑1008轮全部实验。\n'
    generated[MIG/'README.md']=text.encode('utf8')
    return generated


def finish(data):
    """Resume from verified original/current bytes, never overwrite unknown changes."""
    entries=data['entries']
    mapping={(ROOT/e['original']).resolve():(ROOT/e['destination']).resolve() for e in entries}
    for entry in entries:
        src=safe(ROOT/entry['original'],ROOT);dst=safe(ROOT/entry['destination'],ROOT)
        if src!=dst:
            safe(src,OLD);safe(dst)
            assert not(src.exists() and dst.exists()),('Both copies exist',src,dst)
            assert src.exists() or dst.exists(),('Missing source and destination',src,dst)
            if src.exists():
                assert sha(src.read_bytes())==entry['original_sha256'],src
                dst.parent.mkdir(parents=True,exist_ok=True);src.rename(dst)
        raw=dst.read_bytes();digest=sha(raw)
        assert digest in (entry['original_sha256'],entry['current_sha256']),('Unrecognized current bytes',dst)
        if digest!=entry['current_sha256']:
            atomic(dst,apply_edits(raw,entry))
        assert sha(dst.read_bytes())==entry['current_sha256'],dst
    for path,raw in generated_documents(data['phases'],entries,mapping).items():
        safe(path);path.parent.mkdir(parents=True,exist_ok=True)
        if path.exists():assert path.read_bytes()==raw,('Generated document changed',path)
        else:atomic(path,raw)
    if OLD.exists():
        for path in sorted((p for p in OLD.rglob('*') if p.is_dir()),key=lambda x:len(x.parts),reverse=True):
            safe(path,OLD);assert not any(path.iterdir()),path;path.rmdir()
        safe(OLD);assert not any(OLD.iterdir());OLD.rmdir()


def main():
    parser=argparse.ArgumentParser()
    choices=parser.add_mutually_exclusive_group(required=True)
    choices.add_argument('--plan',action='store_true');choices.add_argument('--apply',action='store_true');choices.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if args.resume:
        data=json.loads((MIG/'manifest.json').read_text('utf8'));finish(data)
        print('Resumed verified migration without replacing unknown bytes.');return
    phases,mapping,contents,entries,inherited=prepare()
    update_navigation(phases,contents,entries)
    generated=generated_documents(phases,entries,mapping)
    source_hashes={rel(p):sha(raw) for p,(_,raw,_) in contents.items()}
    summary=dict(phases=[f"{p['start']}-{p['end']}" for p in phases],
        moved_files=len(mapping),changed_markdown=sum(e['original_sha256']!=e['current_sha256'] for e in entries),
        mapped_entries=len(entries),formal_reports=245,generated_documents=len(generated),
        source_hashes=source_hashes,spec_sha256=sha(SPEC.read_bytes()))
    if args.plan:
        assert not PLAN.exists(); dump(PLAN,summary)
        print(json.dumps({k:v for k,v in summary.items() if k!='source_hashes'},ensure_ascii=False,indent=2));return
    plan=json.loads(PLAN.read_text('utf8')); assert summary==plan,'Inputs changed since plan'
    for path,(dst,raw,new) in contents.items():
        assert path.read_bytes()==raw,('Concurrent change',path)
        safe(dst,ROOT)
        if path!=dst:
            safe(path,OLD);safe(dst);assert not dst.exists()
    assert all(not p.exists() for p in generated)
    # Write the reversible transaction record before moving; never store full file copies.
    data=dict(version=1,date='2026-10-08',phases=phases,entries=entries,
        scientific_code_and_json_unchanged=True,new_scientific_rounds=0,
        inherited_missing_links=inherited,archive_shell_retained=False,zip_created=False)
    dump(MIG/'manifest.json',data)
    finish(data)
    print(json.dumps({k:v for k,v in summary.items() if k!='source_hashes'},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
