"""Split the already organized early rounds at documented research transitions."""
from collections import defaultdict
import json
import os
from pathlib import Path
import re
import zipfile
import organize_research_231_775 as m
from organize_research_001_230 import atomic, compose_link_edits

ROOT,BASE=m.ROOT,m.RESEARCH
OLD=BASE/'archive_001_222'
LAST=BASE/'archive_764_'
MAINT=LAST/'_migration/layout_001_230_20261004'
PHASES=json.loads((ROOT/'scripts/research_phases_001_222.json').read_text('utf8'))['phases']
assert [n for p in PHASES for n in range(p['start'],p['end']+1)]==list(range(1,223))
for p in PHASES:
    p['directory']=f"archive_{p['start']:03d}_{p['end']:03d}"
END=BASE/PHASES[-1]['directory']
OVERVIEW=END/'001—222阶段总览.md'


def phase_of(n):
    return next(p for p in PHASES if p['start']<=n<=p['end'])


def read(p):
    return json.loads(p.read_text('utf8'))


def apply():
    assert OLD.is_dir() and all(not (BASE/p['directory']).exists() for p in PHASES)
    early=read(MAINT/'manifest.json')
    late=read(LAST/'_migration/manifest.json')
    early_by={(ROOT/e['destination']).resolve():e for e in early['entries']}
    late_by={(BASE/e['destination']).resolve():e for e in late['entries']}
    mapping={}
    moved=[]
    frozen=[]
    for source in sorted(p for p in OLD.rglob('*') if p.is_file()):
        e=early_by.get(source.resolve())
        if e:
            rel=source.relative_to(OLD)
            archive=BASE/phase_of(e['round_owner'])['directory'] if e['round_owner'] else END
            dest=archive/rel
            rewrite=e['rewrite_markdown_links']
        else:
            assert source.name in ('README.md','ROUND_INDEX.md','文件索引.md') and source.parent==OLD
            dest=END/'_history/pre_phase_split'/source.name
            rewrite=False
            frozen.append(dest.relative_to(ROOT).as_posix())
        mapping[source.resolve()]=dest.resolve()
        raw=source.read_bytes()
        moved.append(dict(original=source.relative_to(ROOT).as_posix(),destination=dest.relative_to(ROOT).as_posix(),
                          sha256=m.sha(raw),bytes=len(raw),rewrite_links=rewrite))
    directories=defaultdict(list)
    for src,dst in mapping.items():
        for parent in src.parents:
            if not parent.is_relative_to(OLD) or parent==OLD:
                break
            directories[parent].append(dst.parent)
    for directory,targets in directories.items():
        mapping[directory]=Path(os.path.commonpath(targets))
    mapping[OLD.resolve()]=END.resolve()
    mapping[(OLD/'README.md').resolve()]=OVERVIEW
    mapping[(OLD/'ROUND_INDEX.md').resolve()]=OVERVIEW
    mapping[(OLD/'文件索引.md').resolve()]=OVERVIEW
    # A project-level reference was moved independently; verify its actual
    # identity against the old workspace snapshot before redirecting links.
    external=ROOT/'猜想/认知联合体架构-AGI工程方案.md'
    with zipfile.ZipFile(LAST/'_migration/original_workspace_231_775.zip') as z:
        assert z.read('认知联合体架构-AGI工程方案.md')==external.read_bytes()
    mapping[(ROOT/'认知联合体架构-AGI工程方案.md').resolve()]=external.resolve()

    by_move={(ROOT/e['original']).resolve():e for e in moved}
    docs={}
    for source in ROOT.rglob('*.md'):
        rel=source.relative_to(ROOT)
        if any(t in ('.git','.agents','.codex','.research_runtime','node_modules','navigation_before_final_edit') for t in rel.parts):
            continue
        src=source.resolve()
        if src in by_move and not by_move[src]['rewrite_links']:
            continue
        evidence=early_by.get(src) or late_by.get(src)
        if evidence and not evidence['rewrite_markdown_links']:
            continue
        raw=source.read_bytes()
        dst=(ROOT/by_move[src]['destination']).resolve() if src in by_move else src
        changed,edits=m.rewrite_links(raw,src,dst,mapping)
        if raw!=changed:
            docs[src]=(dst,raw,changed,edits)
    snapshot=set(by_move)|set(docs)|{MAINT/'manifest.json',LAST/'_migration/manifest.json',MAINT/'round_index_001_222.json',BASE/'README.md',BASE/'RESEARCH_STATE.md',BASE/'research_direction.md',ROOT/'README.md'}
    receipt=[]
    zip_path=MAINT/'before_early_phase_split.zip'
    with zipfile.ZipFile(zip_path,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for source in sorted(snapshot):
            raw=source.read_bytes()
            name=source.relative_to(ROOT).as_posix()
            z.writestr(name,raw)
            receipt.append(dict(path=name,sha256=m.sha(raw),bytes=len(raw)))
    with zipfile.ZipFile(zip_path) as z:
        for e in receipt:
            assert m.sha(z.read(e['path']))==e['sha256']
    m.dump(MAINT/'phase_split_snapshot.json',dict(file=zip_path.name,sha256=m.digest(zip_path),entries=receipt))
    m.dump(MAINT/'phase_split_plan.json',dict(phases=PHASES,files=moved,raw_historical_navigation=frozen))
    for e in moved:
        src,dst=ROOT/e['original'],ROOT/e['destination']
        assert src.resolve().is_relative_to(OLD.resolve()) and dst.resolve().is_relative_to(BASE.resolve())
        assert not dst.exists() and m.digest(src)==e['sha256']
    for e in moved:
        src,dst=ROOT/e['original'],ROOT/e['destination']
        dst.parent.mkdir(parents=True,exist_ok=True)
        src.replace(dst)
    for src,(dst,raw,changed,edits) in docs.items():
        assert dst.read_bytes()==raw
        atomic(dst,changed)
    with zipfile.ZipFile(MAINT/'before_layout.zip') as z:
        for e in early['entries']:
            previous=(ROOT/e['destination']).resolve()
            dst=ROOT/by_move[previous]['destination'] if previous in by_move else previous
            raw=dst.read_bytes()
            original=z.read(e['original'])
            e['destination']=dst.relative_to(ROOT).as_posix()
            e.update(current_sha256=m.sha(raw),current_bytes=len(raw),link_edits=compose_link_edits(original,raw) if e['rewrite_markdown_links'] else [])
            if not e['rewrite_markdown_links']:
                assert raw==original
    early.update(version=2,early_phases=PHASES,phase_split_scientific_changes=0)
    atomic(MAINT/'manifest.json',(json.dumps(early,ensure_ascii=False,indent=2)+'\n').encode('utf8'))
    with zipfile.ZipFile(LAST/'_migration/original_workspace_231_775.zip') as z:
        for e in late['entries']:
            dst=BASE/e['destination']
            raw=dst.read_bytes()
            original=z.read('research_cognition_physics/archive_231_/'+e['original'])
            e.update(current_sha256=m.sha(raw),current_bytes=len(raw),link_edits=compose_link_edits(original,raw) if e['rewrite_markdown_links'] else [])
            if not e['rewrite_markdown_links']:
                assert raw==original
    late['version']=5
    atomic(LAST/'_migration/manifest.json',(json.dumps(late,ensure_ascii=False,indent=2)+'\n').encode('utf8'))
    for directory in sorted((p for p in OLD.rglob('*') if p.is_dir()),key=lambda p:len(p.parts),reverse=True):
        assert directory.resolve().is_relative_to(OLD.resolve())
        if not any(directory.iterdir()):
            directory.rmdir()
    assert not any(OLD.iterdir()) and OLD.resolve().parent==BASE.resolve()
    OLD.rmdir()
    for phase in PHASES:
        archive=BASE/phase['directory']
        for n in range(phase['start'],phase['end']+1):
            (archive/str(n)).mkdir(exist_ok=True)
    m.dump(MAINT/'phase_split_link_changes.json',dict(changes=[dict(original=s.relative_to(ROOT).as_posix(),destination=d.relative_to(ROOT).as_posix(),before_sha256=m.sha(raw),after_sha256=m.sha(changed),link_edits=edits) for s,(d,raw,changed,edits) in docs.items()],external_reference_identity_verified=True))
    write_navigation(early)
    print('Split early rounds into',len(PHASES),'phases; moved',len(moved),'files; link-only edits',len(docs),flush=True)


def write_navigation(early):
    navigation=[]
    def write(p,text):
        raw=text.encode('utf8')
        if p.exists():
            prior=p.read_bytes()
            backup=MAINT/'navigation_before_phase_summary'/p.relative_to(ROOT)
            assert not backup.exists()
            backup.parent.mkdir(parents=True,exist_ok=True);backup.write_bytes(prior)
        atomic(p,raw)
        navigation.append(dict(path=p.relative_to(ROOT).as_posix(),sha256=m.sha(raw)))
    index=read(MAINT/'round_index_001_222.json')['rounds']
    by_round={row['round']:row for row in index}
    entries=early['entries']
    notes={e['round_owner']:ROOT/e['destination'] for e in entries if e['attribution']=='formal numbered report'}
    for row in index:
        n=row['round']
        row['note']=notes[n].relative_to(BASE).as_posix()
        row['current_note_sha256']=m.digest(notes[n])
        row['indexed_artifacts']=[(BASE/phase_of(n)['directory']/a).relative_to(BASE).as_posix() for a in row['indexed_artifacts']]
    atomic(MAINT/'round_index_001_222.json',(json.dumps(dict(count=222,rounds=index,path_base='research_cognition_physics'),ensure_ascii=False,indent=2)+'\n').encode('utf8'))
    summary='# 001—222轮阶段总览\n\n[研究总目录](../README.md) · [阶段论文](../可组合认知结构与复量子状态空间_阶段论文.md)。按问题转折划分11段；这是阅读组织，不改变原结论、轮次与历史适用范围。\n\n|轮次|阶段|转入下一阶段的原因|\n|---|---|---|\n'
    for i,phase in enumerate(PHASES):
        archive=BASE/phase['directory']
        summary+=f"|{phase['start']:03d}—{phase['end']:03d}|[{phase['title']}](../{phase['directory']}/README.md)|{phase['transition']}|\n"
        text=f"# 第{phase['start']:03d}—{phase['end']:03d}轮：{phase['title']}\n\n[研究总目录](../README.md) · [001—222总览](../{END.name}/{OVERVIEW.name}) · [文件索引](文件索引.md)\n\n正式编号报告直接放本目录；代码、结果、检查和补充稿放相应轮次的数字目录。\n\n## 核心问题\n\n{phase['question']}\n\n## 阶段成果\n\n"
        text+=''.join(f'- {result}\n' for result in phase['results'])
        text+=f"\n## 适用范围\n\n{phase['limits']}\n\n## 关键报告\n\n"+'、'.join(f'[{n}]({notes[n].name})' for n in phase['keys'])+'\n\n## 前后衔接\n\n'
        if i:
            p=PHASES[i-1];text+=f"前一阶段：[第{p['start']:03d}—{p['end']:03d}轮](../{p['directory']}/README.md)。\n\n"
        text+=phase['transition']+'\n\n'
        next_dir=PHASES[i+1]['directory'] if i+1<len(PHASES) else 'archive_223_230'
        text+=f'[后一阶段](../{next_dir}/README.md)。\n\n'
        if phase['start']==217:
            text+='本目录同时保存001—222阶段的公共背景（`_shared/`）、原导航和缓存（`_history/`）及222轮结项材料（`222/archive_closure/`）；这些公共原件不表示均产生于217轮以后。[219修订前原稿](219/history/round219_v1/archive_manifest.json)继续保留。\n\n'
        text+='## 全部轮次\n\n|轮次|报告|原索引结论|范围|\n|---|---|---|---|\n'
        for n in range(phase['start'],phase['end']+1):
            row=by_round[n]
            clean=lambda s:s.replace('|','\\|').replace('\n',' ')
            text+=f"|{n}|[{clean(row['title'])}]({notes[n].name})|{clean(row['indexed_conclusion'])}|{clean(row['indexed_scope'])}|\n"
        text+='\n## 复算和原件\n\n代码、JSON结果与历史收据保持原字节，阅读版报告仅调整链接。通过[统一早期复算入口](../archive_764_/_migration/layout_001_230_20261004/README.md)恢复原布局，避免历史跨轮导入和固定路径失效；[本次原测试记录](../archive_764_/_migration/layout_001_230_20261004/replay_checks.json)与迁移核验分开保存。\n'
        write(archive/'README.md',text)
        text='# 配套文件索引\n\n[阶段README](README.md)。以下是实际迁入本阶段的原文件；共用文件按最早已确认轮次归属，不复制科研证据。缓存保留但不逐项列入。\n\n|轮次|文件|迁移前路径|\n|---|---|---|\n'
        for e in entries:
            p=ROOT/e['destination']
            if not p.is_relative_to(archive) or '/runtime_cache/' in e['destination']:
                continue
            text+=f"|{e['round_owner'] or '公共/历史'}|[{p.name}]({m.relative(p,archive)})|`{e['original']}`|\n"
        write(archive/'文件索引.md',text)
    summary+='\n## 保留的结项范围\n\n前期定理和反例的输入合同保持原样。指定经典受限模型、参考辅助实模型以及条件复闭包属于不同适用范围；不能由目录阶段标题抹去前提。219修订稿和222校准—层析桥接作为阶段论文继承关系的一部分。\n\n论文：[可组合认知结构与复量子状态空间](../可组合认知结构与复量子状态空间_阶段论文.md)。后续有限维量子过程见[223—230](../archive_223_230/README.md)，当前继续研究见[archive_764_](../archive_764_/README.md)。\n'
    write(OVERVIEW,summary)
    p=BASE/'README.md';text=p.read_text('utf8')
    block='## 001—222阶段目录\n\n|轮次|主题|\n|---|---|\n'+''.join(f"|{p['start']:03d}—{p['end']:03d}|[{p['title']}]({p['directory']}/README.md)|\n" for p in PHASES)
    block+=f'\n[前期阶段总览]({END.name}/{OVERVIEW.name})。223—230仍为[有限维量子理论阶段](archive_223_230/README.md)。\n\n'
    text=text.replace('## 231—775阶段目录',block+'## 231—775阶段目录',1)
    text=text.replace('001—222档案','001—222分期总览')
    write(p,text)
    p=BASE/'RESEARCH_STATE.md';text=p.read_text('utf8')
    text=text.replace('001—230新增完成按轮次整理；','001—222进一步按问题转折分为11个阶段，223—230保留原阶段；')
    text+=f'\n001—222的[分期依据与阅读入口]({END.name}/{OVERVIEW.name})。原 `archive_001_222` 已实际拆分，不保留第二套重复材料。\n'
    write(p,text)
    p=BASE/'research_direction.md';text=p.read_text('utf8')
    text=text.replace('本次完成001—230轮格式统一，并将末阶段更名为archive_764_','本次完成001—230轮格式统一，将001—222按问题转折分成11段，并将末阶段更名为archive_764_')
    write(p,text)
    p=MAINT/'README.md';text=p.read_text('utf8')
    text=text.replace('001—222和223—230两个阶段保留原边界，230篇正式报告放各阶段根目录，配套材料归入对应数字目录。',f'001—222进一步拆分为11个[主题阶段](../../../{END.name}/{OVERVIEW.name})；223—230仍保留原边界。230篇正式报告放相应阶段根目录，配套材料归入数字目录。')
    text+='\n## 追加分期的追溯\n\n[分期前快照](before_early_phase_split.zip)、[分期快照哈希](phase_split_snapshot.json)、[逐文件移动计划](phase_split_plan.json)和[链接变更](phase_split_link_changes.json)记录本轮用户追加要求。旧001—222整包导航仅作为历史版本保留。\n'
    write(p,text)
    m.dump(MAINT/'phase_split_navigation.json',dict(documents=navigation,phases=PHASES,new_scientific_rounds=0))


if __name__=='__main__':
    apply()
