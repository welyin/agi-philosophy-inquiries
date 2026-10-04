"""Apply the user's final layout: phase-root notes and per-round evidence."""
import json
import os
from pathlib import Path
import re
import zipfile
import organize_research_231_775 as m

ROOT,RESEARCH,OLD=m.ROOT,m.RESEARCH,m.OLD
LAST=RESEARCH/'archive_764_775'
MIG=LAST/'_migration'


def atomic(p,raw):
    temp=p.with_name(p.name+'.round_layout_tmp')
    assert not temp.exists()
    p.parent.mkdir(parents=True,exist_ok=True)
    with temp.open('xb') as f:
        f.write(raw)
    os.replace(temp,p)


def refine():
    manifest=json.loads((MIG/'manifest.json').read_text('utf8'))
    assert manifest['version']==2
    entries=manifest['entries']
    destinations={}
    for e in entries:
        n=e['round_owner']
        match=re.fullmatch(r'research_note_(\d+)\.md',Path(e['original']).name)
        if match and '/' not in e['original']:
            n=int(match[1])
            phase=m.phase_of(n)
            dest=Path(f"archive_{phase['start']}_{phase['end']}")/Path(e['original']).name
        elif n is not None and 231<=n<=775:
            phase=m.phase_of(n)
            base=Path(f"archive_{phase['start']}_{phase['end']}")/str(n)
            parts=Path(e['original']).parts
            if parts[0]==f'round{n}_drafts':
                dest=base/'drafts'/Path(*parts[1:])
            else:
                dest=base/Path(e['original'])
        else:
            dest=Path(e['destination'])
        destinations[e['original']]=dest.as_posix()
    assert len(set(destinations.values()))==8920
    mapping={(OLD/e['original']).resolve():(RESEARCH/destinations[e['original']]).resolve() for e in entries}
    moved_mapping={(RESEARCH/e['destination']).resolve():mapping[(OLD/e['original']).resolve()] for e in entries}
    phase_dirs=[RESEARCH/f"archive_{p['start']}_{p['end']}" for p in manifest['phases']]
    # Capture the current navigators before any moves, to rewrite their links once.
    indexed={str((RESEARCH/e['destination']).resolve()) for e in entries}
    docs=[p for d in phase_dirs for p in d.rglob('*.md') if str(p.resolve()) not in indexed]
    docs += [ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',RESEARCH/'archive_223_230/README.md']
    before={p:p.read_bytes() for p in docs}
    m.NEW=RESEARCH
    new_records=[]
    with zipfile.ZipFile(MIG/'original_workspace_231_775.zip') as z:
        for i,e in enumerate(entries):
            src=RESEARCH/e['destination']
            dst=RESEARCH/destinations[e['original']]
            assert m.digest(src)==e['current_sha256'],e['destination']
            raw=z.read('research_cognition_physics/archive_231_/'+e['original'])
            assert m.sha(raw)==e['original_sha256']
            changed,edits=m.rewrite_links(raw,OLD/e['original'],dst,mapping) if e['rewrite_markdown_links'] else (raw,[])
            assert src.resolve().is_relative_to(RESEARCH.resolve()) and dst.resolve().is_relative_to(RESEARCH.resolve())
            if src!=dst:
                assert not dst.exists(),str(dst)
                dst.parent.mkdir(parents=True,exist_ok=True)
                src.replace(dst)
            if dst.read_bytes()!=changed:
                atomic(dst,changed)
            updated=dict(e)
            updated.update(destination=destinations[e['original']],current_sha256=m.sha(changed),current_bytes=len(changed),link_edits=edits)
            new_records.append(updated)
            if (i+1)%2500==0:
                print(f'Round folders: {i+1}/8920',flush=True)
    for p,raw in before.items():
        assert p.read_bytes()==raw
        changed,_=m.rewrite_links(raw,p,p,moved_mapping)
        text=changed.decode('utf-8-sig')
        text=text.replace('每个archive的README说明问题、成果、假设边界、关键报告和全部轮次；原始报告、代码及结果在其子目录中。',
                          '每个archive根目录直接放README和正式research_note_编号.md；该轮代码、结果、核验及草稿统一放在编号目录（如231/）中。')
        text=text.replace('每阶段包含notes、code、results和必要的drafts。',
                          '每阶段根目录直接放正式research_note_编号.md；其它材料按轮次放入231/、232/等编号目录，工作稿在对应编号目录的drafts中。')
        if p.name=='README.md' and p.parent in phase_dirs:
            text=text.replace('## 研究问题', '正式报告直接位于本目录；配套代码、结果和核验在相应轮次目录中。\n\n## 研究问题',1)
        atomic(p,text.encode('utf8'))
    # Rebuild file indexes so both their links and displayed ownership match the final layout.
    for phase in manifest['phases']:
        base=RESEARCH/f"archive_{phase['start']}_{phase['end']}"
        text='# 文件索引\n\n[阶段README](README.md)。正式报告位于本目录，其余材料按轮次归档；历史脚本通过[统一复算入口](../archive_764_775/_migration/README.md)执行。\n\n|轮次|类型|文件|原路径|\n|---|---|---|---|\n'
        for e in new_records:
            if (RESEARCH/e['destination']).is_relative_to(base):
                text+=f"|{e['round_owner'] if e['round_owner'] else '公共/历史'}|{e['kind']}|[{Path(e['destination']).name}]({m.relative(RESEARCH/e['destination'],base)})|`{e['original']}`|\n"
        atomic(base/'文件索引.md',text.encode('utf8'))
    # Remove only now-empty category directories; no recursive file deletion.
    for base in phase_dirs:
        for d in sorted((p for p in base.rglob('*') if p.is_dir()),key=lambda p:len(p.parts),reverse=True):
            assert d.resolve().is_relative_to(base.resolve())
            if not any(d.iterdir()):
                d.rmdir()
    manifest['version']=3
    manifest['entries']=new_records
    manifest['phase_root_notes_and_round_folders']=True
    old_manifest=MIG/'manifest_sibling_category_layout.json'
    assert not old_manifest.exists()
    (MIG/'manifest.json').replace(old_manifest)
    m.dump(MIG/'manifest.json',manifest)
    m.dump(MIG/'round_layout_navigation_checks.json',dict(
        updated_navigation_hashes={p.relative_to(ROOT).as_posix():m.digest(p) for p in docs},
        root_reports=545,per_round_directories=545,old_category_directories_removed=True,
        new_scientific_rounds=0,goal_still_paused=True))
    print('Final user layout completed: root notes, numbered evidence folders.',flush=True)


if __name__=='__main__':
    refine()
