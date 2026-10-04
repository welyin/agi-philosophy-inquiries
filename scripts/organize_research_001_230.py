"""Preserve early evidence, group it by round, and reopen the last phase folder.

Only Markdown link targets may change in research evidence. Original bytes,
path mappings and every navigation revision are saved before mutation.
"""
import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import re
import zipfile
import organize_research_231_775 as m

ROOT, BASE = m.ROOT, m.RESEARCH
FIRST, SECOND = BASE/'archive_001_222', BASE/'archive_223_230'
OLD_LAST, LAST = BASE/'archive_764_775', BASE/'archive_764_'
MAINT_REL = Path('_migration/layout_001_230_20261004')


def atomic(p, raw):
    temporary = p.with_name(p.name+'.layout_tmp')
    assert not temporary.exists()
    p.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open('xb') as stream:
        stream.write(raw)
    os.replace(temporary, p)


def read(p):
    return json.loads(p.read_text('utf-8-sig'))


def plan():
    assert OLD_LAST.is_dir() and not LAST.exists()
    index = read(FIRST/'ROUND_INDEX.json')['rounds']
    ownership = defaultdict(set)
    for row in index:
        for artifact in row['indexed_artifacts']:
            ownership[(FIRST/artifact).resolve()].add(row['round'])
        note = FIRST/row['note']
        for _, _, _, local in m.links(note.read_text('utf-8-sig')):
            ownership[(note.parent/local).resolve()].add(row['round'])
    stems = {223:'reversible_dynamics_bridge', 224:'tensor_process_bridge',
             225:'global_orientation_bridge', 226:'steering_permission_bridge',
             227:'reversible_control_bridge', 228:'instrument_completion_bridge',
             229:'continuous_seed_bridge', 230:'finite_protocol_closure'}
    explicit = {'dynamic_complex_emergence_hypothesis.md':205,
                'resource_sustainability_hypothesis.md':204,
                'reconstruction_backchain_checks.json':219,
                'reconstruction_backchain_ledger.json':219}
    entries = []
    for archive in (FIRST, SECOND):
        for source in sorted(p for p in archive.rglob('*') if p.is_file()):
            relative = source.relative_to(archive)
            owner, reason, rewrite = None, '', source.suffix=='.md'
            note = re.fullmatch(r'research_note_(\d+)\.md', source.name)
            if note and (relative.parent==Path('research_process') or relative.parent==Path('.')):
                owner = int(note[1])
                dest = archive/source.name
                reason = 'formal numbered report'
            elif archive==FIRST and relative.parts[0]=='research_process':
                tail = Path(*relative.parts[1:])
                if tail.parts[0]=='__pycache__':
                    dest = archive/'_history/runtime_cache'/Path(*tail.parts[1:])
                    reason, rewrite = 'preserved generated cache; not scientific evidence', False
                elif tail.parts[0]=='history':
                    owner = 219
                    dest = archive/'219/history'/Path(*tail.parts[1:])
                    reason, rewrite = 'original round 219 version and hash contract', False
                elif source.name in ('README.md','research_direction.md','RESEARCH_STATE.md','research_history_catalog.json'):
                    dest = archive/'_history/indexes'/source.name
                    reason = 'historical stage navigation'
                elif source.name=='requirements_research_optional.txt':
                    dest = archive/'_shared'/source.name
                    reason = 'shared optional runtime requirements'
                else:
                    numbers = re.findall(r'(?:research_round_|research_review_)(\d+)', source.name)
                    result_code = source.with_name(source.name.replace('_results.json','.py'))
                    candidates = sorted(ownership[source.resolve()] or ownership[result_code.resolve()])
                    if numbers:
                        owner, reason = int(numbers[0]), 'explicit round in name; joint checks stored with first round'
                    elif source.name in explicit:
                        owner, reason = explicit[source.name], 'explicit provenance in memorandum or round 219 backchain'
                    elif candidates:
                        owner, reason = min(candidates), 'earliest indexed artifact or direct report reference'
                    else:
                        raise AssertionError('Unattributed research material: '+str(source))
                    dest = archive/str(owner)/tail
            elif archive==FIRST:
                if source.name in ('认知联合体架构-AGI工程方案.md','研究方向讨论：从认知架构到物理世界.md','article8_hypothesis.md'):
                    dest = archive/'_shared'/relative
                    reason = 'unnumbered background reference'
                else:
                    owner = 222
                    dest = archive/'222/archive_closure'/relative
                    reason = 'original first-stage archiving and paper maintenance'
                    if source.name=='README.md' and relative.parent==Path('.'):
                        dest = dest.with_name('README_before_round_layout.md')
                    rewrite = False
            else:
                numbers = re.findall(r'research_round_(\d+)', source.name)
                owner = next((n for n,stem in stems.items() if source.name in (stem+'.py',stem+'_results.json')), None)
                if numbers:
                    owner = int(numbers[0])
                if owner:
                    dest = archive/str(owner)/relative
                    reason = 'numbered round or experiment module'
                elif source.name=='round229_initial_index.md':
                    owner, dest = 229, archive/'229'/relative
                    reason = 'round 229 provisional index'
                else:
                    owner = 228 if source.name in ('STAGE2_MANIFEST.json','stage2_closure_checks.json','stage2_math_review_checks.json','review_stage2_math.mjs','verify_stage2.py') else 230
                    dest = archive/str(owner)/'archive_closure'/relative
                    reason = 'stage closure and paper versions; original hash contracts retained'
                    if source.name=='README.md' and relative.parent==Path('.'):
                        dest = dest.with_name('README_before_round_layout.md')
                    if relative.parts[0]=='paper_versions' or source.name=='README.md':
                        rewrite = False
            raw = source.read_bytes()
            entries.append(dict(original=source.relative_to(ROOT).as_posix(),
                                destination=dest.relative_to(ROOT).as_posix(),
                                round_owner=owner,attribution=reason,rewrite_markdown_links=rewrite,
                                original_sha256=m.sha(raw),bytes=len(raw)))
    assert len({e['destination'] for e in entries})==len(entries)
    assert sorted(e['round_owner'] for e in entries if e['attribution']=='formal numbered report')==list(range(1,231))
    return dict(version=1,rounds=[1,230],entries=entries,round_index=index,
                active_directory='research_cognition_physics/archive_764_',next_round=776,
                science_content_changes=0,new_research_rounds=0)


def add_directory_maps(mapping):
    groups = defaultdict(list)
    for source, dest in list(mapping.items()):
        for directory in source.parents:
            if directory==BASE or not directory.is_relative_to(BASE):
                break
            groups[directory].append(dest.parent)
    for directory, targets in groups.items():
        mapping.setdefault(directory, Path(os.path.commonpath(targets)))
    mapping[FIRST.resolve()] = FIRST.resolve()
    mapping[SECOND.resolve()] = SECOND.resolve()
    mapping[OLD_LAST.resolve()] = LAST.resolve()
    # Commonpath works for the flattened early layouts; preserve subdirectories
    # under the unchanged late archive by their exact prefix mapping.
    for directory in [OLD_LAST]+[p for p in OLD_LAST.rglob('*') if p.is_dir()]:
        mapping[directory.resolve()] = (LAST/directory.relative_to(OLD_LAST)).resolve()


def compose_link_edits(original, current):
    enc = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
    before, after = original.decode(enc), current.decode(enc)
    left, right = list(m.links(before)), list(m.links(after))
    assert len(left)==len(right)
    edits = [dict(start=a,end=b,before=t,after=r[2]) for (a,b,t,_),r in zip(left,right) if t!=r[2]]
    rebuilt = before
    for edit in reversed(edits):
        rebuilt = rebuilt[:edit['start']]+edit['after']+rebuilt[edit['end']:]
    assert rebuilt.encode(enc)==current, 'Non-link evidence modification'
    return edits


def apply():
    proposal = plan()
    entries = proposal['entries']
    old_manifest = read(OLD_LAST/'_migration/manifest.json')
    original_late = {(BASE/e['destination']).resolve():e for e in old_manifest['entries']}
    originals = {(ROOT/e['original']).resolve():e for e in entries}
    mapping = {source:(ROOT/e['destination']).resolve() for source,e in originals.items()}
    # Include every last-phase file, including navigation and replay utilities.
    for source in OLD_LAST.rglob('*'):
        if source.is_file():
            mapping[source.resolve()] = (LAST/source.relative_to(OLD_LAST)).resolve()
    mapping[(OLD_LAST/'next_round_776/STATUS.md').resolve()] = (LAST/'776/drafts/STATUS.md').resolve()
    add_directory_maps(mapping)
    # Incoming navigation should open the replacement live indexes, rather than
    # the preserved pre-migration copies of those indexes.
    mapping[(FIRST/'README.md').resolve()] = FIRST/'README.md'
    mapping[(SECOND/'README.md').resolve()] = SECOND/'README.md'
    mapping[(FIRST/'ROUND_INDEX.md').resolve()] = FIRST/'ROUND_INDEX.md'
    mapping[(OLD_LAST/'next_round_776').resolve()] = (LAST/'776/drafts').resolve()
    # Historical stage-2 relative alias was already documented by its verifier.
    mapping[(BASE/'archive_223_/README.md').resolve()] = SECOND/'README.md'

    docs = {}
    for p in ROOT.rglob('*.md'):
        rel = p.relative_to(ROOT)
        if any(part in ('.git','.agents','.codex','.research_runtime','node_modules') for part in rel.parts):
            continue
        src = p.resolve()
        if src in originals and not originals[src]['rewrite_markdown_links']:
            continue
        if src in original_late and not original_late[src]['rewrite_markdown_links']:
            continue
        # All raw migration snapshots/backups keep their original link semantics.
        if '_migration' in rel.parts and src not in original_late and p.name!='README.md':
            continue
        raw = p.read_bytes()
        dst = mapping.get(src,src)
        changed, edits = m.rewrite_links(raw,src,dst,mapping)
        if changed!=raw:
            docs[src] = (dst,raw,changed,edits)

    maintenance = OLD_LAST/MAINT_REL
    assert not maintenance.exists()
    maintenance.mkdir(parents=True)
    snapshot_items = set(originals) | set(docs)
    snapshot_items.update([OLD_LAST/'_migration/manifest.json', OLD_LAST/'_migration/replay.py',
                           ROOT/'scripts/verify_research_phase_migration.py', ROOT/'scripts/research_phases_231_775.json'])
    snapshot_items.update([FIRST/'README.md',SECOND/'README.md',OLD_LAST/'README.md',
                           OLD_LAST/'文件索引.md',BASE/'README.md',BASE/'research_direction.md',BASE/'RESEARCH_STATE.md',ROOT/'README.md'])
    checks = []
    with zipfile.ZipFile(maintenance/'before_layout.zip','x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(snapshot_items):
            raw = p.read_bytes()
            key = p.relative_to(ROOT).as_posix()
            z.writestr(key,raw)
            checks.append(dict(path=key,sha256=m.sha(raw),bytes=len(raw)))
    with zipfile.ZipFile(maintenance/'before_layout.zip') as z:
        for item in checks:
            assert m.sha(z.read(item['path']))==item['sha256']
    m.dump(maintenance/'snapshot.json',dict(file='before_layout.zip',sha256=m.digest(maintenance/'before_layout.zip'),entries=checks))
    m.dump(maintenance/'plan.json',proposal)

    # Absolute boundaries and destination collisions are checked before any move.
    for src,e in originals.items():
        dst = ROOT/e['destination']
        assert src.is_relative_to(BASE.resolve()) and dst.resolve().is_relative_to(BASE.resolve())
        assert m.digest(src)==e['original_sha256']
        assert src==dst or not dst.exists(),str(dst)
    assert OLD_LAST.resolve().parent==BASE.resolve() and LAST.resolve().parent==BASE.resolve()
    OLD_LAST.rename(LAST)
    maintenance = LAST/MAINT_REL
    for src,e in originals.items():
        dst = ROOT/e['destination']
        if src!=dst:
            dst.parent.mkdir(parents=True,exist_ok=True)
            src.replace(dst)
    old_status, new_status = LAST/'next_round_776/STATUS.md', LAST/'776/drafts/STATUS.md'
    new_status.parent.mkdir(parents=True,exist_ok=True)
    old_status.replace(new_status)

    changed_documents = []
    for src,(dst,raw,changed,edits) in docs.items():
        assert dst.read_bytes()==raw,str(dst)
        atomic(dst,changed)
        changed_documents.append(dict(original=src.relative_to(ROOT).as_posix(),destination=dst.relative_to(ROOT).as_posix(),
                                      original_sha256=m.sha(raw),current_sha256=m.sha(changed),link_edits=edits))
    for e in entries:
        p = ROOT/e['destination']
        raw = p.read_bytes()
        e.update(current_sha256=m.sha(raw),current_bytes=len(raw),link_edits=docs.get((ROOT/e['original']).resolve(),(None,None,None,[]))[3])
    proposal.pop('round_index')
    proposal.update(original_file_count=len(entries),link_navigation_changes=changed_documents,
                    original_code_results_and_receipts_unchanged=True)
    m.dump(maintenance/'manifest.json',proposal)

    # Preserve the late archive's original scientific manifest contract.
    frozen = LAST/'_migration/original_workspace_231_775.zip'
    with zipfile.ZipFile(frozen) as z:
        for e in old_manifest['entries']:
            previous = (BASE/e['destination']).resolve()
            dst = mapping.get(previous,previous)
            raw = dst.read_bytes()
            original = z.read('research_cognition_physics/archive_231_/'+e['original'])
            e['destination'] = dst.relative_to(BASE).as_posix()
            e.update(current_sha256=m.sha(raw),current_bytes=len(raw),link_edits=compose_link_edits(original,raw) if e['rewrite_markdown_links'] else [])
            if not e['rewrite_markdown_links']:
                assert raw==original,e['destination']
        for repair in old_manifest.get('preexisting_draft_link_repairs',[]):
            for key in ('file','resolved_to'):
                if key in repair:
                    repair[key] = repair[key].replace('archive_764_775/','archive_764_/')
    old_manifest.update(version=4,active_directory='archive_764_',early_layout_manifest='layout_001_230_20261004/manifest.json')
    for phase in old_manifest['phases']:
        phase['directory'] = 'archive_764_' if phase['start']==764 else f"archive_{phase['start']}_{phase['end']}"
    atomic(LAST/'_migration/manifest.json',(json.dumps(old_manifest,ensure_ascii=False,indent=2)+'\n').encode('utf8'))
    for archive in (FIRST,SECOND,LAST):
        for directory in sorted((p for p in archive.rglob('*') if p.is_dir()),key=lambda p:len(p.parts),reverse=True):
            assert directory.resolve().is_relative_to(archive.resolve())
            if not any(directory.iterdir()):
                directory.rmdir()
    for archive,lo,hi in ((FIRST,1,222),(SECOND,223,230)):
        for number in range(lo,hi+1):
            (archive/str(number)).mkdir(exist_ok=True)
    print(json.dumps(dict(early_files=len(entries),updated_markdown=len(docs),active_directory=str(LAST),reports=230),ensure_ascii=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply',action='store_true')
    args=parser.parse_args()
    if args.apply:
        apply()
    else:
        p=plan()
        print(json.dumps(dict(files=len(p['entries']),rounds=p['rounds'],unowned=[e['destination'] for e in p['entries'] if e['round_owner'] is None and '/runtime_cache/' not in e['destination']]),ensure_ascii=False,indent=2))
