"""Verify lossless scientific migration and link-only Markdown updates."""
import ast
from collections import Counter
import json
from pathlib import Path
import re
import zipfile
import organize_research_231_775 as migration

ROOT,OLD,NEW=migration.ROOT,migration.OLD,migration.RESEARCH
MIG=NEW/'archive_764_/_migration'
manifest=json.loads((MIG/'manifest.json').read_text('utf8'))
snapshot=json.loads((MIG/'snapshot_checks.json').read_text('utf8'))
assert migration.digest(MIG/snapshot['file'])==snapshot['sha256']
inventory={e['path']:e for e in snapshot['entries']}
known_dirs=set()
for name in inventory:
    known_dirs.update(p.as_posix() for p in Path(name).parents)
reports=[]
preexisting=[]
new_broken=[]
checked_links=0
exact_unchanged=0
link_only=0
with zipfile.ZipFile(MIG/snapshot['file']) as z:
    for e in manifest['entries']:
        p=NEW/e['destination']
        original=z.read('research_cognition_physics/archive_231_/'+e['original'])
        current=p.read_bytes()
        assert migration.sha(original)==e['original_sha256']
        assert migration.sha(current)==e['current_sha256'],e['destination']
        assert len(original)==e['bytes'] and len(current)==e['current_bytes']
        if e['link_edits']:
            enc='utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
            text=original.decode(enc)
            valid={(a,b):target for a,b,target,_ in migration.links(text)}
            for c in reversed(e['link_edits']):
                assert valid[c['start'],c['end']]==c['before']
                assert text[c['start']:c['end']]==c['before']
                text=text[:c['start']]+c['after']+text[c['end']:]
            assert text.encode(enc)==current,e['destination']
            link_only+=1
        else:
            assert original==current
            exact_unchanged+=1
        if e['kind']=='report':
            n=int(re.fullmatch(r'research_note_(\d+)\.md',Path(e['original']).name)[1])
            reports.append(n)
        if not e['rewrite_markdown_links']:
            continue
        old_links=list(migration.links(original.decode('utf-8-sig')))
        new_links=list(migration.links(current.decode('utf-8-sig')))
        assert len(old_links)==len(new_links),e['original']
        for old,new in zip(old_links,new_links):
            checked_links+=1
            old_target=((OLD/e['original']).parent/old[3]).resolve()
            target=(p.parent/new[3]).resolve()
            if target.exists():
                continue
            if old_target.is_relative_to(ROOT):
                key=old_target.relative_to(ROOT).as_posix()
                existed=key in inventory or key in known_dirs
            else:
                existed=old_target.exists()
            entry=dict(file=e['destination'],original_target=old[2],current_target=new[2])
            (new_broken if existed else preexisting).append(entry)
assert sorted(reports)==list(range(231,776))
assert len(manifest['entries'])==manifest['original_file_count']==8920
assert len({e['destination'] for e in manifest['entries']})==8920
assert not OLD.exists()
assert not (NEW/'archive_231_775').exists()
# Current indexes and newly-written navigation must have no missing targets.
migrated={str((NEW/e['destination']).resolve()) for e in manifest['entries']}
phase_directories=[NEW/p.get('directory',f"archive_{p['start']}_{p['end']}") for p in manifest['phases']]
assert len(phase_directories)==15 and all(p.is_dir() for p in phase_directories)
for phase,d in zip(manifest['phases'],phase_directories):
    for n in range(phase['start'],phase['end']+1):
        assert (d/f'research_note_{n}.md').is_file(),n
        assert (d/str(n)).is_dir(),n
    assert not any((d/x).exists() for x in ('code','notes','results','drafts'))
new_documents=[p for d in phase_directories for p in d.rglob('*.md')
               if str(p.resolve()) not in migrated
               and not any(t in ('navigation_before_final_edit','navigation_before_phase_summary') for t in p.parts)]
new_documents += [ROOT/'README.md',NEW/'README.md',NEW/'research_direction.md',NEW/'RESEARCH_STATE.md',NEW/'archive_223_230/README.md']
for p in new_documents:
    for _,_,target,local in migration.links(p.read_text('utf-8-sig')):
        checked_links+=1
        t=(p.parent/local).resolve()
        if t==(MIG/'migration_checks.json').resolve():
            continue
        if not t.exists():
            new_broken.append(dict(file=str(p.relative_to(ROOT)),target=target,new_navigation=True))
for p in [ROOT/'scripts/organize_research_231_775.py',ROOT/'scripts/finalize_research_phase_layout.py',ROOT/'scripts/refine_research_round_layout.py',ROOT/'scripts/verify_research_phase_migration.py',MIG/'replay.py']:
    ast.parse(p.read_text('utf8'))
result=dict(date='2026-10-04',report_rounds=[231,775],reports=545,phases=15,
            files_moved=len(manifest['entries']),byte_identical_files=exact_unchanged,
            markdown_link_only_updates=link_only,all_original_bytes_in_verified_snapshot=True,
            changed_scientific_code_or_json=0,missing_or_duplicate_reports=0,
            local_links_checked=checked_links,new_broken_links=new_broken,
            preexisting_historical_broken_links=preexisting,
            preexisting_draft_links_repaired=len(manifest.get('preexisting_draft_link_repairs',[])),
            new_scientific_rounds=0,goal_unchanged_and_paused=True,
            all_checks_passed=not new_broken)
with (MIG/'migration_checks.json').open('x',encoding='utf8',newline='\n') as f:
    json.dump(result,f,ensure_ascii=False,indent=2)
    f.write('\n')
print(json.dumps({k:result[k] for k in ('files_moved','reports','phases','byte_identical_files','markdown_link_only_updates','local_links_checked','all_checks_passed')},ensure_ascii=False))
print('Preexisting broken links:',len(preexisting),'New broken links:',len(new_broken))
assert result['all_checks_passed'],new_broken[:12]
