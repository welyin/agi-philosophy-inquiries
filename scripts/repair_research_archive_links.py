"""Resolve directory links and unambiguous root-relative links in frozen drafts."""
from collections import defaultdict
import json
import os
from pathlib import Path
import zipfile
import organize_research_231_775 as m

ROOT,BASE,OLD=m.ROOT,m.RESEARCH,m.OLD
MIG=BASE/'archive_764_775/_migration'
manifest=json.loads((MIG/'manifest.json').read_text('utf8'))
snapshot=json.loads((MIG/'snapshot_checks.json').read_text('utf8'))
mapping={(OLD/e['original']).resolve():(BASE/e['destination']).resolve() for e in manifest['entries']}
groups=defaultdict(set)
for src,dst in mapping.items():
    for depth,parent in enumerate(src.parents):
        if parent==OLD:
            break
        groups[parent].add(dst.parents[depth])
for directory,targets in groups.items():
    mapping[directory]=Path(os.path.commonpath(list(targets)))
known={(ROOT/e['path']).resolve() for e in snapshot['entries']}
known.update(parent for p in list(known) for parent in p.parents if parent.is_relative_to(ROOT))
repairs=[]
changed_files=[]
m.NEW=BASE
with zipfile.ZipFile(MIG/snapshot['file']) as z:
    for e in manifest['entries']:
        if not e['rewrite_markdown_links']:
            continue
        src=OLD/e['original']
        dst=BASE/e['destination']
        assert m.digest(dst)==e['current_sha256']
        raw=z.read('research_cognition_physics/archive_231_/'+e['original'])
        localmap=dict(mapping)
        if e['kind']=='draft_or_support':
            for _,_,target,local in m.links(raw.decode('utf-8-sig')):
                resolved=(src.parent/local).resolve()
                if resolved in known or resolved.exists():
                    continue
                candidate=(OLD/local).resolve()
                if candidate in mapping:
                    localmap[resolved]=mapping[candidate]
                elif candidate in known or candidate.is_file():
                    localmap[resolved]=candidate
                else:
                    continue
                repairs.append(dict(file=e['destination'],original_target=target,
                                    resolved_to=str(localmap[resolved].relative_to(ROOT)),
                                    reason='frozen draft retained link relative to the original archive root'))
        changed,edits=m.rewrite_links(raw,src,dst,localmap)
        if changed!=dst.read_bytes():
            temp=dst.with_name(dst.name+'.link_repair_tmp')
            assert not temp.exists()
            temp.write_bytes(changed)
            os.replace(temp,dst)
            changed_files.append(e['destination'])
        e.update(current_sha256=m.sha(changed),current_bytes=len(changed),link_edits=edits)
old=MIG/'manifest_before_directory_link_repair.json'
assert not old.exists()
(MIG/'manifest.json').replace(old)
manifest['preexisting_draft_link_repairs']=repairs
m.dump(MIG/'manifest.json',manifest)
m.dump(MIG/'link_repair_checks.json',dict(repaired_files=len(changed_files),
    preexisting_draft_links_repaired=len(repairs),link_destinations_only=True,
    original_scientific_prose_code_results_unchanged=True,files=changed_files))
previous=MIG/'migration_checks.json'
if previous.exists():
    previous.replace(MIG/'migration_checks_before_link_repair.json')
print('Updated files:',len(changed_files),'Preexisting draft links repaired:',len(repairs))
