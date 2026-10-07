"""Verify physical ownership, reversible references and unchanged scientific bytes."""
from pathlib import Path
import argparse
import hashlib
import json
import re

from organize_research_231_775 import links
from research_layout import Layout
from split_research_764_1008 import ROOT, BASE, OLD, MIG, documents, rel


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def restore(raw, entry):
    assert sha(raw)==entry['current_sha256'],entry['destination']
    if not entry['edits']:
        assert sha(raw)==entry['original_sha256']
        return raw
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    text=raw.decode(enc); offset=0; previous=0; pieces=[]
    for e in entry['edits']:
        start=e['start']+offset; end=start+len(e['after'])
        assert start>=previous and text[start:end]==e['after'],entry['destination']
        assert e['end']-e['start']==len(e['before'])
        pieces.extend((text[previous:start],e['before']));previous=end
        offset+=len(e['after'])-len(e['before'])
    pieces.append(text[previous:]); original=''.join(pieces).encode(enc)
    assert sha(original)==entry['original_sha256'],entry['original']
    return original


def verify(writing=False):
    data=json.loads((MIG/'manifest.json').read_text('utf8'))
    assert not OLD.exists()
    assert len(data['phases'])==10
    originals=set(); destinations=set(); science=md_changed=0
    for e in data['entries']:
        assert e['original'] not in originals and e['destination'] not in destinations
        originals.add(e['original']);destinations.add(e['destination'])
        target=(ROOT/e['destination']).resolve()
        assert target.is_relative_to(ROOT.resolve())
        raw=target.read_bytes(); before=restore(raw,e)
        if target.suffix in ('.py','.json'):
            assert raw==before and not e['edits'],e['destination']
            science+=1
        if raw!=before:
            assert target.suffix.lower()=='.md';md_changed+=1
        if e['original']!=e['destination']:
            assert not (ROOT/e['original']).exists(),e['original']
    phase_reports=[]
    for p in data['phases']:
        folder=BASE/f"archive_{p['start']}_{p['end']}"
        numbers=sorted(int(re.fullmatch(r'research_note_(\d+)\.md',f.name)[1])
            for f in folder.iterdir() if re.fullmatch(r'research_note_(\d+)\.md',f.name))
        assert numbers==list(range(p['start'],p['end']+1)),folder
        for n in numbers:
            assert (folder/str(n)).is_dir(),n
        assert (folder/'README.md').is_file() and (folder/'文件索引.md').is_file()
        for f in folder.iterdir():
            if f.is_dir():assert f.name.isdigit() and int(f.name) in numbers,f
            else:assert f.suffix=='.md' and (f.name in ('README.md','文件索引.md') or f.name.startswith('research_note_')),f
        phase_reports+=numbers
    assert phase_reports==list(range(764,1009))
    allnotes=[]
    for f in BASE.glob('archive_*/research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+)\.md',f.name)
        if m:allnotes.append(int(m[1]))
    assert sorted(allnotes)==list(range(1,1009))
    inherited={tuple(item) for item in data['inherited_missing_links']}
    checked=0; old_missing=[]; new_missing=[]
    imminent={MIG/'checks.json',MIG/'replay_checks.json'} if writing else set()
    for doc in documents():
        prose=doc.read_text('utf-8-sig')
        for _,_,target,local in links(prose):
            path=(doc.parent/local.replace('\\','/')).resolve()
            checked+=1
            if path.exists() or path in imminent:continue
            pair=(rel(doc),target)
            (old_missing if pair in inherited else new_missing).append(pair)
    assert not new_missing,('New missing links',new_missing[:20],len(new_missing))
    layout=Layout()
    historic=layout.verify()
    relocation=layout.verify_relocation()
    assert historic==dict(mapped_files=9965,checked_files=7146,unchanged_files=5512,
        verified_link_only_changes=1634,optional_local_files_skipped=2819,zip_files_read=0)
    return dict(date='2026-10-08',all_migration_checks_passed=True,phase_count=10,
        archived_rounds=[764,1008],formal_reports_in_split=245,all_formal_reports=1008,
        mapped_entries=len(data['entries']),unchanged_python_and_json= science,
        markdown_files_with_reversible_edits=md_changed,local_links_checked=checked,
        preexisting_missing_links=old_missing,new_missing_links=[],
        historical_layout=historic,relocation_checks=relocation,
        source_directory_removed=True,zip_created=False,duplicate_tree_created=False,
        new_scientific_rounds=0,goal_created_or_changed=False,visual_checks_performed=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=verify(args.write)
    if args.write:
        with (MIG/'checks.json').open('x',encoding='utf8') as out:
            json.dump(result,out,ensure_ascii=False,indent=2);out.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='preexisting_missing_links'},ensure_ascii=False,indent=2))
