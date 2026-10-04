"""Check all numbered reports, untouched evidence, and the live local links."""
import argparse
import ast
from collections import Counter
import json
from pathlib import Path
import re
import zipfile
import organize_research_231_775 as m

ROOT,BASE=m.ROOT,m.RESEARCH
LAST=BASE/'archive_764_'
MAINT=LAST/'_migration/layout_001_230_20261004'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--report',type=Path,help='Write a new receipt instead of replacing the historical layout report')
args=parser.parse_args()


def read(p):
    return json.loads(p.read_text('utf8'))


def evidence_check(manifest, snapshot, prefix, path_base):
    assert m.digest(snapshot['path'])==snapshot['sha256']
    identical=links_only=0
    documents=[]
    raw_docs=set()
    with zipfile.ZipFile(snapshot['path']) as z:
        for e in manifest['entries']:
            path=path_base/e['destination']
            original=z.read(prefix+e['original'])
            current=path.read_bytes()
            assert m.sha(original)==e['original_sha256'] and len(original)==e['bytes']
            assert m.sha(current)==e['current_sha256'] and len(current)==e['current_bytes'],e['destination']
            if original!=current:
                assert e['rewrite_markdown_links'] and path.suffix=='.md',e['destination']
                enc='utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
                text=original.decode(enc)
                spans={(a,b):t for a,b,t,_ in m.links(text)}
                for edit in reversed(e['link_edits']):
                    assert spans[edit['start'],edit['end']]==edit['before']
                    text=text[:edit['start']]+edit['after']+text[edit['end']:]
                assert text.encode(enc)==current,e['destination']
                links_only+=1
            else:
                identical+=1
            if path.suffix=='.md':
                if e['rewrite_markdown_links']:
                    documents.append(path)
                else:
                    raw_docs.add(path.resolve())
    return dict(files=len(manifest['entries']),byte_identical=identical,link_only=links_only),documents,raw_docs


early,late=read(MAINT/'manifest.json'),read(LAST/'_migration/manifest.json')
early_snapshot=read(MAINT/'snapshot.json')
late_snapshot=read(LAST/'_migration/snapshot_checks.json')
ec,ed,er=evidence_check(early,dict(path=MAINT/early_snapshot['file'],sha256=early_snapshot['sha256']),'',ROOT)
lc,ld,lr=evidence_check(late,dict(path=LAST/'_migration'/late_snapshot['file'],sha256=late_snapshot['sha256']),'research_cognition_physics/archive_231_/',BASE)
split=read(MAINT/'phase_split_plan.json')
frozen=er|lr|{(ROOT/p).resolve() for p in split['raw_historical_navigation']}
all_archives=[BASE/p['directory'] for p in early['early_phases']]+[BASE/'archive_223_230']+[BASE/p['directory'] for p in late['phases']]
assert len(all_archives)==27 and len(set(all_archives))==27
reports=[]
for archive in all_archives:
    assert archive.is_dir() and (archive/'README.md').is_file() and (archive/'文件索引.md').is_file()
    for path in archive.glob('research_note_*.md'):
        match=re.fullmatch(r'research_note_(\d+)\.md',path.name)
        assert match,path
        n=int(match[1]);reports.append(n)
        assert (archive/str(n)).is_dir(),n
assert sorted(reports)==list(range(1,776)),Counter(reports)
assert all(not (BASE/name).exists() for name in ('archive_001_222','archive_764_775','archive_231_','archive_231_775'))
assert sorted(p.name for p in BASE.iterdir() if p.is_file())==sorted(['README.md','RESEARCH_STATE.md','research_direction.md','可组合认知结构与复量子状态空间_阶段论文.md','可组合认知结构与有限维量子理论_阶段论文.md'])
for row in read(MAINT/'round_index_001_222.json')['rounds']:
    assert (BASE/row['note']).is_file(),row
    assert m.digest(BASE/row['note'])==row['current_note_sha256']
    for artifact in row['indexed_artifacts']:
        assert (BASE/artifact).is_file(),(row['round'],artifact)
documents=set(ed+ld)
for archive in all_archives:
    for p in archive.rglob('*.md'):
        if p.resolve() in frozen or any(part in ('navigation_before_final_edit','navigation_before_phase_summary') for part in p.parts):
            continue
        documents.add(p)
documents.update([ROOT/'README.md']+list(BASE.glob('*.md')))
broken=[];total=0
out=args.report.resolve() if args.report else MAINT/'layout_checks.json'
assert out.resolve().is_relative_to(ROOT.resolve())
for path in sorted(documents):
    for _,_,target,local in m.links(path.read_text('utf-8-sig')):
        dest=(path.parent/local).resolve()
        total+=1
        if dest==out.resolve():
            continue
        if not dest.exists():
            broken.append(dict(file=path.relative_to(ROOT).as_posix(),target=target))
for p in [ROOT/'scripts/organize_research_001_230.py',ROOT/'scripts/finish_research_early_navigation.py',ROOT/'scripts/split_research_001_222_phases.py',ROOT/'scripts/verify_research_final_layout.py',MAINT/'replay_early.py',LAST/'_migration/replay.py']:
    ast.parse(p.read_text('utf8'))
replay=read(MAINT/'replay_checks.json')
assert replay['all_tests_passed'] and sum(j['tests_run'] for j in replay['jobs'])==1985
result=dict(date='2026-10-04',early_phases=11,total_archive_directories=27,numbered_reports=775,
            report_numbers_complete_and_unique=True,early_evidence=ec,later_evidence=lc,
            unchanged_scientific_code_results_and_receipts=True,markdown_changes_only_link_targets=True,
            checked_markdown_documents=len(documents),local_links_checked=total,broken_links=broken,
            original_bytes_in_verified_snapshots=True,prior_historical_contracts_preserved=True,
            early_unittest_cases_passed=1985,active_directory='archive_764_',next_round=776,
            new_scientific_rounds=0,all_checks_passed=not broken)
if 'major_stage' in early:
    stage=early['major_stage']
    assert stage['rounds']==[1,230]
    assert [n for p in stage['subphases'] for n in range(p['rounds'][0],p['rounds'][1]+1)]==list(range(1,231))
    assert len(stage['subphases'])==12 and (BASE/stage['overview']).is_file()
    removed=early.get('canonical_reference_changes',[])
    assert len(removed)==2
    for reference in removed:
        assert not (ROOT/reference['deleted_copy']).exists()
        assert m.digest(ROOT/reference['canonical'])==reference['canonical_sha256']
    result.update(major_stage_rounds=[1,230],major_stage_subphases=12,
                  deleted_reference_duplicates=2,canonical_reference_files_unchanged=True,
                  scientific_tests_rerun=False,earlier_test_report_retained=True)
assert not out.exists(), 'Preserve previous receipt before rerunning'
m.dump(out,result)
print(json.dumps(result,ensure_ascii=False,indent=2))
assert not broken,broken[:20]
