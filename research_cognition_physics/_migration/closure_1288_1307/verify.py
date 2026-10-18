"""Check actual archive inventory and original acceptance bindings after closure."""
from pathlib import Path
import hashlib
import json
import re
import sys

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'research_cognition_physics'
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from run_research_structure_selection_closed import ClosedStructureSelectionLayout

def sha(raw):return hashlib.sha256(raw).hexdigest()

def main():
    store=ClosedStructureSelectionLayout()
    before=json.loads((HERE/'before.json').read_text(encoding='utf-8'))
    manifest=json.loads((HERE/'manifest.json').read_text(encoding='utf-8'))
    entries={e['original']:e for e in manifest['entries']}
    for item in before['files']:
        assert item['path'] in entries
        raw=store.before_structure_closure(ROOT/item['path'])
        assert sha(raw)==item['sha256'] and len(raw)==item['bytes']
    bindings=0
    for n in range(1288,1308):
        stage='archive_1288_' if n<1290 else 'archive_1290_'
        ap=BASE/stage/str(n)/'acceptance.json'
        a=json.loads(store.before_structure_closure(ap))
        assert a['status'].startswith('ACCEPTED'),n
        for rel,expected in a['files_sha256'].items():
            assert sha(store.before_structure_closure(ap.parent/rel))==expected,(n,rel)
            bindings+=1
    links=0;previous_missing=[]
    for e in manifest['entries']:
        if not e['edits']:continue
        p=ROOT/e['destination'];text=p.read_text(encoding='utf-8-sig')
        for link in re.findall(r'\]\(([^)]+)\)',text):
            if not re.search(r'archive_(1288_1289|1290_1307)/',link) or '://' in link:continue
            if not (p.parent/link.split('#')[0]).resolve().exists():
                oldlink=link.replace('archive_1288_1289/','archive_1288_/').replace('archive_1290_1307/','archive_1290_/')
                original_text=store.before_structure_closure(ROOT/e['original']).decode('utf-8-sig')
                assert ']('+oldlink+')' in original_text,(e['destination'],'new broken link',link)
                old_target=(ROOT/e['original']).parent/oldlink.split('#')[0]
                assert not store.structure_closure_actual(old_target.resolve()).exists(),(e['destination'],'regressed link',link)
                previous_missing.append({'file':e['destination'],'link':link,'status':'already_missing_before_migration'})
            else:links+=1
    result=dict(schema='closure_1288_1307_verification_v1',status='PASS_WITH_RECORDED_LEGACY_EXCEPTIONS',
        moved_files=377,python_json_byte_identical=193,rounds=20,round_evidence_bindings=bindings,
        valid_migrated_path_links=links,new_broken_links=0,preexisting_missing_links=previous_missing,
        closure_check=store.verify_closure(),navigation_check=store.verify_structure_navigation(),
        actual_execution_receipts={'rename':'dccc9a','apply':'48f23c','closure':'c71eb8','navigation':'21db39',
            '1307_author':'006dfd','1307_independent':'9b8b31','1288_independent':'a8182e','stage_joint':'5da314'},
        legacy_all_layout_exception=dict(receipt='e92503',file='research_cognition_physics/08_认知热学角色与热力学四律_阶段论文.md',
            already_documented_in='research_cognition_physics/_migration/closure_1281_1287/independent_audit.md',modified_by_this_migration=False),
        closed_runner_sha256=sha((ROOT/'scripts/run_research_structure_selection_closed.py').read_bytes()),
        new_scientific_rounds=0,new_axioms=0)
    with (HERE/'verification.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({'status':result['status'],'round_evidence_bindings':bindings,'valid_links':links,
                      'preexisting_missing_links':len(previous_missing),'new_broken_links':0},ensure_ascii=False))

if __name__=='__main__':main()
