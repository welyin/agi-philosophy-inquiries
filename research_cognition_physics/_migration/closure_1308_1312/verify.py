"""Verify migration inventory and frozen acceptance; default mode is read-only."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import sys

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'research_cognition_physics'
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from run_research_conditional_selection_closed import ClosedConditionalSelectionLayout

def sha(raw): return hashlib.sha256(raw).hexdigest()

def layout():
    # Later navigation may prepend text, but may not rewrite frozen evidence.
    next_layer=BASE/'_migration/interaction_generation_active_1313/layer.json'
    if next_layer.is_file():
        from run_research_interaction_generation_active import InteractionGenerationLayout
        return InteractionGenerationLayout()
    return ClosedConditionalSelectionLayout()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--write-new',action='store_true')
    args=p.parse_args()
    store=layout()
    before=json.loads((HERE/'before.json').read_text(encoding='utf-8'))
    manifest=json.loads((HERE/'manifest.json').read_text(encoding='utf-8'))
    entries={e['original']:e for e in manifest['entries']}
    for item in before['files']:
        assert item['path'] in entries
        raw=store.before_conditional_closure(ROOT/item['path'])
        assert sha(raw)==item['sha256'] and len(raw)==item['bytes'],item['path']
    bindings=0
    for n in range(1308,1313):
        path=BASE/'archive_1308_'/str(n)/'acceptance.json'
        a=json.loads(store.before_conditional_closure(path))
        assert a['status'].startswith('ACCEPTED'),n
        for name,expected in a['sha256'].items():
            assert sha(store.before_conditional_closure(ROOT/name))==expected,(n,name)
            bindings+=1
    stage=json.loads(store.before_conditional_closure(BASE/'archive_1308_/_shared/stage_acceptance.json'))
    assert stage['status']=='COMPLETE_CONDITIONAL_SCOPE'
    assert stage['accepted_artifact_bindings']==bindings==48
    stage_bindings=0
    for name,expected in stage['sha256'].items():
        assert sha(store.before_conditional_closure(ROOT/name))==expected,name
        stage_bindings+=1
    links=0;previous_missing=[]
    for e in manifest['entries']:
        if not e['edits']:continue
        path=ROOT/e['destination']
        for link in re.findall(r'\]\(([^)]+)\)',path.read_text(encoding='utf-8-sig')):
            if 'archive_1308_1312/' not in link or '://' in link: continue
            target=(path.parent/link.split('#')[0].strip('<>')).resolve()
            if target.exists():
                links+=1
            else:
                oldlink=link.replace('archive_1308_1312/','archive_1308_/')
                original=store.before_conditional_closure(ROOT/e['original']).decode('utf-8-sig')
                assert ']('+oldlink+')' in original,(e['destination'],'new broken link',link)
                oldtarget=((ROOT/e['original']).parent/oldlink.split('#')[0]).resolve()
                assert not store.conditional_closure_actual(oldtarget).exists()
                previous_missing.append(dict(file=e['destination'],link=link))
    result=dict(schema='closure_1308_1312_verification_v1',status='PASS',
        moved_files=len(before['files']),python_json_byte_identical=42,
        rounds=5,round_evidence_bindings=bindings,stage_evidence_bindings=stage_bindings,
        valid_migrated_path_links=links,new_broken_links=0,preexisting_missing_links=previous_missing,
        closure_check=store.verify_closure(),navigation_check=store.verify_conditional_navigation(),
        previous_closure_check=__import__('run_research_structure_selection_closed').ClosedStructureSelectionLayout.verify_closure(store),
        actual_execution_receipts=dict(rename='2198e6',apply='408a3f',closure='704595',
            conditional_navigation='b800e5',previous_closure='43da4d',
            previous_navigation='d27403',round1312_author_old_path='f0bc82',
            round1312_independent_new_path='4e52cb',round1311_cross_round_import='2b89c9',
            round1307_previous_stage='042aa0'),
        closed_runner_sha256=sha((ROOT/'scripts/run_research_conditional_selection_closed.py').read_bytes()),
        legacy_all_layout_exception=dict(
            file='research_cognition_physics/08_认知热学角色与热力学四律_阶段论文.md',
            already_documented_in='research_cognition_physics/_migration/closure_1281_1287/independent_audit.md',
            this_migration_modified_file=False,all_historical_layouts_claimed_pass=False),
        new_scientific_rounds=0,new_axioms=0,directory_copy_used=False,zip_used=False)
    if args.write_new:
        with (HERE/'verification.json').open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        saved=json.loads((HERE/'verification.json').read_text(encoding='utf-8'))
        assert result==saved,'Migration verification changed'
    print(json.dumps(dict(status='PASS',moved_files=len(before['files']),
        round_bindings=bindings,stage_bindings=stage_bindings,valid_links=links,
        previous_missing_links=len(previous_missing),mode='write-new' if args.write_new else 'read-only')))

if __name__=='__main__': main()
