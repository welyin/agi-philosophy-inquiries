"""Read-only validation of the reversible 1313--1322 archive closure."""
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
from run_research_interaction_generation_closed import ClosedInteractionGenerationLayout
from run_research_conditional_selection_closed import ClosedConditionalSelectionLayout

def sha(raw): return hashlib.sha256(raw).hexdigest()

def layout():
    if (BASE/'_migration/interaction_rules_active_1323/layer.json').is_file():
        from run_research_interaction_rules_active import InteractionRulesLayout
        return InteractionRulesLayout()
    return ClosedInteractionGenerationLayout()

def migration_current(original,entry):
    enc='utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf-8'
    text=original.decode(enc)
    for edit in reversed(entry['edits']):
        assert text[edit['start']:edit['end']]==edit['before']
        text=text[:edit['start']]+edit['after']+text[edit['end']:]
    raw=text.encode(enc)
    assert sha(raw)==entry['current_sha256']
    return raw

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-new',action='store_true')
    args=parser.parse_args()
    store=layout()
    before=json.loads((HERE/'before.json').read_text(encoding='utf-8'))
    manifest=json.loads((HERE/'manifest.json').read_text(encoding='utf-8'))
    entries={e['original']:e for e in manifest['entries']}
    scientific=0
    for item in before['files']:
        entry=entries[item['path']]
        raw=store.before_interaction_closure(ROOT/item['path'])
        assert sha(raw)==item['sha256'] and len(raw)==item['bytes'],item['path']
        if Path(item['path']).suffix.lower() in ('.py','.json'):
            assert (ROOT/entry['destination']).read_bytes()==raw
            scientific+=1
    assert scientific==before['scientific_python_json']==119
    bindings=0
    for number in range(1313,1323):
        old=BASE/'archive_1313_'/str(number)/'acceptance.json'
        acceptance=json.loads(store.before_interaction_closure(old))
        assert acceptance['status'].startswith('ACCEPTED'),number
        for name,expected in acceptance['sha256'].items():
            assert sha(store.before_interaction_closure(ROOT/name))==expected,(number,name)
            bindings+=1
    assert bindings==before['round_bindings']==91
    stage_raw=store.before_interaction_closure(BASE/'archive_1313_/_shared/stage_acceptance_after1322.json')
    assert sha(stage_raw)==before['stage_acceptance_sha256']
    stage=json.loads(stage_raw)
    assert stage['status']=='ACCEPTED_DECLARED_STAGE_SCOPE'
    assert stage['scientific_rounds']==list(range(1313,1323))
    for name,expected in stage['sha256'].items():
        assert sha(store.before_interaction_closure(ROOT/name))==expected,('stage',name)
    assert len(stage['sha256'])==before['stage_bindings']==105
    links=0
    previous_missing=[]
    for entry in manifest['entries']:
        if not entry['edits']:continue
        original=store.before_interaction_closure(ROOT/entry['original'])
        current=migration_current(original,entry).decode('utf-8-sig')
        path=ROOT/entry['destination']
        for link in re.findall(r'\]\(([^)]+)\)',current):
            if 'archive_1313_1322/' not in link or '://' in link:continue
            target=(path.parent/link.split('#')[0].strip('<>')).resolve()
            if target.exists():
                links+=1
            else:
                oldlink=link.replace('archive_1313_1322/','archive_1313_/')
                assert ']('+oldlink+')' in original.decode('utf-8-sig'),(entry['destination'],'new broken link',link)
                oldtarget=((ROOT/entry['original']).parent/oldlink.split('#')[0].strip('<>')).resolve()
                assert not store.interaction_closure_actual(oldtarget).exists()
                previous_missing.append(dict(file=entry['destination'],link=link))
    result=dict(schema='closure_1313_1322_verification_v1',status='PASS',
        moved_files=len(before['files']),python_json_byte_identical=scientific,
        rounds=10,round_evidence_bindings=bindings,stage_evidence_bindings=len(stage['sha256']),
        valid_migrated_path_links=links,new_broken_links=0,preexisting_missing_links=previous_missing,
        closure_check=ClosedInteractionGenerationLayout.verify_closure(store),
        navigation_check=store.verify_interaction_navigation(),
        previous_closure_check=ClosedConditionalSelectionLayout.verify_closure(store),
        actual_execution_receipts=dict(rename='d71a97',apply='d9d02d',closure='f364cc',
            interaction_navigation='6cbd34',previous_closure='a27c61',
            round1322_author_old_path='cd678c',round1322_independent_new_path='6b7cb8',
            round1321_author_old_path='5cf2e1'),
        closed_runner_sha256=sha((ROOT/'scripts/run_research_interaction_generation_closed.py').read_bytes()),
        legacy_all_layout_exception=dict(
            file='research_cognition_physics/08_认知热学角色与热力学四律_阶段论文.md',
            already_documented_in='research_cognition_physics/_migration/closure_1281_1287/independent_audit.md',
            this_migration_modified_file=False,all_historical_layouts_claimed_pass=False),
        new_scientific_rounds=0,new_axioms=0,directory_copy_used=False,zip_used=False)
    if args.write_new:
        with (HERE/'verification.json').open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write('\n')
    else:
        saved=json.loads((HERE/'verification.json').read_text(encoding='utf-8'))
        assert result==saved,'Migration verification changed'
    print(json.dumps(dict(status='PASS',moved_files=len(before['files']),
        python_json=scientific,round_bindings=bindings,stage_bindings=len(stage['sha256']),
        valid_links=links,previous_missing_links=len(previous_missing),
        mode='write-new' if args.write_new else 'read-only')))

if __name__=='__main__':main()
