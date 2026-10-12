"""Read-only asset, historical hash, and reversible navigation check."""
from pathlib import Path
import ast, hashlib, json, re, sys
ROOT=Path(__file__).resolve().parents[4]
BASE=ROOT/'research_cognition_physics'
STAGE=BASE/'archive_1120_'
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ledger=json.loads((STAGE/'_shared/dependency_ledger.json').read_text(encoding='utf-8'))
    assert ledger['latest_accepted_round']==1156
    assert ledger['new_accepted_rounds']==list(range(1120,1157))
    assert ledger['stage_accepted_rounds']==37
    assert ledger['stage_numerical_check_groups']==144
    assert ledger['in_progress_rounds']==[] and ledger['next_round'] is None
    assert ledger['scientific_calibration_increment']==0
    assert ledger['inherited_scientific_calibration_total']==3860
    prior=0; current=0; reviews=0
    for n in range(1120,1157):
        folder=STAGE/str(n)
        value=json.loads((folder/'acceptance.json').read_text(encoding='utf-8'))
        assert value['status']=='PASS',n
        for name,expected in value.get('core_sha256',{}).items():
            path=(folder/name).resolve()
            assert path.is_relative_to(STAGE)
            assert sha(path)==expected,(n,name)
            if n==1156: current+=1
            else: prior+=1
        rv=value.get('review_sha256')
        if isinstance(rv,str):
            assert rv in {sha(q) for q in folder.glob('*.md')},('review',n)
            reviews+=1
    acc=json.loads((STAGE/'1156/acceptance.json').read_text(encoding='utf-8'))
    assert acc['status']=='PASS'
    assert sha(STAGE/'1156'/acc['review_path'])==acc['review_sha256']
    src=json.loads((STAGE/'1156/source_hashes.json').read_text(encoding='utf-8'))['files']
    for entry in src:
        assert sha(ROOT/entry['path'])==entry['sha256'],entry['path']
    closure=json.loads((HERE/'completion_record.json').read_text(encoding='utf-8'))
    for name,expected in closure['core_sha256'].items():
        assert sha(HERE/name)==expected,name
    assert sha(STAGE/'1156/acceptance.json')==closure['last_round_acceptance_sha256']
    bindings=0
    review=(HERE/'closure_review.md').read_text(encoding='utf-8')
    for path,expected in re.findall(r'^\|([^|\r\n]+)\|([0-9a-f]{64})\|$',review,re.M):
        assert sha(STAGE/path)==expected,path
        bindings+=1
    protected={
      'scripts/run_research_merged.py':'5b7914ab5cf6ceaa795ea19bdb2b1285d0025845e8e6b665662b992594c157ec',
      'research_cognition_physics/_migration/merge_1106_1119/manifest.json':'b39ac02f94a4902f1f6fa9508962496812e25ae726a36c29eeeed77a5735a18b'}
    for name,expected in protected.items(): assert sha(ROOT/name)==expected,name
    sys.path.insert(0,str(ROOT/'scripts'))
    from research_layout import relocation_original_bytes
    layer=json.loads((BASE/'_migration/field_active_1120/layer.json').read_text(encoding='utf-8'))
    frozen=json.loads((BASE/'_migration/merge_1106_1119/manifest.json').read_text(encoding='utf-8'))
    prior_manifest={e['destination']:e for e in frozen['entries']}
    for entry in layer['entries']:
        current_bytes=(ROOT/entry['destination']).read_bytes()
        assert hashlib.sha256(current_bytes).hexdigest()==entry['current_sha256']
        original=relocation_original_bytes(current_bytes,entry)
        assert hashlib.sha256(original).hexdigest()==prior_manifest[entry['original']]['current_sha256']
    docs=list((STAGE/'1156').glob('*.md'))+list(HERE.glob('*.md'))+[
        STAGE/'research_note_1156.md',STAGE/'README.md',STAGE/'文件索引.md',
        STAGE/'_shared/after1155_joint_participant_model_screen.md',
        STAGE/'_shared/after1155_axiom_implication_screen.md']
    links=0
    for path in docs:
        text=path.read_text(encoding='utf-8')
        assert text.count('$$')%2==0,('display math delimiters',path)
        for target in re.findall(r'\]\(([^)\r\n]+)\)',text):
            if target.startswith(('http:','https:','mailto:','#')): continue
            target=target.strip('<>').split('#')[0]
            assert (path.parent/target).resolve().exists(),(path,target)
            links+=1
    for name in ['1156/check.py','_shared/update_navigation.py','_shared/stage_closure/verify_assets.py']:
        ast.parse((STAGE/name).read_text(encoding='utf-8'))
    return dict(status='PASS',previous_core_hashes=prior,current_core_hashes=current,
         formal_review_hashes=reviews,current_review_hash=True,frozen_sources=len(src),closure_review_bindings=bindings,
         closure_core_hashes=len(closure['core_sha256']),protected_files=len(protected),
         reversible_navigation_files=len(layer['entries']),local_links_checked=links,
         markdown_files_checked=len(docs),python_ast_checks=3,
         accepted_rounds=37,numerical_check_groups=144,new_scientific_calibrations=0,
         scope='Asset checks only; analytic and numeric claims use signed round documents.')
if __name__=='__main__': print(json.dumps(main(),ensure_ascii=False,indent=2))
