"""913 working evidence preservation, not physical accuracy certification."""
from pathlib import Path
import argparse,ast,hashlib,json,re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'drafts/working_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text('utf-8-sig'))
def run(writing=False):
    formal=read(STAGE/'912/research_round_912_checks.json')
    for key in ('frozen_inputs','new_scientific_and_entry_files'):
        for p,s in formal[key].items():assert sha(ROOT/p)==s,p
    data=read(HERE/'receiver_response_readout_results.json')
    for p,s in data['source_hashes'].items():assert sha(STAGE/p)==s,p
    assert data['formal_previous']==912 and data['cumulative_previous']==3697
    assert [r['N'] for r in data['rows']]==[17,25]
    for row in data['rows']:
        assert row['original_source_and_preparation_unchanged']
        assert row['full_pairing']['jet_pairing_identity_residual']<1e-16
        for a,b,c in zip(row['full_pairing']['total'],row['full_pairing']['direct'],row['full_pairing']['reference']):assert abs(a-b-c)<1e-16
        assert abs(row['full_pairing']['direct'][2])<1e-16
        assert row['source_stats']['time_min']>-row['T'] and row['source_stats']['time_max']<row['T']
        assert row['finite_response_modes']==row['N']**3
        assert not row['physical_A_certified'] and not row['actual_continuous_pairing_error_certified']
        assert not row['full_vA_field_computed'] and not row['full872_response_computed']
    diffs=[abs(a-b) for a,b in zip(data['rows'][0]['full_pairing']['total'],data['rows'][1]['full_pairing']['total'])]
    assert diffs==data['mesh_difference_not_error_bound']
    files=[Path(__file__),HERE/'receiver_response_readout.py',HERE/'receiver_response_readout_results.json',HERE/'drafts/receiver_readout_working.md',HERE/'drafts/STATUS.md']
    docs=[HERE/'drafts/receiver_readout_working.md',HERE/'drafts/STATUS.md']
    docs += [STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    docs += [STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    links=0
    for doc in docs:
        body=doc.read_text('utf-8-sig');assert body.count('$$')%2==0
        for link in re.findall(r'\]\(([^)]+)\)',re.sub(r'\$\$.*?\$\$','',body,flags=re.S)):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    for p in files:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    hashes={str(p.relative_to(ROOT)):sha(p) for p in files}
    if not writing:assert read(TARGET)['files']==hashes
    return dict(date='2026-10-06',checks_passed=True,formal_rounds=912,cumulative_numbered_groups=3697,new_scientific_groups=0,
        current_round_status='913_working',actual_common_stress_response_paired=True,
        actual_finite_observable_error_certified=False,full_vA_field_computed=False,full872_response_computed=False,
        local_links_checked=links,files=hashes,preserved_formal912_receipt_sha256=sha(STAGE/'912/research_round_912_checks.json'),full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not TARGET.exists()
    result=run(args.write)
    if args.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='files'},ensure_ascii=False,indent=2))
