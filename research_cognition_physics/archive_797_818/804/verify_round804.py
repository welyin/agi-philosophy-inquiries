"""Reproduce 804 diagnostics, verify its links and preserve frozen science."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import covariance_contract_probe as covariance
import compact_response_certificate as certificate
RECEIPT=HERE/'research_round_804_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    for ex in (covariance,certificate):assert ex.run()==json.loads(ex.TARGET.read_text('utf-8'))
    for n in range(776,804):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for section in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[section].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_804.md';body=note.read_text('utf-8')
    assert body.count('$$')==22
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,12)]
    assert not any(ord(c)<32 and c not in '\r\n' for c in body)
    docs=[note,STAGE/'805/drafts/STATUS.md',HERE/'drafts/research_note_804_working.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    links=0
    for doc in docs:
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text('utf-8-sig')):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            path=(doc.parent/target.split('#')[0].strip('<>')).resolve()
            assert path.exists() or (writing and path==RECEIPT.resolve()),(str(doc),target)
            links+=1
    numbers=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if p.parent.name.startswith('archive_') and m:numbers.append(int(m[1]))
    assert sorted(n for n in numbers if n<=804)==list(range(1,805))
    frozen=[STAGE/'803/research_round_803_checks.json',STAGE/'research_note_803.md',
        STAGE.parent/'archive_702_741/research_note_734.md',STAGE.parent/'archive_702_741/research_note_741.md',HERE/'drafts/STATUS.md']
    fresh=[note,STAGE/'805/drafts/STATUS.md',*sorted(HERE.glob('*.py')),
        *sorted(p for p in HERE.glob('*.json') if p!=RECEIPT),HERE/'drafts/research_note_804_working.md']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=804,date='2026-10-05',all_checks_passed=True,fresh_test_groups=2,
        cumulative_numbered_test_groups_from_803=3581,saved_results_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_803_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=804,local_links_checked=links,display_equations=11,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='Reuse 734/741 smooth terminal history variation to prove that 803 ordered cross response is smoothing and its local sterile compression compact self-adjoint. A fixed finite bosonic menu has no uniform positive response bound over all permitted normalized smooth local record modes. Existence of a nonzero mode is equivalent to a nonzero compressed operator; a certified finite continuous block suffices for a witness, whereas a global optimum also needs a tail bound. The occupied/greater covariance dictionary preserves the state.',
        fixed_original_state_preserved=True,compact_record_response_proven_under_inherited_contract=True,
        continuous_original_response_value_evaluated=False,original_nonzero_response_proven=False,
        diagnostic_is_original_PDE=False,finite_coupling_or_autonomous_apparatus_proven=False,
        graph_to_continuum_map_proven=False,spacetime_group_action_remain_inputs=True,
        original_all_physics_unification_complete=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite a frozen receipt.'
    r=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert old[k]==r[k],k
    print(json.dumps({k:v for k,v in r.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
