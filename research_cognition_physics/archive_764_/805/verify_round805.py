"""Reproduce 805, validate navigation, preserve all prior frozen science."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import original_past_covariance_action as past
import galerkin_residual_probe as residual
RECEIPT=HERE/'research_round_805_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    for ex in (past,residual):assert ex.run()==json.loads(ex.TARGET.read_text('utf-8'))
    for n in range(776,805):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for section in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[section].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_805.md';body=note.read_text('utf-8')
    assert body.count('$$')==20
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,11)]
    docs=[note,STAGE/'806/drafts/STATUS.md',HERE/'drafts/research_note_805_working.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=805)==list(range(1,806))
    frozen=[STAGE/'804/research_round_804_checks.json',STAGE/'research_note_804.md',
            STAGE.parent/'archive_702_741/research_note_730.md',HERE/'drafts/STATUS.md']
    fresh=[note,STAGE/'806/drafts/STATUS.md',*sorted(HERE.glob('*.py')),
        *sorted(p for p in HERE.glob('*.json') if p!=RECEIPT),
        HERE/'drafts/research_note_805_working.md',HERE/'drafts/working_checks.json']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=805,date='2026-10-05',all_checks_passed=True,fresh_test_groups=2,
        cumulative_numbered_test_groups_from_804=3583,saved_results_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_804_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=805,local_links_checked=links,display_equations=10,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='The original 730 auxiliary-past covariance is an explicit full 64-dimensional Fourier multiplier. Its action on fixed modes is backward propagation, multiplier action, and forward propagation in the original auxiliary background. Initial-state, full PDE residual, coefficient and quadrature errors propagate to the original finite response matrix. Future joint background and bosonic covariance values have not been numerically supplied.',
        original_past_multiplier_implemented=True,instantaneous_target_reference_substituted=False,
        all_original_fermion_coefficients_retained=True,
        continuous_original_response_value_evaluated=False,original_nonzero_response_proven=False,
        original_future_background_numerically_solved=False,
        autonomous_apparatus_or_finite_coupling_proven=False,
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
