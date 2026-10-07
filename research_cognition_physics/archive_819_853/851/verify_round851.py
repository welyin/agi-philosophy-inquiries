"""Reproduce 851 and verify frozen history, citations and scope receipts."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import full_source_mean_content as experiment
RECEIPT=HERE/'research_round_851_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    assert experiment.run()==json.loads(experiment.TARGET.read_text('utf-8'))
    for n in range(776,851):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    audit=json.loads((HERE/'drafts/common_model_stage_audit_receipt.json').read_text('utf-8'))
    assert sha(HERE/'drafts/common_model_stage_audit.md')==audit['audit_sha256']
    for name,digest in audit['reference_files'].items():assert sha(ROOT/name)==digest,name
    draft=json.loads((HERE/'drafts/stage_paper_working_checks.json').read_text('utf-8'))
    assert sha(ROOT/draft['draft_path'])==draft['draft_sha256']
    history=Layout().verify()
    note=STAGE/'research_note_851.md';body=note.read_text('utf-8')
    assert body.count('$$')==20 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,11)]
    docs=[note,STAGE/'852/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    docs += list((HERE/'drafts').glob('*.md'))
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
    assert sorted(n for n in numbers if n<=851)==list(range(1,852))
    frozen=[STAGE/'850/research_round_850_checks.json',STAGE/'research_note_850.md',
        STAGE/'research_note_796.md',STAGE/'research_note_844.md',STAGE/'research_note_848.md',
        STAGE/'research_note_765.md',STAGE/'769/drafts/research_note_769_working.md',
        HERE/'drafts/STATUS.md',*sorted((HERE/'drafts').glob('*.json')),
        HERE/'drafts/common_model_stage_audit.md',ROOT/draft['draft_path']]
    fresh=[note,STAGE/'852/drafts/STATUS.md',HERE/'drafts/research_note_851_working.md',
        *sorted(HERE.glob('*.py')),*sorted(p for p in HERE.glob('*.json') if p!=RECEIPT)]
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=851,date='2026-10-05',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_850=3636,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_850_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=851,local_links_checked=links,display_equations=10,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='At the original classical on-shell background, matching compatible complete-source first-order Cauchy means by a common full homogeneous displacement preserves the same thermal organized CAR input, physical formal positivity and Ward, and the original nonzero hbar-cubed joint-field content cumulant. It does not re-quantize about an off-shell mean background or establish finite coupling, autonomous apparatus or full feedback.',
        common_displacement_independent_of_record_content=True,
        original_nonzero_joint_content_coefficient_preserved_analytically=True,
        free_fermionic_input_entropy_domain_unchanged=True,
        complete_source_or_original_continuum_PDE_computed_numerically=False,
        arbitrary_local_source_extended_globally_without_proof=False,
        sourced_response_itself_used_as_homogeneous_displacement=False,
        generic_semiclassical_balance_implies_free_gauge_complex=False,
        quantum_theory_rebuilt_on_nonclassical_background=False,
        all_higher_content_coefficients_unchanged_claimed=False,
        spacetime_group_action_remain_inputs=True,full_goal_completed=False,
        visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite a frozen receipt.'
    result=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert old[k]==result[k],k
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
