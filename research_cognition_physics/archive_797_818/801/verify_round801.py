"""Replay 801 and check frozen science, numbering and the live research entry."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import original_scattering_records as experiment
RECEIPT=HERE/'research_round_801_checks.json'
RESULT=HERE/'original_scattering_records_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    actual=experiment.run()
    assert actual==json.loads(RESULT.read_text(encoding='utf-8'))
    assert actual['retained_state_blocks']==json.loads((HERE/'native_scattering_blocks_probe_results.json').read_text(encoding='utf-8'))
    for n in range(776,801):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text(encoding='utf-8'))
        for section in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[section].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_801.md';body=note.read_text(encoding='utf-8')
    assert body.count('$$')==24
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,13)]
    assert not any(ord(c)<32 and c not in '\r\n' for c in body)
    docs=[note,STAGE/'802/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    links=0
    for doc in docs:
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text(encoding='utf-8-sig')):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            path=(doc.parent/target.split('#')[0].strip('<>')).resolve()
            assert path.exists() or (writing and path==RECEIPT.resolve()),(str(doc),target)
            links+=1
    numbers=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if p.parent.name.startswith('archive_') and m:numbers.append(int(m[1]))
    assert sorted(n for n in numbers if n<=801)==list(range(1,802))
    frozen=[STAGE/'800/research_round_800_checks.json',STAGE/'research_note_800.md',
            STAGE/'research_note_799.md',STAGE/'research_note_798.md',STAGE/'research_note_796.md',
            STAGE/'research_note_795.md',STAGE/'research_note_793.md',HERE/'drafts/STATUS.md']
    fresh=[note,HERE/'original_scattering_records.py',RESULT,HERE/'native_scattering_blocks_probe.py',
           HERE/'native_scattering_blocks_probe_results.json',HERE/'drafts/research_note_801_working.md',
           Path(__file__),STAGE/'802/drafts/STATUS.md']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
    return dict(round=801,date='2026-10-05',all_checks_passed=True,
        fresh_test_groups=1,cumulative_numbered_test_groups_from_800=3575,
        saved_result_reproduced=True,historical_manifest_evidence=history,
        previous_776_through_800_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=801,local_links_checked=links,display_equations=12,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='Within the original real compact QME branch, its fixed S-matrix determines an additional outgoing Cauchy record menu. Causal factorization preserves the full old source family and contacts in one joint generating function. Old interior records are retained through the common output-picture dictionary, not identified with the new outgoing records. Internal matrix blocks and the original joint-state matrix functional preserve all source blocks without assuming product preparation. This is a formal input-output connection, not a continuum distinguishing-signal calculation or autonomous terminal apparatus.',
        original_S_used_instead_of_designed_Cayley_dilation=True,
        old_source_contacts_and_mixed_menu_coefficients_retained=True,
        old_interior_and_new_outgoing_records_distinguished=True,
        original_joint_state_and_full_source_blocks_retained=True,
        product_preparation_not_assumed=True,
        continuum_nonzero_readout_signal_proven=False,
        abstract_local_BV_faithfulness_proven=False,
        autonomous_preparation_or_terminal_proven=False,
        finite_coupling_probability_proven=False,
        original_graph_to_continuum_map_proven=False,
        finite_probes_are_continuum_proof=False,
        spacetime_group_and_action_remain_inputs=True,
        original_all_physics_unification_complete=False,
        independent_agent_review=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite a frozen receipt.'
    report=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads(RECEIPT.read_text(encoding='utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert saved[key]==report[key],key
    print(json.dumps({k:v for k,v in report.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
