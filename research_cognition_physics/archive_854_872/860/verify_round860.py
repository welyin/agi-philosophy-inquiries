"""Reproduce 860 and verify publication scope, frozen history and references."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import compensated_receiver_bridge as experiment
RECEIPT=HERE/'research_round_860_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    assert experiment.run()==json.loads(experiment.TARGET.read_text('utf-8'))
    for n in range(776,860):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_860.md';body=note.read_text('utf-8')
    assert body.count('$$')==28 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,15)]
    docs=[note,STAGE/'861/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    docs += list((HERE/'drafts').glob('*.md'))
    links=0
    for doc in docs:
        prose=re.sub(r'\$\$.*?\$\$', '', doc.read_text('utf-8-sig'), flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            path=(doc.parent/target.split('#')[0].strip('<>')).resolve()
            assert path.exists() or (writing and path==RECEIPT.resolve()),(str(doc),target)
            links+=1
    numbers=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if p.parent.name.startswith('archive_') and m:numbers.append(int(m[1]))
    assert sorted(n for n in numbers if n<=860)==list(range(1,861))
    frozen=[STAGE/'859/research_round_859_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/compensated_shared_probe_working.md']
    frozen += [STAGE/f'research_note_{n}.md' for n in (764,765,766,767,793,794,795,796,801,803,842,844,847,850,851,852,853,859)]
    frozen += [STAGE.parent/'archive_554_584'/f'research_note_{n}.md' for n in (573,574,578,579)]
    frozen += [STAGE.parent/'archive_585_628'/f'research_note_{n}.md' for n in (589,590,598,623)]
    frozen += [STAGE.parent/'archive_742_763'/f'research_note_{n}.md' for n in (752,753)]
    frozen += [STAGE/'765/joint_covariant_gauge_complex.py',STAGE/'853/native_receiver_content.py',STAGE/'859/probe_reference_cone_compatibility_results.json']
    fresh=[note,STAGE/'861/drafts/STATUS.md',HERE/'compensated_receiver_bridge.py',experiment.TARGET,Path(__file__),HERE/'compensated_probe_kernel_probe.py',HERE/'compensated_probe_kernel_probe_results.json']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=860,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_859=3645,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_859_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=860,local_links_checked=links,display_equations=14,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='For sufficiently small nonzero probe preparations, a declared internal baseline action connects shared local references, the full mixed positive free branch, common background linearized cones and a nonzero original six-coupling leading record coefficient. A smooth family of whole bosonic states is not required for that coefficient. Finite-coupling, nonlinear-cone and graph-continuum conclusions remain open.',
        original_mixed_gauge_complex_extended=True,
        classical_constraint_and_reference_background_from_859_retained=True,
        positive_local_formal_branch_and_leading_record_connection='analytic proof in research_note_860.md',
        all_compensation_source_derivatives_retained=True,
        fixed_full_quantum_state_equivalence_claimed=False,
        nonlinear_or_quantum_exact_common_cones_proved=False,
        full_reduced_canonical_dictionary_completed=False,
        residual_579_divergence_cancelled=False,original_graph_to_continuum_proved=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite a frozen receipt.'
    result=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert old[k]==result[k],k
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
