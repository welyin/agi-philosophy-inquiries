"""Reproduce 861 and verify publication scope, frozen history and references."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import magnetic_reduced_hamiltonian as experiment
RECEIPT=HERE/'research_round_861_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    assert experiment.run()==json.loads(experiment.TARGET.read_text('utf-8'))
    for n in range(776,861):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_861.md';body=note.read_text('utf-8')
    assert body.count('$$')==30 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,16)]
    docs=[note,STAGE/'862/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=861)==list(range(1,862))
    frozen=[STAGE/'860/research_round_860_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/magnetic_reference_reduction_working.md']
    frozen += [STAGE/f'research_note_{n}.md' for n in (765,766,767,812,852,853,856,857,858,859,860)]
    frozen += [STAGE.parent/'archive_554_584'/f'research_note_{n}.md' for n in (573,574,578,579)]
    frozen += [STAGE.parent/'archive_585_628'/f'research_note_{n}.md' for n in (590,598,623)]
    frozen += [STAGE.parent/'archive_742_763'/'research_note_753.md',STAGE/'859/shared_probe_reference_probe.py']
    fresh=[note,STAGE/'862/drafts/STATUS.md',HERE/'magnetic_reduced_hamiltonian.py',experiment.TARGET,Path(__file__),HERE/'magnetic_coordinate_canonical_probe.py',HERE/'magnetic_coordinate_canonical_probe_results.json']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=861,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_860=3646,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_860_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=861,local_links_checked=links,display_equations=15,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='Original two-derivative zero-fermion Einstein constraints admit a local h,s,p,magnetic-scalar canonical reduction with explicit physical Hamiltonian, volume and source pullbacks. The inverse-volume candidate is connected at first effective-action order, not as a proved full quantum or off-gauge constraint completion.',
        original_color_and_metric_momenta_transported=True,
        original_859_classical_constraint_clock_root_checked=True,
        local_classical_reduction_proved=True,
        full_high_derivative_EFT_Hamiltonian_completed=False,
        reduced_quantization_or_original_graph_equivalence_proved=False,
        residual_579_divergence_cancelled=False,full_goal_completed=False,
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
