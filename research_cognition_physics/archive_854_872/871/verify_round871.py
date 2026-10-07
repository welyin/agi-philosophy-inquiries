"""Reproduce871 and verify frozen research, numbering and links."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import receiver_history_influence as experiment
RECEIPT=HERE/'research_round_871_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def verify(writing=False):
    assert experiment.run()==json.loads(experiment.TARGET.read_text('utf-8'))
    for n in range(776,871):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_871.md';body=note.read_text('utf-8')
    assert body.count('$$')==18 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,10)]
    for formula in re.findall(r'\$\$(.*?)\$\$',body,re.S):
        assert not re.search(r'(?<!\\)\b(quad|qquad|left|right)\b',formula),formula
    docs=[note,STAGE/'872/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=871)==list(range(1,872))
    frozen=[STAGE/'870/research_round_870_checks.json',HERE/'drafts/STATUS.md']
    frozen += [STAGE/f'research_note_{n}.md' for n in (785,796,802,803,810,844,845,846,853,860,863,866,867,868,869,870)]
    frozen += [STAGE/'870/receiver_four_point_backreaction.py',STAGE/'870/receiver_four_point_backreaction_results.json',STAGE.parent/'archive_629_652/research_note_643.md',STAGE.parent/'archive_585_628/research_note_622.md']
    fresh=[note,STAGE/'872/drafts/STATUS.md',HERE/'receiver_history_influence.py',experiment.TARGET,Path(__file__),HERE/'drafts/total_source_working.md',HERE/'drafts/conditional_receiver_influence_working.md']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=871,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_870=3655,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_870_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=871,local_links_checked=links,display_equations=9,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='The original870 receiver preparation is exactly a positive mixture of four quasifree receiver inputs, evaluated with the same formal S, records and sources. Its full finite-CAR conditional history weight retains exterior-mode excursions and has an explicit nonGaussian correction with nonsingular cofactor source insertions. Equal classical histories cannot distinguish r through quadratic conditional means; original870 quantum variance remains different. No continuous finite-coupling determinant or evaluated full mean-source sum is claimed.',
        exact_positive_receiver_input_branches=4,
        full_finite_CAR_history_and_source_identity=True,
        exterior_mode_excursions_retained=True,
        same_fixed_formal_S_records_and_sources=True,
        arbitrary_classical_history_quadratic_means_equal=True,
        original870_nonzero_noise_coefficient_recovered=True,
        continuous_finite_coupling_determinant_proved=False,
        original_four_leg_total_source_evaluated=False,
        complete_graph_continuum_dynamics_matching_proved=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite frozen receipt.'
    result=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert old[k]==result[k],k
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
