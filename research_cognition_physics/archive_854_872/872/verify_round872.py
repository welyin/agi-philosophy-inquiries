"""Reproduce872, audit history and freeze the new source construction."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import read_vertex_source_closure as experiment
RECEIPT=HERE/'research_round_872_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    science=experiment.run()
    assert science==json.loads(experiment.TARGET.read_text('utf-8'))
    assert science['cumulative_numbered_groups']==3657
    assert science['fresh_numbered_groups']==1
    correction=STAGE/'871/research_round_871_count_correction.json'
    fix=json.loads(correction.read_text('utf-8'))
    prior=STAGE/'871/research_round_871_checks.json'
    assert sha(prior)==fix['preserved_receipt_sha256']
    assert json.loads(prior.read_text('utf-8'))[fix['corrected_field']]==fix['recorded_value']
    oldscience=json.loads((STAGE/'871/receiver_history_influence_results.json').read_text('utf-8'))
    assert fix['correct_value']==fix['previous_count']+fix['fresh_groups']==oldscience['cumulative_numbered_groups']==3656
    assert fix['correct_value']+science['fresh_numbered_groups']==science['cumulative_numbered_groups']
    for n in range(776,872):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_872.md';body=note.read_text('utf-8')
    assert body.count('$$')==16 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,9)]
    for formula in re.findall(r'\$\$(.*?)\$\$',body,re.S):
        assert not re.search(r'(?<!\\)\b(quad|qquad|left|right)\b',formula),formula
    docs=[note,STAGE/'873/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    docs += list((HERE/'drafts').glob('*.md'))
    links=0
    for doc in docs:
        prose=re.sub(r'\$\$.*?\$\$', '', doc.read_text('utf-8-sig'), flags=re.S)
        prose=re.sub(r'`[^`\n]*`', '', prose)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            path=(doc.parent/target.split('#')[0].strip('<>')).resolve()
            assert path.exists() or (writing and path==RECEIPT.resolve()),(str(doc),target)
            links+=1
    numbers=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if p.parent.name.startswith('archive_') and m:numbers.append(int(m[1]))
    assert sorted(n for n in numbers if n<=872)==list(range(1,873))
    frozen=[prior,correction,HERE/'drafts/STATUS.md']
    frozen += [STAGE/f'research_note_{n}.md' for n in (785,791,802,803,810,845,846,859,860,869,870,871)]
    frozen += [STAGE/'846/constrained_six_leg_recursion.py',STAGE/'859/probe_reference_cone_compatibility_results.json']
    fresh=[note,STAGE/'873/drafts/STATUS.md',HERE/'read_vertex_source_closure.py',experiment.TARGET,Path(__file__),HERE/'drafts/read_vertex_total_source_working.md']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=872,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_871=science['cumulative_numbered_groups'],saved_result_reproduced=True,
        previous871_count_correction_verified=True,
        historical_manifest_evidence=history,previous_776_through_871_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=872,local_links_checked=links,display_equations=8,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='The new869 read action fixes its density Hessian and Dirac mass variation. Together with the old full recursion they construct the lambda_alpha^2 hbar^2 four-leg projected source and relative mean under the inherited regular local causal/Cauchy/Noether contract. Exact nilpotent Euler and original-metric derivative checks calibrate the formula; no original continuous source value or nonzero mean geometry is inferred.',
        extra_independent_second_read_vertex=False,
        full_projected_source_and_relative_mean_defined=True,
        original_continuous_total_source_numerically_evaluated=False,
        original_total_mean_difference_proved_nonzero=False,
        complete_graph_continuum_dynamics_matching_proved=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite frozen receipt.'
    result=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope','cumulative_numbered_test_groups_from_871'):assert old[k]==result[k],k
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
