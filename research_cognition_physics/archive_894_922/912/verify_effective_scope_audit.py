"""Incremental scope/evidence audit; no new physics or convergence certificate."""
from pathlib import Path
import argparse, ast, hashlib, json, re
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'drafts/effective_scope_reaudit_checks.json'
NOTE=HERE/'drafts/effective_scope_reaudit_912_working.md'
NAV=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]
NAV += [STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]

def read(p): return json.loads(p.read_text('utf-8-sig'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def audit(writing=False):
    frozen={}
    def incorporate(table):
        for name,digest in table.items():
            assert name not in frozen or frozen[name]==digest, name
            frozen[name]=digest
    for n in range(776,912):
        receipt=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            incorporate(receipt[key])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json',
        '884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json',
        '896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json',
        '900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json',
        '909/drafts/working_checks.json','912/drafts/working_checks.json'):
        incorporate(read(STAGE/path)['files'])
    old=read(STAGE/'909/drafts/scope_reaudit_checks.json')
    incorporate(old['preserved_files']); incorporate(old['evidence_and_audit_hashes'])
    incorporate(read(STAGE/'911/drafts/finite_scope_checkpoint_checks.json')['evidence_and_audit_hashes'])
    for name,digest in frozen.items(): assert sha(ROOT/name)==digest, name
    a=read(HERE/'receiver_initial_completion_results.json')
    b=read(HERE/'receiver_momentum_curl_results.json')
    w=read(HERE/'drafts/working_checks.json')
    f=read(STAGE/'911/research_round_911_checks.json')
    assert f['formal_reports']==911 and f['cumulative_numbered_test_groups_from_910']==3696
    assert f['independent_full_Euler_comparison_passed'] and f['full_four_mode_CAR_identity_reproduced']
    assert not f['actual_finite_observable_error_certified']
    assert w['formal_rounds']==911 and w['current_round_status']=='912_working'
    assert w['preserved_formal911_receipt_sha256']==sha(STAGE/'911/research_round_911_checks.json')
    assert not w['independent_canonical_constraints_numerically_certified']
    assert not w['actual_A_time_evolution_computed']
    assert not a['actual_finite_observable_error_certified'] and not a['full872_numerical_feedback_completed']
    assert b['exact_continuum_total_momentum_zero'] and not b['extra_physical_matter_counterflow_needed_for_this_source']
    for data in (a,b):
        for name,digest in data['source_hashes'].items(): assert sha(STAGE/name)==digest, name
    assert [r['N'] for r in a['rows']]==[17,25]
    assert [r['N'] for r in b['actual_original_source_rows']]==[17,25]
    rows=[]
    for x,y in zip(a['rows'],b['actual_original_source_rows']):
        assert x['total_coordinate_momentum_before']==y['raw_point_sample_total_momentum']
        assert x['conformal_linear_equation_residual']<1e-12
        assert x['independent_canonical_linear_constraints_with_source']['H']>1e-4
        assert x['independent_canonical_linear_constraints_with_source']['M']>1e-4
        assert y['discrete_curl_divergence_max']<1e-13
        assert not y['source_approximation_error_certified']
        rows.append(dict(N=x['N'],conformal_residual=x['conformal_linear_equation_residual'],
            canonical_constraints=x['independent_canonical_linear_constraints_with_source'],
            zero_A_H_defect=x['zero_A_constraint_defect']['H'],
            curl_divergence=y['discrete_curl_divergence_max'],
            source_approximation_difference_not_error_bound=y['pointwise_difference_between_two_source_approximations']))
    links=0
    for p in [NOTE]+NAV:
        body=p.read_text('utf-8-sig')
        assert body.count('$$')%2==0, p
        if p in NAV: assert 'effective_scope_reaudit_912_working.md' in body, p
        body=re.sub(r'\$\$.*?\$\$','',body,flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',body):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'): continue
            target=(p.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()), (p,link)
            links+=1
    ast.parse(Path(__file__).read_text('utf-8'))
    evidence=[NOTE,Path(__file__),HERE/'receiver_initial_completion_results.json',
        HERE/'receiver_momentum_curl_results.json',HERE/'drafts/working_checks.json',
        STAGE/'research_note_911.md',STAGE/'911/research_round_911_checks.json',
        STAGE/'899/drafts/effective_scope_joint_audit.md',STAGE/'911/drafts/finite_scope_checkpoint_910.md']
    hashes={str(p.relative_to(ROOT)):sha(p) for p in evidence}
    digest=hashlib.sha256(json.dumps(frozen,sort_keys=True,ensure_ascii=False).encode('utf-8')).hexdigest()
    judgement=dict(controlled_effective_description_admissible=True,stage_completed=False,
        microscopic_continuity_assumed=False,physical_minimum_scale_assumed=False,
        uncut_completion_required=False,finite_prediction_and_constraint_errors_required=True,
        nonzero_total_response_required=False,full_pointwise_solution_universally_required=False,
        tolerance_binding_complete=False)
    result=dict(date='2026-10-06',audit_checks_passed=True,formal_rounds=911,
        cumulative_numbered_groups=3696,new_scientific_groups=0,round912_status='working',
        new_physics_error_certificate=False,scientific_experiments_rerun=False,app_goal_changed=False,
        historical_frozen_files_checked=len(frozen),frozen_history_digest=digest,
        local_links_checked=links,rows=rows,judgement=judgement,evidence_and_audit_hashes=hashes,
        next_priority='Bind finite prediction budget; transport the same structure-preserving source through compatible initial data and evolution; control remaining911 terms via899.',
        stop_rule='Stop refinement when the predeclared joint observable/constraint/backreaction budget is met; no uncut limit is required.')
    if not writing:
        previous=read(TARGET)
        for key in ('evidence_and_audit_hashes','frozen_history_digest','rows','judgement'):
            assert previous[key]==result[key], key
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write: assert not TARGET.exists(), 'Preserve prior receipt'
    result=audit(args.write)
    if args.write: TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='evidence_and_audit_hashes'},ensure_ascii=False,indent=2))
