"""Scope/evidence audit, not a new physical simulation or error certificate."""
from pathlib import Path
import argparse, ast, hashlib, json, re
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'drafts/finite_scope_checkpoint_checks.json'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text('utf-8'))

def run(writing=False):
    frozen={}
    def incorporate(table):
        for name,digest in table.items():
            assert name not in frozen or frozen[name]==digest, name
            frozen[name]=digest
    for n in range(776,911):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            incorporate(old[key])
    for name in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json',
                 '884/drafts/working_checks.json','887/drafts/working_checks.json',
                 '894/drafts/working_checks.json','896/drafts/working_checks.json',
                 '897/drafts/working_checks.json','899/drafts/working_checks.json',
                 '900/drafts/working_checks.json','904/drafts/working_checks.json',
                 '907/drafts/working_checks.json','909/drafts/working_checks.json'):
        incorporate(read(STAGE/name)['files'])
    prior=read(STAGE/'909/drafts/scope_reaudit_checks.json')
    incorporate(prior['preserved_files'])
    incorporate(prior['evidence_and_audit_hashes'])
    for name,digest in frozen.items(): assert sha(ROOT/name)==digest, name
    a=read(STAGE/'909/weak_propagation_validation_results.json')
    b=read(STAGE/'910/physical_direction_validation_results.json')
    c=read(HERE/'future_test_projection_results.json')
    assert a['physical_observer_accuracy_accepted'] is False
    assert a['complete872_quantum_feedback'] is False
    assert b['implementation_checks_passed'] is True
    assert b['actual_continuous_response_error_certified'] is False
    assert b['complete872_response'] is False
    assert b['cumulative_numbered_groups']==3695
    assert c['status']=='working' and [r['N'] for r in c['rows']]==[16,24]
    for key in ('full872_response_computed','actual_net_tail_error_evaluated','full_goal_completed'):
        assert c[key] is False, key
    for r in c['rows']:
        assert r['all_other_projector_components_computed'] is False
        assert r['uniform_reference_support_certified'] is False
    for name,digest in c['source_hashes'].items(): assert sha(STAGE/name)==digest, name
    for name,digest in b['source_hashes'].items(): assert sha(STAGE/'910'/name)==digest, name
    for name,digest in a['source_file_sha256'].items(): assert sha(STAGE/'909'/name)==digest, name
    q=read(STAGE/'909/observer_quadrature_check_results.json')
    probes=[r['relational_probe'] for r in q['rows']]
    assert min(probes)<0<max(probes)
    assert abs(max(probes)-min(probes)-a['differences']['observer_quadrature_range']['relational_probe'])<1e-15
    assert abs(abs(b['conditional_future_read_coefficients'][0]-b['conditional_future_read_coefficients'][1])-b['N16_N24_conditional_read_difference'])<1e-25
    ledger=read(HERE/'drafts/finite_scope_acceptance_910.json')
    assert ledger['formal_rounds']==910 and ledger['new_scientific_groups']==0
    assert ledger['stage_judgement']['controlled_effective_description_admissible'] is True
    for k in ('stage_completed','microscopic_continuity_assumed','physical_minimum_scale_assumed','uncut_completion_required'):
        assert ledger['stage_judgement'][k] is False
    assert not (STAGE/'research_note_911.md').exists(), 'Revisit the checkpoint if formal911 has since completed.'
    note=HERE/'drafts/finite_scope_checkpoint_910.md'
    work=HERE/'drafts/future_test_projection_working.md'
    docs=[note,work]
    docs += [STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    docs += [STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    links=0
    for doc in docs:
        body=doc.read_text('utf-8-sig')
        assert body.count('$$')%2==0, doc
        text=re.sub(r'\$\$.*?\$\$','',body,flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',text):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'): continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()), (doc,link)
            links+=1
    ast.parse(Path(__file__).read_text('utf-8'))
    ast.parse((HERE/'future_test_projection.py').read_text('utf-8'))
    evidence=[note,work,Path(__file__),HERE/'drafts/finite_scope_acceptance_910.json',
              HERE/'future_test_projection.py',HERE/'future_test_projection_results.json']
    evidence += [STAGE/f'research_note_{n}.md' for n in (773,785,799,872,894,898,899,909,910)]
    evidence += [STAGE/x for x in ('899/drafts/effective_scope_joint_audit.md',
        '909/drafts/effective_scope_reaudit_908.md','909/drafts/scope_reaudit_checks.json',
        '909/weak_propagation_validation_results.json','909/observer_quadrature_check_results.json',
        '910/physical_direction_validation_results.json','910/research_round_910_checks.json')]
    hashes={str(x.relative_to(ROOT)):sha(x) for x in evidence}
    frozen_hash=hashlib.sha256(json.dumps(frozen,sort_keys=True,ensure_ascii=False).encode('utf-8')).hexdigest()
    if not writing:
        old=read(TARGET)
        assert old['evidence_and_audit_hashes']==hashes
        assert old['frozen_history_digest']==frozen_hash
    return dict(date='2026-10-06',audit_checks_passed=True,formal_rounds=910,
        cumulative_numbered_groups=3695,new_scientific_groups=0,
        new_physics_certificate=False,scientific_experiments_rerun_by_audit=False,
        round911_status='working',full_goal_completed=False,
        historical_frozen_files_checked=len(frozen),frozen_history_digest=frozen_hash,
        local_links_checked=links,evidence_and_audit_hashes=hashes,
        diagnostic909_probe_range=[min(probes),max(probes)],
        conditional910_mesh_difference=b['N16_N24_conditional_read_difference'],
        working911_sample_p_commutator_maxima=[max(t['max_abs_full_projector_p_commutator'] for t in row['rows']) for row in c['rows']],
        judgement=ledger['stage_judgement'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write: assert not TARGET.exists(), 'Preserve prior checkpoint.'
    result=run(args.write)
    if args.write: TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='evidence_and_audit_hashes'},ensure_ascii=False,indent=2))
