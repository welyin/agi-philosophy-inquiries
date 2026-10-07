"""925 delivery: same autonomous record transfer and resource response, scoped spectrum selection."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_925_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,925):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'record_resource_exchange_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
    assert result['all_scientific_checks_passed']
    for row in result['exact_classification']:
        assert row['exact_rational_rank']==row['reference_levels']+row['payload_dimension']-2
        assert row['free_parameters']==2 and row['affine_spectrum_and_additive_endpoint_gap_span_nullspace']
    assert result['model_dimension']==60 and result['positive_total_H_minimum']>0
    assert max(result['conservation_residuals'].values())<1e-13
    for key in ('same_transfer_current_identity_error','all_input_isometry_error','relational_recovery_error','nonlinear_probability_formula_max_error'):assert result[key]<1e-12
    assert abs(result['bare_receiver_entanglement_fidelity']-1/3)<1e-13
    assert abs(result['bare_receiver_trace_distance_from_ideal']-2/3)<1e-13
    flow=result['resource_energy_flow']
    assert abs(flow['reference_energy_change']-.7)<2e-13 and abs(flow['endpoint_energy_change']+.7)<2e-13
    assert abs(sum(flow.values()))<2e-13
    assert abs(result['nonlinear_spectrum_counterexample']['maximum_transfer_probability']-25/61)<1e-13
    assert abs(result['fixed_dictionary_without_shift_ledger_violation']-.8)<1e-13
    for key in ('total_energy_conservation_alone_forces_resonance','gauge_group_or_angular_reference_derived','universal_mass_law_derived','selected_three_dimensional_spacetime','GR_generated','continuum_or_unlimited_reference_required','whole_stage_completed','full_goal_completed'):assert not result[key]
    note=STAGE/'research_note_925.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==18 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,10)]
    newdocs=[note,STAGE/'926/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    links=0
    for doc in newdocs+nav:
        text=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',text):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for q in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',q.name)
        if m and q.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=925)==list(range(1,926))
    assert '001—925轮共925份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—925的695份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:assert '925：同一交换连接未知记录' in q.read_text('utf-8-sig')
    oldfiles=[STAGE/'924/research_round_924_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'record_resource_exchange.py',HERE/'record_resource_exchange_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    out=dict(round=925,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=925,
      cumulative_numbered_test_groups_from_924=3710,historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      same_record_resource_operator=True,spectrum_selection_only_in_declared_architecture=True,
      exact_rational_spectral_classification=True,full_unknown_input_isometry_checked=True,
      finite_detuning_counterexample=True,total_energy_alone_not_mistaken_for_bare_energy_conservation=True,
      known_battery_channel_and_programming_results_reused=True,battery_size_not_optimized=True,
      physical_gauge_group_or_GR_generated=False,full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
      frozen_inputs={str(q.relative_to(ROOT)):sha(q) for q in oldfiles},new_scientific_and_entry_files={str(q.relative_to(ROOT)):sha(q) for q in newfiles})
    if not writing:
        before=read(TARGET)
        for k in ('frozen_inputs','new_scientific_and_entry_files'):assert before[k]==out[k],k
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    d=run(a.write)
    if a.write:TARGET.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in d.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
