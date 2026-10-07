"""920 delivery: same-action joint Legendre defect and source dictionary."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_920_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,920):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'common_legendre_defect_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
    assert result['N']==17 and result['source_points']==512
    assert result['original_material_source_points'] and result['same912_preparation_and_104_variable_flow']
    assert max(result['algebraic_checks'].values())<1e-14
    assert result['actual_momentum_defects']['dp']['max']>2e-6
    assert result['actual_tangent_momentum_defects']['dp']['max']>2e-5
    assert result['action_defect']['minimum']>0.
    assert result['action_defect_smallness_is_not_source_or_observable_bound']
    assert result['mixed_action_is_same_theory_offshell_dictionary_not_new_physical_counterterm']
    for k in ('full_relation_observable_error_certified','quantum_auxiliary_measure_equivalence_proved','full872_response_computed','full_goal_completed'):assert not result[k]
    note=STAGE/'research_note_920.md';text=note.read_text('utf-8')
    assert text.count('$$')==12 and re.findall(r'\\tag\{(\d+)\}',text)==[str(i) for i in range(1,7)]
    docs=[note,STAGE/'921/drafts/STATUS.md']+[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    links=0
    for doc in docs:
        prose=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',prose):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for q in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',q.name)
        if m and q.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=920)==list(range(1,921))
    oldfiles=[STAGE/'919/research_round_919_checks.json',HERE/'drafts/STATUS.md']+[STAGE/f'research_note_{n}.md' for n in (764,874,904,905,919)]
    newfiles=[note,Path(__file__),STAGE/'921/drafts/STATUS.md',HERE/'common_legendre_defect.py',HERE/'common_legendre_defect_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    out=dict(round=920,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=920,
      cumulative_numbered_test_groups_from_919=3705,historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      same_action_defect_and_joint_sources_connected=True,auxiliary_momentum_equation_retained=True,
      original_material_source_points_tested=512,additional_physical_counterterm_introduced=False,
      physical_observable_error_certified=False,quantum_auxiliary_measure_equivalence_proved=False,
      full872_response_computed=False,full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
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
