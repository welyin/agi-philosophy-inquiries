"""Read-only evidence audit; writes only its new local receipt, no physics rerun."""
from pathlib import Path
import hashlib,json,re
BASE=Path(__file__).resolve().parents[1]
ROOT=BASE.parent.parent
HERE=Path(__file__).resolve().parent
REPORT=HERE/'drafts/finite_scope_decision_after_913.md'
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
receipt=read(BASE/'913/research_round_913_checks.json')
assert receipt['formal_reports']==913 and receipt['all_delivery_checks_passed']
checked={}
for field in ('frozen_inputs','new_scientific_and_entry_files'):
 for relative,expected in receipt[field].items():
  path=ROOT/relative
  assert sha(path)==expected,relative
  checked[relative]=expected
f=read(BASE/'912/receiver_forced_response_results.json')
a=read(BASE/'913/discrete_record_tangent_results.json')
b=read(BASE/'913/receiver_response_readout_results.json')
row25=next(x for x in f['rows'] if x['N']==25 and x['Tend']>0)
assert not f['actual_finite_observable_error_certified']
assert not receipt['actual_finite_observable_error_certified']
assert not receipt['full_goal_completed']
paired={x['N']:x['full_pairing']['total'][2] for x in b['rows']}
path_gaps={str(x['segments']):x['direct_minus_curvature_source'][2] for x in a['rows']}
links=[]
for target in re.findall(r'\]\(([^)]+)\)',REPORT.read_text(encoding='utf-8')):
 if target.startswith(('https:','http:','#')): continue
 p=(REPORT.parent/target.split('#')[0]).resolve()
 assert p.is_file() or p == (HERE/"drafts/finite_scope_decision_checks.json").resolve(),target
 links.append(target)
nav=[BASE.parent/x for x in ('README.md','research_direction.md','RESEARCH_STATE.md')]
nav += [BASE/x for x in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
for p in nav: assert '## 截至913：有限共同预测的验收裁定' in p.read_text(encoding='utf-8-sig')
summary={
 'date':'2026-10-06','audit_only':True,'new_scientific_groups':0,'formal_reports':913,
 'cumulative_numbered_test_groups':3698,'all_evidence_checks_passed':True,
 'controlled_effective_description_admissible':True,'actual_joint_budget_certified':False,
 'stage_completed':False,'microscopic_continuity_assumed':False,'minimum_physical_scale_assumed':False,
 'uncut_completion_universal_gate':False,'full_field_reconstruction_universal_gate':False,
 'finite_conservation_and_feedback_errors_still_required':True,
 'source25_final_constraint_diagnostics':row25['final_total_constraints'],
 'path_derivative_differences_diagnostic_not_error_bounds':path_gaps,
 'grid_pairing_difference_diagnostic_not_error_bound':abs(paired[17]-paired[25]),
 'historical_receipt_entries_checked':len(checked),'frozen_inputs_verified':checked,
 'local_report_links_checked':len(links),'navigation_entries_verified':len(nav),
 'evidence_hashes':{str(p.relative_to(ROOT)):sha(p) for p in (BASE/'912/receiver_forced_response_results.json',BASE/'913/discrete_record_tangent_results.json',BASE/'913/receiver_response_readout_results.json',BASE/'913/research_round_913_checks.json')},
 'new_document_and_verifier_hashes':{str(p.relative_to(ROOT)):sha(p) for p in (REPORT,Path(__file__))},
 'scientific_experiments_rerun':False,'images_checked':False,'app_goal_changed':False,
}
(HERE/'drafts/finite_scope_decision_checks.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in summary.items() if k not in ('frozen_inputs_verified','evidence_hashes','new_document_and_verifier_hashes')},ensure_ascii=False,indent=2))
