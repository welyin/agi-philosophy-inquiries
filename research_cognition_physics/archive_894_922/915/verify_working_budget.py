"""915 delivery checks, without changing the formal research count."""
from pathlib import Path
import json,hashlib,re,ast
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'drafts/working_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
old=read(STAGE/'914/research_round_914_checks.json');frozen={}
for k in ('frozen_inputs','new_scientific_and_entry_files'):frozen.update(old[k])
for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
r=read(HERE/'weighted_source_budget_results.json')
for rel,digest in r['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
assert r['status']=='working' and not r['actual_physical_input_errors_bounded']
assert not r['actual_finite_observable_error_certified'] and not r['full_goal_completed']
for row in r['rows']:
 assert len(row['finite_input_error_amplification'])==8
 assert all(v is None for v in row['unknown_actual_physical_input_errors'].values())
 assert row['independent_formula_error']<1e-20
 assert row['diagnostic_random_error']<=row['diagnostic_bound']
 assert abs(row['diagnostic_aligned_complex_error']-row['diagnostic_bound'])<1e-18
assert r['constant_weight_extra_coefficients_max']==0
files=[HERE/'weighted_source_budget.py',HERE/'weighted_source_budget_results.json',HERE/'drafts/research_note_915_working.md',Path(__file__),STAGE/'914/research_round_914_checks.json']
for p in files:
 if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
doc=HERE/'drafts/research_note_915_working.md';links=0
for link in re.findall(r'\]\(([^)]+)\)',doc.read_text('utf-8')):
 target=(doc.parent/link).resolve();assert target.is_file() or target==TARGET.resolve(),link;links+=1
result=dict(round=915,status='working',formal_reports=914,cumulative_numbered_test_groups=3699,
 all_delivery_checks_passed=True,prior_frozen_entries_verified=len(frozen),local_working_report_links_checked=links,
 numerical_source_decomposition_checked=True,actual_physical_input_errors_bounded=False,
 full_goal_completed=False,files={str(p.relative_to(ROOT)):sha(p) for p in files},visual_checks_performed=False,app_goal_changed=False)
assert not TARGET.exists()
TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='files'},ensure_ascii=False,indent=2))
