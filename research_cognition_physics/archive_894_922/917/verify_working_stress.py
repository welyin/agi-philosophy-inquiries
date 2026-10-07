"""Preserve916 and verify the executed917 working result, without closing917."""
from pathlib import Path
import json,hashlib,re,ast
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'drafts/working_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
old=read(STAGE/'916/research_round_916_checks.json');frozen={}
for k in ('frozen_inputs','new_scientific_and_entry_files'):frozen.update(old[k])
for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
r=read(HERE/'receiver_stress_derivative_results.json')
for rel,digest in r['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
assert r['status']=='working' and r['samples']==8
assert r['original_receiver_and_full_offshell_stress_retained']
assert r['independent_original_stress_error']<1e-12 and r['original_frame_error']<1e-12
assert 8e-6<r['actual_stress_divergence_sample_max']<9e-6
assert not r['uniform_source_divergence_bound_certified'] and not r['actual_finite_observable_error_certified'] and not r['full_goal_completed']
assert max(r['independent_central_difference_diagnostics'][-1]['coordinate_derivative_differences'])<1e-9
files=[HERE/'receiver_stress_derivative.py',HERE/'receiver_stress_derivative_results.json',HERE/'drafts/receiver_stress_working.md',Path(__file__),STAGE/'916/research_round_916_checks.json']
for p in files:
 if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
doc=HERE/'drafts/receiver_stress_working.md';prose=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8'),flags=re.S);links=0
for link in re.findall(r'\]\(([^)]+)\)',prose):
 target=(doc.parent/link).resolve();assert target.exists() or target==TARGET.resolve(),link;links+=1
result=dict(round=917,status='working',formal_reports=916,cumulative_numbered_test_groups=3701,all_delivery_checks_passed=True,
 prior_frozen_entries_verified=len(frozen),working_report_links_checked=links,
 full_joint_constraint_driver_computed=False,actual_finite_observable_error_certified=False,full_goal_completed=False,
 files={str(p.relative_to(ROOT)):sha(p) for p in files},visual_checks_performed=False,app_goal_changed=False)
assert not TARGET.exists();TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='files'},ensure_ascii=False,indent=2))
