"""Audit saved evidence and scope only; do not rerun or modify frozen science."""
from pathlib import Path
import hashlib,json,re
BASE=Path(__file__).resolve().parents[1]
ROOT=BASE.parent.parent
HERE=Path(__file__).resolve().parent
REPORT=HERE/'drafts/effective_scope_decision_after_918.md'
OUTPUT=HERE/'drafts/effective_scope_after_918_checks.json'
NAV_BASELINES={'research_cognition_physics\\README.md': {'old_sha256': '9559b3549691993dc86102f06ca66d52ee5b56c45aa7d44c6a0bf1959ce52c94', 'inserted_bytes': 664, 'after_heading_offset': 52}, 'research_cognition_physics\\research_direction.md': {'old_sha256': '0019a19755dce33c5e8049875daa503f852a28589c10abd0a72001ab5690530c', 'inserted_bytes': 664, 'after_heading_offset': 64}, 'research_cognition_physics\\RESEARCH_STATE.md': {'old_sha256': '424f5dc49c4a72ae2345d4f00a665e93446b7937677955020923b6b3cc17dc9e', 'inserted_bytes': 664, 'after_heading_offset': 16}, 'research_cognition_physics\\archive_764_\\README.md': {'old_sha256': 'ff8e8d343f759a0c575d2b34d1b5e47d4bd785cbfd04f65053911dc3989e673b', 'inserted_bytes': 638, 'after_heading_offset': 73}, 'research_cognition_physics\\archive_764_\\文件索引.md': {'old_sha256': 'fb1fac1a3bf2924055d810f57b63ef00e106b367614aca176ea129c5f624065b', 'inserted_bytes': 638, 'after_heading_offset': 16}, 'research_cognition_physics\\archive_764_\\阶段成果总览.md': {'old_sha256': '0708d5680347ff3685dbf40683ce997f288d7b8ab82fba5e791466dfc3e8518a', 'inserted_bytes': 638, 'after_heading_offset': 34}, 'research_cognition_physics\\archive_764_\\跨阶段主题索引.md': {'old_sha256': 'ba9b1337c53d947c99a1481aee100ec37ac13da81ce7ae440fba4f4701c049ca', 'inserted_bytes': 638, 'after_heading_offset': 25}, 'research_cognition_physics\\archive_764_\\_shared\\notes\\unified_physics_condition_ledger_current.md': {'old_sha256': 'd586625da65d0b87596d79e1971183568be0b8dea8b608c205c89938adfab8e8', 'inserted_bytes': 650, 'after_heading_offset': 49}}
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
checked={}
for n in range(915,919):
 receipt=read(BASE/f'{n}/research_round_{n}_checks.json')
 assert receipt['formal_reports']==n and receipt['all_delivery_checks_passed']
 for key in ('frozen_inputs','new_scientific_and_entry_files'):
  for name,expected in receipt[key].items():
   assert sha(ROOT/name)==expected, name
   if name in checked: assert checked[name]==expected,name
   checked[name]=expected
r=read(BASE/'918/euler_density_contacts_results.json')
assert r['same_action_background_response_receiver']
assert r['physical_error_not_removed_by_representation']
assert not r['actual_finite_observable_error_certified']
assert not r['full872_response_computed'] and not r['full_goal_completed']
assert not r['small_internal_driver_resolved_by_coordinate_difference']
assert 5.70e-5 < r['same_Gauss_variation_sample_max'] < 5.71e-5
assert 1.03e-12 < r['corrected_density_internal_driver_max'] < 1.05e-12
assert r['same_Gauss_variation_identity_error']<1e-18
v=read(BASE/'916/covariant_residual_validation_results.json')
assert v['projected_response_on_approximate_background_requires_correction']
assert not v['physical_solution_error_certified']
links=[]
for name in re.findall(r'\]\(([^)]+)\)',REPORT.read_text(encoding='utf-8')):
 if name.startswith(('https:','http:','#')): continue
 target=(REPORT.parent/name.split('#')[0]).resolve()
 assert target.is_file() or target==OUTPUT.resolve(),name
 links.append(name)
nav_hashes={}
for name,old in NAV_BASELINES.items():
 data=(ROOT/name).read_bytes()
 at=old['after_heading_offset']; count=old['inserted_bytes']
 assert hashlib.sha256(data[:at]+data[at+count:]).hexdigest()==old['old_sha256'],name
 added=data[at:at+count].decode('utf-8')
 assert '## 截至918：有限预测与连续外推分开验收' in added
 for target in re.findall(r'\]\(([^)]+)\)',added):
  p=((ROOT/name).parent/target).resolve()
  assert p.is_file() or p==OUTPUT.resolve(),target
 nav_hashes[name]=sha(ROOT/name)
new={'date':'2026-10-06','audit_only':True,'new_scientific_groups':0,
 'formal_reports':918,'cumulative_numbered_test_groups':3703,
 'controlled_effective_description_admissible':True,
 'microscopic_continuity_assumed':False,'physical_minimum_scale_assumed':False,
 'uncut_completion_universal_gate':False,'full_field_arbitrary_precision_required':False,
 'finite_prediction_conservation_feedback_budget_still_required':True,
 'physical_observable_error_certified':False,'stage_completed':False,
 'actual_918_sample_diagnostics_not_error_bounds':{k:r[k] for k in (
  'raw_tensor_internal_driver_times_volume_max','corrected_density_internal_driver_max',
  'same_Gauss_variation_sample_max','corrected_density_diffeomorphism_driver_max',
  'small_internal_driver_resolved_by_coordinate_difference')},
 'historical_files_checked':len(checked),'frozen_inputs_verified':checked,
 'local_report_links_checked':len(links),'navigation_files_checked':len(NAV_BASELINES),
 'preexisting_navigation_content_preserved':True,'navigation_hashes_at_audit':nav_hashes,
 'new_document_and_verifier_hashes':{str(p.relative_to(ROOT)):sha(p) for p in (REPORT,Path(__file__))},
 'all_evidence_checks_passed':True,'scientific_experiments_rerun':False,
 'images_checked':False,'app_goal_changed':False}
OUTPUT.write_text(json.dumps(new,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:new[k] for k in ('all_evidence_checks_passed','historical_files_checked',
 'local_report_links_checked','navigation_files_checked','preexisting_navigation_content_preserved',
 'formal_reports','new_scientific_groups','physical_observable_error_certified')},ensure_ascii=False))
