"""Round 1143 complete reproducible evidence; default read-only.
Finite CAR and controlled continuous-mode quadrature are not a numerical simulation
of the actual continuum BT process. Its domain/response proof is in proof.md.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,math
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
FROZEN={'1143/proof.md': '5ff71d0adfaf73b5402e8c1e7cf68ccc6adddc9b07fbe06f9155c88bc56601a4', '_shared/round1143_admission.md': '344c9c8e15aac516fd6760e47a99f78d4b7ddf06986f002eb777ba21096f8e40', '_shared/after1142_active_response_preflight.md': 'f31eb3b2d4a4ee5fee5fb03f1eb69d0512058457493ad6d091a464fe68e46576', '_shared/after1142_active_response_preflight.py': '1b92c0968fe573ea52403fe7943f1208a4fe20ccd236851cd966fd7088c23f46', '_shared/after1142_active_response_preflight.json': '74d578dbfc5b3ef217edea1cd0ef5f95d3d57fd7d2809690501c6605a68a80e1', '_shared/after1142_active_actual_check.py': 'cddde07058c6ea2d6d5d3fbca647d8b625960d28507b82f613566ca34e13b7c3', '_shared/after1142_active_actual_results.json': 'a26f1c76b1973cc38dbe563488255d4afe1f14a6ca53b6e225676ba3c1607b34', '_shared/after1142_active_actual_independent_review.md': 'e34f61d02276a259ec0a76f953384035980e54c3b5ed33e9cc77916073f2f283'}
def load(name):
 path=ROOT/'_shared'/name
 spec=importlib.util.spec_from_file_location(path.stem,path)
 mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 return mod
def calculate():
 hashes={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in FROZEN}
 assert hashes==FROZEN
 ideal=load('after1142_active_response_preflight.py')
 actual=load('after1142_active_actual_check.py')
 i=ideal.run();a=actual.calculate()
 ideal.compare(i,json.loads((ROOT/'_shared/after1142_active_response_preflight.json').read_text(encoding='utf-8')))
 assert a==json.loads((ROOT/'_shared/after1142_active_actual_results.json').read_text(encoding='utf-8'))
 return dict(round=1143,status='PASS',mathematical_check_groups=5,groups={
  'frozen_identity':hashes,
  'original_continuous_g_tail':i['tail'],
  'full_CAR_response':i['full_CAR'],
  'derivative_and_graph_arithmetic':a['derivative_and_graph_bounds'],
  'actual_full_output_budget':dict(rows=a['full_output_rows'],conclusions=a['conclusions'])},
  new_axioms=0,scientific_calibration_increment=0,stage='active',
  scope=dict(actual_autonomous_regional_response=True,
   selector_b=0,T=4,original_parameters_changed=False,
   original_Q_and_passive_R_full_outputs=True,
   witness='original rail1 with C=Z, fixed bare n(g1)',
   normal_source_sea_Y_O_Z_and_program_retained=True,
   positive_program_mass=1,relational_weight_radii=[1.1,1.25],
   finite_response_gap_lower=.0335,isometry_error_upper=.0015,
   new_material_permissions=True,physical_continuum_simulation=False,
   full_local_energy_momentum_balance=False,
   all_axiom_counterexample=False,fixed_lab_energy_instrument=False,
   arbitrary_boost_response=False,terminal_detector_autonomous_energy_closed=False))
def compare(a,b,path=()):
 if isinstance(a,dict):
  assert a.keys()==b.keys(),path
  for key in a:compare(a[key],b[key],path+(key,))
 elif isinstance(a,list):
  assert len(a)==len(b),path
  for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+(str(i),))
 elif isinstance(a,float):
  tol=0. if ('actual_full_output_budget' in path or 'derivative_and_graph_arithmetic' in path) else 1e-12
  assert math.isclose(a,b,rel_tol=1e-10,abs_tol=tol),(path,a,b)
 else:assert a==b,(path,a,b)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
 data=calculate();out=HERE/'results.json'
 if args.write_results:
  if out.exists():raise FileExistsError(out)
  out.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 else:compare(data,json.loads(out.read_text(encoding='utf-8')))
 print('1143 PASS: 5 evidence groups; actual full-output error <0.0015 and fixed-bare-effect gap >0.0335 under reviewed finite contract; new axioms/science 0; stage active.')
