"""Nonnumbered finite-block capacity and action-scale admission; read-only default."""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json
HERE=Path(__file__).resolve().parent
S=HERE.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
old=json.loads((S/'1249/quantum_bounds_results.json').read_text('utf-8'))
ledger=old['group_constructive_finite_resource_bound']['ledger']
N=int(ledger['integer_n_power_exponent_decimal'])
EF=int(ledger['Fourier_transport_exponent_decimal'])
EH=int(ledger['Fourier_H_exponent'])
assert 3*N-EF>42 and 6*N-EH-EF>42
ell=F(1,100); V0=ell**3
sqrt3_upper=F(7,4)
assert sqrt3_upper**2>3
dmax_upper=sqrt3_upper*ell/2
roundtrip_upper=6*dmax_upper
assert roundtrip_upper==F(21,400)<F(1,8)<F(1,4)
reserve_slack=F(1,8)-roundtrip_upper
hold_slack=F(1,4)-roundtrip_upper
assert reserve_slack==F(29,400) and hold_slack==F(79,400)
# s=1000*m*2^(2N), m>=max(1,hbar); V0*s^2=m^2*2^(4N).
assert V0*1000**2==1
lin_margin=4*N-EF
quad_margin=8*N-EH-EF
assert min(lin_margin,quad_margin)>42
assert F(2,2**42)<F(1,2**40)
rho0=F(111,100)
inventory_coeff=rho0*V0
assert inventory_coeff==F(111,10**8)
sources={
 'self':HERE/'check.py','proof':HERE/'proof.md',
 'quantum_proof':S/'1249/quantum_bridge.md',
 'quantum_results':S/'1249/quantum_bounds_results.json',
 'hold_proof':S/'1248/proof.md',
 'complete_action_similarity':HERE.parent/'full_action_similarity_screen_after1252.md',
 'locality_scope':HERE.parent/'record_locality_screen_after1252.md',
}
out={
 'status':'PASS_NONNUMBERED_CAUSAL_CAPACITY_AND_ACTION_SCALE',
 'new_scientific_groups':0,
 'rational':{'ell':str(ell),'V0':str(V0),'sqrt3_upper':str(sqrt3_upper),
  'distance_to_fixed_port_upper':str(dmax_upper),'geometric_roundtrip_time_upper':str(roundtrip_upper),
  'reserve_slack':str(reserve_slack),'hold_slack':str(hold_slack),
  'initial_complete_classical_energy_per_s':str(inventory_coeff)},
 'finite_resource':{'s':'1000*m*2^(2N), m=ceil(max(1,hbar))',
  'N_decimal':str(N),'epsilon_upper':'2^(-4N)',
  'linear_error_exponent_margin':str(lin_margin),'quadratic_error_exponent_margin':str(quad_margin),
  'same_bounded_menu_error':'<2^-41<2^-40',
  'physical_time_and_period':'t_s=s*t, period=s/100; real material parameters also scaled'},
 'not_claimed':['actual quantum access or instrument','full-field invariant quantum sector',
  'original per-label K numerical precision','added physical source certificate','GR completion'],
 'sha256':{k:sha(v) for k,v in sources.items()}}
ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
path=HERE/'results.json'
if args.write:
 with path.open('x',encoding='utf-8',newline='\n') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
else:
 assert json.loads(path.read_text('utf-8'))==out,'read-only result mismatch'
print(json.dumps({'status':out['status'],'mode':'FIRST_WRITE' if args.write else 'READ_ONLY_MATCH',
 'new_scientific_groups':0,'geometric_roundtrip_upper':str(roundtrip_upper),
 'reserve_slack':str(reserve_slack),'quantum_error':out['finite_resource']['same_bounded_menu_error']}))
