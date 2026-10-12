"""Active F3 actual-process arithmetic qualification. Default read-only.
The analytic domain and Duhamel arguments are in the sibling Markdown and review;
this script does not simulate the continuum BT process or prove those arguments.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
FROZEN={
 '1133/proof.md':'476dc1d1258b819a72e2650684e1096cfa800d7c7c885c77cc1d41a07ebbde8d',
 '1136/proof.md':'360c450f85f40798c6fb9524b00124364b604434d814c770503e0e97e2196462',
 '1137/proof.md':'b665fe6ac4a8f326b5289b93321206666c1c24d5684e7277d5b683f0ebad1beb',
 '_shared/after1142_active_response_preflight.md':'f31eb3b2d4a4ee5fee5fb03f1eb69d0512058457493ad6d091a464fe68e46576',
 '_shared/after1142_active_response_preflight.py':'1b92c0968fe573ea52403fe7943f1208a4fe20ccd236851cd966fd7088c23f46',
 '_shared/after1142_active_response_preflight.json':'74d578dbfc5b3ef217edea1cd0ef5f95d3d57fd7d2809690501c6605a68a80e1'}
def pack(d):
 return {'fractions':{k:str(v) for k,v in d.items()},'decimal':{k:float(v) for k,v in d.items()}}
def calculate():
 hashes={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in FROZEN}
 assert hashes==FROZEN
 # Conservative C3 radial weight graph bounds; no empirical fit.
 width=F(3,20)
 grad=F(22)
 second=sum(F(abs(a)*j*(j-1)) for j,a in [(4,35),(5,-84),(6,70),(7,-20)])/width**2
 lap=second+2*grad/F(11,10)
 h2=1800+2*grad*20+lap
 assert second==224000 and lap<=224040 and h2<230000
 coeff=F(11,40)
 dV=coeff*16; dY=coeff*16*(68+124)
 dO1=coeff*16*264; dO2=coeff*16*96000
 assert dV==F(22,5) and dY<850 and dO1<1200 and dO2<500000
 deriv=pack(dict(delta_V_norm=dV,weighted_mode_gradient=F(42),
   weight_second=second,weight_Laplacian=lap,weighted_mode_H2=h2,
   weighted_mode_h_norm=F(43),kernel_Y_derivative=F(124),
   added_Y_derivative=dY,added_O_first=dO1,added_O_second=dO2,
   Z_first=F(4),Z_second=F(2560)))
 M=I=F(10**80); A=F(10**34); eps=F(1,10**5)
 sig=F(1,10**30); T=F(4); a=F(1,10**9); dc=F(1,10**6)
 lam=F(11,7); kap=F(11,14); C=F(33536000000000000000)
 vall=16*lam*C+F(558,7)+dV
 mstar=3*M+1-vall
 assert mstar>2*M
 CY=F(24780800000000000000000,49)+16316+850
 C1=F(37200); C2=F(10500000); CZ=F(4); CZZ=F(2560)
 ly=CY*(2*T*sig+T*T/(2*M*sig))
 dcp=28*T/(2*M*M*sig*sig)
 ds=3*T/(eps*A)*(lam*C*512+42*8/dc+(lam+kap)*512+coeff*512)
 dr=(T*F(2*10**21)+2*C1*T*T*F(4*10**10)+F(3,2)*C2*T*T+C1*C1*T**3)/(2*I)
 dz=(175*T+10*CZ*T*T+CZZ*T*T/2+CZ*CZ*T**3/3)/(2*I)
 dd=T*C1*a; dm=F(11,8750); capture=42*dc
 assert dz<F(12,10**77)
 rows={}
 for key,b1,b2 in [('original_code',268,28000),('old_far_menu',600,400000)]:
  p1=b1+94*A; p2=b2+188*b1*A+12000*A*A
  a1=1/sig+p1; a2=1/sig**2+2*p1/sig+p2; bmix=p1/sig+p2
  bt=T*a2/(2*mstar)+T*(a2+2*bmix)/(2*M)+CY*a1*T*T/(2*M)
  total=dm+capture+ds+ly+dcp+bt+dr+dd+dz
  causal=2*(ly+dcp+bt)
  assert total<F(3,2000) and causal<F(1681,10**10)
  rows[key]=pack(dict(initial_K2_upper=a2,delta_material=dm,delta_capture=capture,
    delta_source=ds,delta_Y=ly,delta_DCP=dcp,delta_BT=bt,
    delta_rotor=dr,delta_orientation=dd,delta_Z=dz,
    complete_isometry_error=total,far_menu_probability_error=causal,
    margin_below_public_error=F(3,2000)-total))
 E=F(rows['original_code']['fractions']['complete_isometry_error'])
 gain=(F(1,4)-F(121,800))*F(2,5)
 actual=gain-4*E; public=F(67,2000)
 K2=F(1201*10**69)
 cost=T*K2/(2*(2*M)**2)
 assert cost<F(601,10**91) and actual>public
 extra=pack(dict(ideal_probability_gap_lower=gain,
   actual_probability_gap_from_exact_budget=actual,public_actual_probability_gap=public,
   off_dictionary_isometry_error_upper=cost,
   added_controller_force_upper=2*CZ,
   old_far_menu_target=F(1,10**6)))
 return dict(status='PASS_ARITHMETIC_QUALIFICATION',numbered_round=None,
   frozen_sources=hashes,derivative_and_graph_bounds=deriv,
   full_output_rows=rows,conclusions=extra,
   analytic_domain_proof_verified_by_code=False,
   new_scientific_calibrations=0,continuum_field_simulation=False,
   full_local_energy_momentum_balance_proved=False)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
 data=calculate();dest=HERE/'after1142_active_actual_results.json'
 if args.write_results:dest.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 else:assert data==json.loads(dest.read_text(encoding='utf-8'))
 print('PASS_ARITHMETIC_QUALIFICATION: frozen inputs, rational domain/derivative constants and actual-process budgets; analytic proof reviewed separately; stage active.')
