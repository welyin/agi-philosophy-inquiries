"""Positive-spectrum finite-window screen; --write creates results, default read-only."""
from fractions import Fraction as F
from pathlib import Path
import json,argparse,random
HERE=Path(__file__).resolve().parent
def show(x):return {"exact":str(x),"display":float(x)}
def run():
 A=F(1); B=F(4); mu=F(1,100)
 denom=7*A*B+(4*B+3*A)*mu
 delta=(A*B+(4*B-3*A)*mu)/denom
 c=6*(A+mu)*(B+mu)/denom
 x=lambda q:c*q/(q+mu)
 assert delta==F(413,2819)
 assert x(A)==1-delta and F(4,3)*x(B)==1+delta
 # Every grid point is diagnostic; monotonicity in proof establishes entire interval.
 for j in range(101):
  xx=x(A+(B-A)*F(j,100)); yy=F(4,3)*xx
  assert max(abs(xx-1),abs(yy-1))<=delta
 # Finite positive rational mixtures: independent sanity checks of endpoint inequality.
 rng=random.Random(1186); tests=0
 for unused in range(80):
  masses=[mu+F(rng.randrange(0,100),13) for k in range(5)]
  weights=[F(rng.randrange(1,100),17) for k in range(5)]
  scalar=[F(rng.randrange(0,30),41) for k in range(5)]
  xx=lambda q:sum(w*q/(q+m) for w,m in zip(weights,masses))
  yy=lambda q:F(4,3)*xx(q)+2*sum(w*q/(q+m) for w,m in zip(scalar,masses))
  ratio=B*(A+mu)/(A*(B+mu))
  assert xx(B)>=ratio*xx(A)
  # Free normalization chosen optimally for the worst two endpoints of this sample.
  norm=2/(xx(A)+yy(B))
  err=max(abs(norm*xx(A)-1),abs(norm*yy(B)-1))
  assert err>=delta;tests+=1
 tol=F(1,100);fraction=(1-7*tol)/(1-tol)
 assert fraction==F(31,33)
 small=F(1,10**6); single_error=small/(A+small)
 assert single_error<tol
 # Negative scalar removes mismatch; non-universal but PSD Gram also removes it.
 gram=[[F(3,4),F(1)],[F(1),F(4,3)]]
 assert gram[0][0]*gram[1][1]-gram[0][1]**2==0
 assert gram[0][0]>0 and gram[1][1]>0
 assert F(4,3)*gram[0][0]==gram[0][1]==1
 massive=F(7,9);signed_scalar=-massive/6
 assert F(4,3)*massive+2*signed_scalar==massive
 return {"schema":"positive_spectrum_joint_response_screen_v1","date":"2026-10-14",
 "status":"PASS_RESTRICTED_KERNEL_CLASS","point_or_common_average_floor":"1/7",
 "band":{"A":"1","B":"4","mass_squared_gap":"1/100",
         "optimal_error":show(delta),"saturating_weight":show(c),
         "x_A":show(x(A)),"y_B":show(F(4,3)*x(B)),
         "rational_mixture_diagnostics":tests,"saturation_grid_diagnostics":101},
 "one_percent_massless_response_fraction_lower":"31/33",
 "negative_controls":{"single_channel_error":show(single_error),
   "signed_scalar_cancels_trace_excess":True,"nonuniversal_positive_Gram_cancels_mismatch":True},
 "scope":{"Lorentz_linear_positive_universal_kernel_is_input":True,
          "no_independent_contact_or_medium_terms":True,
          "finite_band_optimum_is_analytic":True,"actual_material_instruments_certified":False,
          "arbitrary_R_channel_bound_certified":False,"nonlinear_screening_excluded_by_theorem":False,
          "dark_matter_or_energy_predicted":False,"GR_complete":False}}
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--write",action="store_true");args=p.parse_args()
 result=run();dest=HERE/"results.json"
 if args.write:
  with dest.open("x",encoding="utf-8") as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
  print("SAVED_POSITIVE_SPECTRUM_SCREEN")
 else:
  assert json.loads(dest.read_text(encoding="utf-8"))==result
  print("READ_ONLY_MATCH; restricted finite-band optimum 413/2819")
