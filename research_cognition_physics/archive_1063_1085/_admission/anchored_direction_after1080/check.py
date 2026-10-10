"""Counted-zero anchored direction bridge: exact budgets and instrument algebra.
Reuse frozen certificates, do not re-count their proofs or claim new physical data.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
BASE=Path(__file__).resolve().parent
STAGE=BASE.parents[1]

def run():
    groups=[]
    def check(name,c):
        assert c,name
        groups.append(name)
    old=json.loads((STAGE/"1077/results.json").read_text(encoding="utf-8"))
    m=F(old["uniform_bloch_per_parameter_lower"])
    r0=F(old["reused_480_control_radius"])
    B=F(30000);L=F(2052);r=r0/2;delta=r/(4*B);a=delta/4
    check("reused actual-source inverse and embedding certificate",
          old["bloch_inverse_infinity_certified_upper"]==2000 and m==F(59,120000))
    q=B*L*r0/4;image=q*r0/4+B*a
    check("entire anchored ball lies strictly inside original source chart",
          q==F(1,16) and image==3*r0/64 and image<r0/4)
    check("all rotations and antipodes fit one common original request budget",
          a==F(1,236390400000000) and 2*a==delta/2<delta)
    # Vector basis: I, E_x, E_o, c_x I, c_o I.
    basis=[]
    for s in (0,1):
        for reset in (0,1):
            for z in (0,1):
                v=[F(0)]*5
                k=(1+s) if not reset else (3+s)
                v[k]=F(1 if z==0 else -1,4)
                if z==1:v[0]=F(1,4)
                basis.append((s,reset,z,v))
    coarse=[sum(v[i] for s,r,z,v in basis if (s^r^z)==0) for i in range(5)]
    total=[sum(v[i] for s,r,z,v in basis) for i in range(5)]
    check("eight fine effects sum to identity and coarse report is physical contrast",
          total==[F(1),F(0),F(0),F(0),F(0)]
          and coarse==[F(1,2),F(1,4),F(-1,4),F(-1,4),F(1,4)])
    # Exact pi/2 rotations, with normalised source direction e_z.
    rz=((0,-1,0),(1,0,0),(0,0,1));rx=((1,0,0),(0,0,-1),(0,1,0))
    def apply(mat,v):return tuple(sum(x*y for x,y in zip(row,v)) for row in mat)
    ez=(0,0,1);ab=apply(rz,apply(rx,ez));ba=apply(rx,apply(rz,ez))
    check("same-anchor order witness and antipodal separation",
          ab==(1,0,0) and ba==(0,-1,0) and a/2>0)
    C=1+38/m
    check("complete source and passive-reference canonical error constant",
          C==F(4560059,59))
    # A finite positive target error, including an explicitly nonzero root residual.
    zeta=a/100;eps=min(r/(16*B),zeta/(2*C));eta_b=zeta*m/38
    bound=eps+F(19)/m*(2*eps+eta_b)
    check("joint root and recurrence budget fits a specified positive error",
          eps>0 and eta_b>0 and bound<=zeta and eps<=r/(16*B))
    report_error=a/100
    check("finite per-effect error retains antipodal margin at this scale",
          a/2-2*report_error==12*a/25>0)
    names=("1077/results.json","1077/proof.md","1079/proof.md",
           "_admission/history_aware_direction_after1080.md")
    return {"formal_rounds_added":0,"status":"passed","assertion_groups":len(groups),"groups":groups,
      "r0":str(r0),"old_control_B":str(B),"control_delta":str(delta),"anchored_radius":str(a),
      "parameter_contraction":str(q),"source_image_bound":str(image),
      "coarse_effect_coefficients_I_Ex_Eo_cxI_coI":list(map(str,coarse)),
      "effect_formula":"F=I/2+(b_x-b_o).sigma/4",
      "uniform_antipodal_probability_sup_gap":str(a/2),
      "order_gap_squared":str(a*a/8),
      "canonical_source_half_trace_error_factor":str(C),
      "finite_error_example":{"zeta":str(zeta),"recurrence_epsilon":str(eps),
        "root_residual_eta_b":str(eta_b),"certified_bound":str(bound)},
      "source_hashes":{n:hashlib.sha256((STAGE/n).read_bytes()).hexdigest() for n in names},
      "scope":{"fixed_anchor_and_routing_are_added_permissions":True,
        "exact_task_shell_realized_under_old_ideal_control_contract":True,
        "arbitrary_actual_Q_history_quotient_identified":False,
        "all_fine_outputs_covariant":False,"microscopic_exactness_required_for_all_cognition":False,
        "spatial_generation_goal_completed":False,"new_accepted_cognitive_axioms":0,"empirical_results":0}}
if __name__=="__main__":print(json.dumps(run(),ensure_ascii=False,indent=2))
