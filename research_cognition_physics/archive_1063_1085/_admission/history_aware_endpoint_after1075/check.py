"""Read-only rational checks of a mature 481 corollary; not a new spatial model."""
from fractions import Fraction as F
import json
from pathlib import Path
import hashlib

def run():
    B,L=30000,2052
    r0=F(1,246240000); r=r0/2
    eps=r/(16*B); delta=r/(4*B)
    assert B*L*r0==F(1,4)
    assert eps==F(1,236390400000000)
    assert delta==F(1,59097600000000)
    assert 4*(9+9+1)==76
    lam=F(1,4)+76*B*eps
    assert lam==F(1,4)+19*r/4 < F(25000001,100000000)
    assert B*delta+r/4+2*B*eps==5*r/8<r
    # Finite iteration/error budget: only the proven scalar upper bounds.
    budget=delta/16
    k=0
    while 38*r*lam**k>budget/4:k+=1
    nu=budget*(1-lam)/(4*38)
    assert 38*(r*lam**k+nu/(1-lam))<=budget/2
    eta=F(1,1000)
    dimension_checks=[]
    for d in (1,2,3,4,7):
        x=tuple(F(i,10*d) for i in range(d))
        a=tuple(F((-1)**i,100*d) for i in range(d))
        b=tuple(F(i+1,1000*d) for i in range(d))
        add=lambda u,v: tuple(s+t for s,t in zip(u,v))
        assert add(add(x,a),b)==add(x,add(a,b))
        assert add(add(x,a),tuple(-u for u in a))==x
        # Saturating arithmetic errors; these are not quantum trajectories.
        e=tuple(eta for _ in range(d))
        ideal=add(x,add(a,b))
        two=add(add(add(add(x,a),e),b),e)
        direct=add(ideal,tuple(-v for v in e))
        assert max(abs(s-t) for s,t in zip(two,direct))==3*eta
        dimension_checks.append(d)
    root=Path(__file__).resolve().parents[4]
    assets=[
        "research_cognition_physics/archive_467_530/research_note_481.md",
        "research_cognition_physics/archive_467_530/481/positive_time_coordinate_reaccess.py",
        "research_cognition_physics/archive_467_530/481/positive_time_coordinate_reaccess_results.json",
    ]
    return {
        "status":"passed","scientific_round_increment":0,"latest_formal_round":1075,
        "groups":3,"analytic_proof":"proof.md",
        "r":str(r),"epsilon":str(eps),"delta":str(delta),
        "contraction_lambda":str(lam),"contraction_float":float(lam),
        "fixed_branch_target_lipschitz_upper":str(F(B)/(1-lam)),
        "scalar_budget_iteration_count":k,"scalar_endpoint_error_budget":str(budget),
        "arithmetic_translation_dimensions_checked":dimension_checks,
        "two_step_vs_direct_error_factor":3,
        "historical_assets_sha256":{p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in assets},
        "scope":{
            "historical_six_tree_recomputed_separately":True,
            "scalar_checks_are_not_execution_of_recurrence":True,
            "no_new_control_permission":True,
            "no_new_universal_cognitive_axiom":True,
            "actual_spatial_three_dimension_derived":False,
            "empirical_results":0
        }
    }

if __name__=="__main__":
    print(json.dumps(run(),ensure_ascii=False,indent=2))
