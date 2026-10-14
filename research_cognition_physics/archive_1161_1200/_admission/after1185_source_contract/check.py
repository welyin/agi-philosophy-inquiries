"""Read-only default; --write creates results.json once. Exact conditional bounds only."""
from fractions import Fraction as Q
from pathlib import Path
import argparse, json
HERE = Path(__file__).resolve().parent
def value(x):
    return {"exact": str(x), "display": float(x)}
def compute():
    T=Q(22,7); L=Q(300); w=Q(60000); g=Q(1,10**15)
    C=Q(22700475,13958643712)
    A=2*L*T; B=A*A+w*T
    ev=evv=ew=Q(1,10**9); ef=Q(1,10**12)
    d1=T*ev; d2=T*T*evv/2+T*ew
    ec=2*((2*A+d1)*d1+2*d2+ef*((A+d1)**2+2*(B+d2)))
    u=2*g*301*T; v=g*g*60001*T; d=u*u+v
    r=u**3/(1-u)+u*v/((1-u)*(1-v))+v*v/(1-v)
    tail=2*(2*r*(1+u+d)+r*r+2*u*d+d*d)
    bound=(ec+tail/(g*g))/C
    assert 0<u<1 and 0<v<1
    assert bound<Q(8,125)  # .064, rational comparison
    assert ec<C/20 and tail/(g*g)<C/20
    low0=Q(1,10**12)/(g*g); low1=Q(1,10**12)/g
    assert low0>C and low1>C
    sigma=Q(1,100); m=Q(10**9)
    assert 15*7500**2<30000**2
    assert 945*7500**4<1800000000**2
    com=T*(3*Q(30000)/m**2+Q(1800000000)/(8*m**3))
    rms=sigma+T*150/m
    assert com<Q(3,10**12)
    assert rms<Q(10001,10**6)
    return {
        "schema":"conditional_source_coefficient_budget_v1",
        "status":"PASS_CONDITIONAL_ONLY",
        "date":"2026-10-14",
        "coefficient":value(C), "A":value(A), "B":value(B),
        "delta1":value(d1), "delta2":value(d2),
        "coefficient_error":value(ec), "actual_pair_tail_over_g2":value(tail/(g*g)),
        "xi_radius":value(bound), "certified_radius_upper":"8/125",
        "low_order_negative_controls":{"e0_over_g2":value(low0),"e1_over_g":value(low1)},
        "omitted_preparation_counterexample":{"a":"1","zero_V_W":True,"p2":"1","CC_countermodel":False},
        "com_free_only":{"evolution_upper":value(com),"rms_upper":value(rms),
                        "mass":"1000000000","does_not_share_old_W_bound":True},
        "scope":{
            "full_insertion_matching_is_added_input":True,
            "local_sources_to_full_W_matching_proved":False,
            "actual_source_class_nonempty_proved":False,
            "time_windows_feasible_proved":False,
            "low_order_cancellation_is_required":True,
            "actual_Dyson_bounds_are_independent_inputs":True,
            "old_fock_tail_reused_as_actual_certificate":False,
            "physical_xi_selected":False,
            "new_numbered_round":None,
            "new_adopted_axioms":0,"complete_GR_deliverables":0}}
def main():
    p=argparse.ArgumentParser();p.add_argument("--write",action="store_true"); a=p.parse_args()
    result=compute(); dest=HERE/"results.json"
    if a.write:
        with dest.open("x",encoding="utf-8") as f: json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
        print("SAVED_CONDITIONAL_RESULTS")
    else:
        assert json.loads(dest.read_text(encoding="utf-8"))==result
        print("READ_ONLY_MATCH; conditional xi radius < 0.064; actual material existence not certified")
if __name__=="__main__":main()
