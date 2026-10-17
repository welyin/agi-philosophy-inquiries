"""Round 1263, one author source group. Default replay is read-only.

The result binds the same-round scattering code/result and historical 1262.
Only leading material T0 is transported; no actual Ricci instrument is certified.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json
import numpy as np

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/"1262"/"proof.md"
OLD_SHA="657a7d8ba0521bf105cbf6ec95d1bf806571e482cf993c739164cadba16f7516"
assert hashlib.sha256(OLD.read_bytes()).hexdigest()==OLD_SHA
SOURCES={"check.py":HERE/"check.py","results.json":HERE/"results.json",
         "../1262/proof.md":OLD}
source_hashes={name:hashlib.sha256(path.read_bytes()).hexdigest()
               for name,path in SOURCES.items()}
saved=json.loads((HERE/"results.json").read_text(encoding="utf-8-sig"))
assert saved["status"]=="PASS_AUTHOR_WITH_SCOPE"
eps=F(saved["rational_authoritative_bounds"]["state_error_upper"])
b=F(saved["parameters"]["b"])
k=F(saved["parameters"]["k"])
sigma=F(saved["parameters"]["sigma"])
L=b
assert b==10**100 and k==10**40 and sigma==10**65
assert eps<F(32,10**18)
fineA=F(9,26)/(L*L)
null_trace=F(2)/L
B=32+fineA+null_trace
transport=2*B*eps
reference=48*sigma*sigma/(L*L)+(128*k/b)**2/52+F(3,208)/(sigma*sigma)
total=transport+reference+fineA+null_trace
relative=total/26
conditional_R=F(176,7)*total
assert total<F(21,10**16)
assert relative<F(1,10**16)
assert conditional_R<F(6,10**14)
# Exact arithmetic explaining the independent constants, not new samples-as-proof.
assert 4*F(9,4*26)==F(9,26)
assert F(32)*3/2==48
assert F(3,4)/(2*26)==F(3,208)

# The declared window w=6t(T-t)/T^3 is nonnegative on [0,T] and normalized.
time_integral=6*(F(1,2)-F(1,3))
assert time_integral==1
# A signed L1-normalized weight alone does not fix its ordinary integral.
assert abs(F(-1))==1 and F(-1)!=1

rng=np.random.default_rng(126203)
ratios=[]
for _ in range(400):
    p=rng.normal(size=3)*rng.uniform(.01,20)
    q=rng.normal(size=3)*rng.uniform(.01,20)
    Ep=np.linalg.norm(p); Eq=np.linalg.norm(q)
    s=np.sqrt(Ep*Eq); d=np.linalg.norm(p-q)
    trace=(Ep*Eq-np.dot(p,q))/s
    assert trace>=-1e-12
    assert trace<=d*(1+1e-12)
    ratios.append(float(trace/d))
for p,q in [(np.array([1.,0,0]),np.array([-1.,0,0])),
            (np.zeros(3),np.array([1.,0,0]))]:
    Ep=np.linalg.norm(p); Eq=np.linalg.norm(q); s=np.sqrt(Ep*Eq)
    trace=(Ep*Eq-np.dot(p,q))/s if s else 0.
    d=np.linalg.norm(p-q)
    assert trace<=d*(1+1e-12)
    ratios.append(float(trace/d))
out={
    "round":1263,
    "status":"PASS",
    "source_author_groups":1,
    "formal_acceptance_external_to_this_arithmetic":True,
    "trace_scope":"same H states, declared free/minimal material T0 insertion",
    "b_L":str(b),"k":str(k),"sigma":str(sigma),
    "state_error_exact":str(eps),
    "trace_operator_bound_exact":str(B),
    "dynamic_trace_error_exact":str(transport),
    "reference_trace_to_mass_error_exact":str(reference),
    "total_trace_to_original_mass_exact":str(total),
    "normalized_to_mass_26_error_exact":str(relative),
    "conditional_linear_R_error_exact_pi_upper":str(conditional_R),
    "conditional_geometry_time_weight":{
        "example":"6 t (T-t)/T^3 on [0,T], zero outside",
        "nonnegative":True,
        "ordinary_integral":str(time_integral),
        "L1_norm":"1",
        "signed_weight_center_includes_integral":True,
    },
    "display_only":{
        "state_error":float(eps),"dynamic_trace_error":float(transport),
        "reference_error":float(reference),"total_trace_error":float(total),
        "relative_to_mass26":float(relative),
        "conditional_R_error":float(conditional_R),
        "massless_kernel_samples":len(ratios),"max_ratio_bound_one":max(ratios),
    },
    "source_hashes":source_hashes,
    "claims":{
        "arbitrary_reference_and_original_mass_coherence_retained":True,
        "photon_full_traceless_stress_removed":False,
        "band_required_for_trace_bound":False,
        "actual_Mr_trace_joint_instrument_proved":False,
        "complete_interacting_source_from_same_action_proved":False,
        "independent_physical_metric_proved":False,
        "free_comparison_excluded_by_trace_test":False,
    },
}

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-once",action="store_true",
                        help="Create source_results.json exclusively; never overwrite.")
    args=parser.parse_args()
    target=HERE/"source_results.json"
    if args.write_once:
        with target.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(out,stream,ensure_ascii=False,indent=2,sort_keys=True)
            stream.write("\n")
        print("WRITE_ONCE_PASS_1263_SOURCE: one author source group")
    else:
        expected=json.loads(target.read_text(encoding="utf-8-sig"))
        assert expected==out,"Source result or dependency mismatch; no files changed."
        print("READ_ONLY_MATCH_1263_SOURCE: one author source group; no files changed")
    print(json.dumps(out["display_only"],ensure_ascii=False,indent=2))

