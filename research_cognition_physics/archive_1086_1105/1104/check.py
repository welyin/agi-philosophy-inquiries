"""Finite conditional checks; reuses round481's n=393, no recurrence search."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCES = [
    "research_cognition_physics/archive_370_428/research_note_372.md",
    "research_cognition_physics/archive_370_428/research_note_373.md",
    "research_cognition_physics/archive_467_530/research_note_481.md",
    "research_cognition_physics/archive_1086_/1088/finite_carrier_proof.md",
    "research_cognition_physics/archive_1086_/1092/proof.md",
    "research_cognition_physics/archive_1086_/1103/proof.md",
]


def td(a, b):
    return float(np.abs(np.linalg.eigvalsh(a-b)).sum()/2)


def residual(a):
    return float(np.max(np.abs(a)))


def compute():
    groups = {}
    # Integer arithmetic for the existing round481 rotation, with a new readout.
    n = 393
    a, b = 1, 0
    for _ in range(n):
        a, b = 3*a-4*b, 4*a+3*b
    q = 5**n
    assert a*a+b*b == q*q
    distance = F(abs(b), q)
    velocity = F(-2*a*b, q*q)
    target = F(3**n-1, 3**n+1)
    gap = abs(target-velocity)
    assert distance < F(1, 400)
    assert gap > F(99, 100)
    step_error, read_error = F(1, 10**6), F(1, 10**4)
    margin = gap-2*n*step_error-2*read_error
    assert margin > F(99, 100)
    U_n = np.array([[float(F(a,q)), -float(F(b,q))],
                    [float(F(b,q)), float(F(a,q))]])
    plus = np.array([1.,1.])/np.sqrt(2)
    rho = np.outer(plus, plus)
    evolved = U_n@rho@U_n.T
    assert abs(td(evolved,rho)-float(distance)) < 1e-12
    assert abs(np.trace(np.diag([1.,-1.])@evolved)-float(velocity)) < 1e-12
    groups["reused_481_rational_word_new_velocity_test"] = {
        "steps": n, "trace_distance_to_initial": float(distance),
        "actual_normalized_velocity": float(velocity),
        "target_normalized_velocity": float(target),
        "ideal_gap": float(gap),
        "per_step_half_diamond_budget": float(step_error),
        "each_endpoint_readout_budget": float(read_error),
        "remaining_margin": float(margin),
        "integer_certificate_checks": 4,
        "recurrence_search_repeated": False,
    }
    # Relative eigenphases remove global phase; pi < 22/7 gives a strict bound.
    Q = 64
    bound = F(1,2)-F(4*22,7*Q)
    assert bound == F(17,56)
    groups["analytic_finite_word_bound"] = {
        "Q": Q, "m": 1, "chi_definition": "atanh(1/2)",
        "normalized_velocity_initial": 0,
        "strict_error_lower_bound": str(bound),
        "word_upper_bounds_by_dimension": {
            str(d): Q**(d-1) for d in (2,3,4)
        },
        "whole_word_not_single_step_tolerance": True,
    }
    # Compact invertibility alone is not isometry.
    r = .5
    z = np.linspace(-1,1,17)
    f = lambda x, s: (x+s)/(1+s*x)
    inv_res = residual(f(f(z,r),-r)-z)
    old_sep = 1.
    new_sep = float(f(.5,r)-f(-.5,r))
    assert inv_res < 1e-12
    assert abs(new_sep-old_sep) > .1
    groups["compact_homeomorphism_is_not_reversible_quantum_isometry"] = {
        "inverse_residual": inv_res,
        "original_separation": old_sep, "transformed_separation": new_sep,
        "claimed_CPTP": False,
    }
    # Finite-window sharp trajectory, with fully coherent payload/reference test.
    N = 8
    labels = np.arange(-N,N+1)
    dim = len(labels)
    shift = np.zeros((dim,dim))
    for j in range(dim):
        shift[(j+1)%dim,j] = 1
    beta = lambda k: float(F(3**int(k)-1,3**int(k)+1)) if k >= 0 else -beta(-k)
    v = np.array([beta(k) for k in labels])
    V = np.diag(v)
    e0 = np.eye(dim)[:,N]
    trajectory_error = 0.
    for k in range(-N,N+1):
        state = np.linalg.matrix_power(shift,int(k))@e0
        trajectory_error = max(trajectory_error,abs(float(state@V@state)-beta(k)))
    assert trajectory_error < 1e-12
    wrapped = np.linalg.matrix_power(shift,N+1)@e0
    wrap_gap = abs(float(wrapped@V@wrapped)-beta(N+1))
    assert wrap_gap > 1.9
    support = [-1,0,1]
    displacement = 2
    embed = np.eye(dim)[:,np.array(support)+N]
    target_embed = np.eye(dim)[:,np.array(support)+N+displacement]
    A = np.kron(np.linalg.matrix_power(shift,displacement)@embed,np.eye(4))
    C = np.kron(target_embed,np.eye(4))
    operator_residual = residual(A-C)
    rng = np.random.default_rng(1104)
    W = rng.normal(size=(12,12))+1j*rng.normal(size=(12,12))
    joint = W@W.conj().T
    joint /= np.trace(joint)
    full_residual = residual(A@joint@A.conj().T-C@joint@C.conj().T)
    assert operator_residual < 1e-12 and full_residual < 1e-12
    groups["finite_window_and_coherent_directory_scope"] = {
        "N": N, "directory_dimension": dim,
        "sharp_trajectory_residual": trajectory_error,
        "next_step_wrap_gap": wrap_gap,
        "operator_identity_residual": operator_residual,
        "payload_dimension": 2, "reference_dimension_calibrated": 2,
        "complete_joint_state_residual": full_residual,
        "general_reference_scope": "all passive R by the operator identity",
        "actual_motion_preparation_certified": False,
        "all_state_nonlinear_velocity_law_claimed": False,
    }
    return {
        "round":1104, "status":"PASS_FIXED_CARRIER_SCOPE",
        "groups":groups,
        "source_sha256":{
            p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES
        },
        "new_universal_axioms":0, "new_scientific_experiments":0,
        "scientific_count_increment":0, "scientific_count_total":3860,
        "fixed_carrier_boost_contract_adopted_in_current_theory":False,
        "full_consensus_conjecture_refuted":False,
        "Lorentz_structure_derived_from_cognition":False,
        "goal_complete":False, "empirical_data":False,
    }


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys() == b.keys()
        for k in a:
            compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a) == len(b)
        for x,y in zip(a,b):
            compare(x,y)
    elif isinstance(a,float):
        assert np.isclose(a,b,atol=1e-12,rtol=1e-12), (a,b)
    else:
        assert a == b, (a,b)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=compute()
    if args.write:
        with (HERE/"results.json").open("x",encoding="utf-8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
    else:
        compare(json.loads((HERE/"results.json").read_text(encoding="utf-8")),result)
    print(json.dumps({"status":result["status"],"groups":len(result["groups"]),
                      "saved_results_match":not args.write}))
