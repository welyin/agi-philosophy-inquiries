"""1152: control-source necessity for the fixed actual task; not a field simulator."""
from pathlib import Path
import argparse, hashlib, json, math, runpy
import numpy as np
HERE=Path(__file__).resolve().parent
QUAL=runpy.run_path(str(HERE.parent/"_shared/after1151_source_budget_check.py"))
distance=QUAL["distance"]
outer=QUAL["outer"]
def unitary(H,t=1):
    e,v=np.linalg.eigh(H); return (v*np.exp(-1j*t*e))@v.conj().T
def loop(W,G,rho):
    u=[W,G@W@G.conj().T]
    blocks=[[u[a]@rho@u[b].conj().T/2 for b in range(2)] for a in range(2)]
    plus=(blocks[0][0]+blocks[0][1]+blocks[1][0]+blocks[1][1])/2
    return float(np.trace(plus).real),blocks
def boundaries():
    T=np.diag([1.,-1.])
    charge=np.diag([0.,1.])
    swap=np.eye(4,dtype=complex)[[0,2,1,3]]
    Q=np.kron(charge,np.eye(2))+np.kron(np.eye(2),charge)
    assert np.linalg.norm(swap@Q-Q@swap)<1e-14
    sig=np.eye(2)/2; rhoS=outer(np.array([1.,0.]))
    rho=np.kron(rhoS,sig)
    p_relative,blocks=loop(swap,np.kron(T,np.eye(2)),rho)
    p_joint,_=loop(swap,np.kron(T,T),rho)
    def redE(m): return np.einsum("aibi->ab",m.reshape(2,2,2,2))
    actual=np.block([[redE(m) for m in row] for row in blocks])
    X=np.array([[0.,1.],[1.,0.]])
    V=(np.eye(2)-1j*X)/math.sqrt(2); psi=np.array([1.,0.])
    target=np.concatenate([V@psi,T@V@T@psi])/math.sqrt(2)
    full_error=distance(actual,outer(target))
    marginal=redE(swap@rho@swap.conj().T)
    marginal_error=distance(marginal,outer(V@psi))
    assert abs(p_relative-.5)<1e-14 and abs(p_joint-1)<1e-14
    assert abs(marginal_error-.5)<1e-14 and full_error>=.5-1e-14
    assert distance(sig,T@sig@T)==0
    # Correlations can supply relational information even with invariant marginals.
    H=np.zeros((4,4)); H[1,2]=H[2,1]=1
    W=unitary(H,math.pi/4)
    xi=np.array([0.,1.,1j,0.])/math.sqrt(2)
    rho_corr=outer(xi)
    pcorr,_=loop(W,np.kron(T,np.eye(2)),rho_corr)
    after_E=np.einsum("aiaj->ij",rho_corr.reshape(2,2,2,2))
    joint_shift=np.kron(np.eye(2),T)
    relational_distance=distance(rho_corr,joint_shift@rho_corr@joint_shift)
    assert np.linalg.norm(W@Q-Q@W)<1e-13
    assert abs(pcorr-.5)<1e-13
    assert distance(after_E,T@after_E@T)<1e-13 and abs(relational_distance-1)<1e-13
    # A truly trivial source representation keeps the whole loop trivial,
    # including for correlated S/E input; degeneracy permits nontrivial W.
    rng=np.random.default_rng(1152)
    W0=np.zeros((8,8),complex)
    for ids in [np.arange(4),np.arange(4,8)]:
        m=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        W0[np.ix_(ids,ids)]=unitary((m+m.conj().T)/2)
    Ts=np.diag([1.,1.,np.exp(.7j),np.exp(.7j)])
    G0=np.kron(Ts,np.eye(2))
    eta=rng.normal(size=8)+1j*rng.normal(size=8); eta/=np.linalg.norm(eta)
    trivial_p,_=loop(W0,G0,outer(eta))
    residual=float(np.linalg.norm(G0@W0@G0.conj().T-W0))
    assert residual<1e-13 and abs(trivial_p-1)<1e-13
    return dict(name="probability_only_relative_joint_and_correlated_boundaries",passed=True,
      invariant_source_probability_matches_target=p_relative,
      complete_error_despite_probability_match=full_error,
      system_marginal_error=marginal_error,joint_translation_probability=p_joint,
      correlated_input_relative_probability=pcorr,
      invariant_source_marginal_distance=distance(after_E,T@after_E@T),
      correlated_joint_relational_distance=relational_distance,
      trivial_representation_correlated_probability=trivial_p,
      trivial_representation_operator_residual=residual,
      finite_charge_diagnostic_not_actual_continuous_field=True)
def close(a,b):
    if isinstance(a,float): assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12),(a,b)
    elif isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a: close(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): close(x,y)
    else: assert a==b,(a,b)
def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    sources=json.loads((HERE/"source_hashes.json").read_text(encoding="utf-8"))
    for name,digest in sources.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
    result=dict(round=1152,all_passed=True,numerical_check_groups=3,
      groups=[QUAL["budget"](),QUAL["matrix_check"](),boundaries()],
      frozen_source_count=len(sources),qualification_groups_reused_without_double_count=2,
      new_adopted_axioms=0,scientific_calibration_increment=0,
      necessary_not_sufficient_source_condition=True,full_local_parent_constructed=False,
      common_axiom_counterexample_claimed=False,goal_complete=False)
    dest=HERE/"results.json"
    if args.write_results:
        if dest.exists(): close(result,json.loads(dest.read_text(encoding="utf-8")))
        else: dest.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    else: close(result,json.loads(dest.read_text(encoding="utf-8")))
    print(json.dumps(dict(status="PASS",groups=3,
      branch_distance_lower=result["groups"][0]["target_branch_distance_lower"],
      source_distance_required_gt=result["groups"][0]["with_complete_tester_error_1e_minus3_gt"])))
if __name__=="__main__": main()
