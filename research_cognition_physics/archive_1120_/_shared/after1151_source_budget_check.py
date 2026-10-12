"""Unnumbered 1151 follow-up qualification; finite matrices are not the actual field."""
from pathlib import Path
from fractions import Fraction as Q
import argparse, hashlib, json, math
import numpy as np
HERE=Path(__file__).resolve().parent
def distance(a,b):
    return float(np.abs(np.linalg.eigvalsh((a-b+a.conj().T-b.conj().T)/2)).sum()/2)
def outer(v): return np.outer(v,v.conj())
def budget():
    e=Q(4,39999)+Q(2,10**6)
    rad=1-4*e*e
    a=math.sqrt(float(rad))
    safe=Q(99999997,10**8)
    assert safe*safe<rad
    return dict(name="same_actual_model_budget",passed=True,epsilon=float(e),
      target_branch_distance_lower=a,certified_branch_distance_gt=float(safe),
      with_complete_tester_error_1e_minus3_gt=float(safe-Q(2,1000)),
      stationary_source_complete_error_gt=float(safe/2),
      trivial_representation_effect_error_gt=float(Q(1,2)-e),
      actual_parameters_unchanged=True)
def matrix_check():
    # Only a diagnostic of the source-covariance and dephasing identities.
    # This finite charge spectrum is not the physical K of round 1148.
    n=9
    theta=math.pi/4
    W=np.eye(2*n,dtype=complex)
    for total in range(1,n):
        ids=[total,n+total-1]
        W[np.ix_(ids,ids)]=np.array([[math.cos(theta),-1j*math.sin(theta)],
                                    [-1j*math.sin(theta),math.cos(theta)]])
    charges=np.concatenate([np.arange(n),np.arange(n)+1.])
    assert np.linalg.norm(W.conj().T@W-np.eye(2*n))<1e-13
    assert np.linalg.norm(W*charges[None,:]-charges[:,None]*W)<1e-13
    T=np.diag([1.,-1.])
    Rc=np.diag(np.exp(1j*math.pi*np.arange(n)))
    eta=np.ones(n)/math.sqrt(n)
    sigma=.85*outer(eta)+.15*np.eye(n)/n
    sigma_s=Rc@sigma@Rc.conj().T
    def reduce_c(mat):
        return np.einsum("acbc->ab",mat.reshape(2,n,2,n))
    def channel(rho,source):
        return reduce_c(W@np.kron(rho,source)@W.conj().T)
    rng=np.random.default_rng(115100)
    r=rng.normal(size=2)+1j*rng.normal(size=2); r/=np.linalg.norm(r)
    rho=outer(r)
    lhs=channel(rho,sigma_s)
    rhs=T@channel(T.conj().T@rho@T,sigma)@T.conj().T
    residual=float(np.linalg.norm(lhs-rhs))
    ds=distance(sigma,sigma_s); dout=distance(channel(rho,sigma),lhs)
    assert residual<1e-13 and dout<=ds+1e-13
    # The one-call loop uses the A-flavour |0> witness, no controlled W.
    psi=np.array([1.,0.],complex); rho0=outer(psi)
    X=np.array([[0.,1.],[1.,0.]],complex)
    V=(np.eye(2)-1j*X)/math.sqrt(2)
    target=np.concatenate([V@psi,T@V@T.conj().T@psi])/math.sqrt(2)
    branches=[np.eye(2*n),np.kron(T,np.eye(n))]
    # CT(-s) fixes this particular |0> input. It is still the stated circuit.
    joint_input=np.kron(rho0,sigma)
    candidate=W@joint_input@W.conj().T
    blocks=[[reduce_c(branches[a]@candidate@branches[b].conj().T)/2
             for b in range(2)] for a in range(2)]
    actual=np.block(blocks)
    ideal=outer(target)
    delta=distance(actual,ideal)
    d0=distance(2*blocks[0][0],outer(V@psi))
    d1=distance(2*blocks[1][1],outer(T@V@T.conj().T@psi))
    branch_distance=distance(outer(V@psi),outer(T@V@T.conj().T@psi))
    assert abs(branch_distance-1)<1e-13
    assert d0+d1<=2*delta+1e-13
    assert ds>=branch_distance-d0-d1-1e-13
    assert ds>=1-2*delta-1e-13
    # Invariant source: its marginal output is covariant and cannot meet
    # a small full-loop error, despite allowing a nontrivial source representation.
    stationary=np.eye(n)/n
    cstat=W@np.kron(rho0,stationary)@W.conj().T
    a_stat=np.block([[reduce_c(branches[a]@cstat@branches[b].conj().T)/2
                    for b in range(2)] for a in range(2)])
    delta_stat=distance(a_stat,ideal)
    assert delta_stat>=.5-1e-13
    return dict(name="source_covariance_and_single_input_block_budget",passed=True,
      covariance_residual=residual,source_translate_distance=ds,output_distance=dout,
      complete_tester_error=delta,branch_error_sum=d0+d1,
      branch_distance=branch_distance,stationary_source_error=delta_stat,
      same_actual_field_simulated=False)
def run():
    # Bound to the accepted data; do not recompute old field diagnostics.
    a=json.loads((HERE.parent/"1151/acceptance.json").read_text(encoding="utf-8"))
    source_checks=0
    for name in ["proof.md","results.json"]:
        raw=(HERE.parent/"1151"/name).read_bytes()
        assert hashlib.sha256(raw).hexdigest()==a["core_sha256"][name]
        source_checks+=1
    return dict(status="PASS",formal_round=None,nonformal_check_groups=2,
      frozen_actual_sources=source_checks,groups=[budget(),matrix_check()],
      new_formal_or_scientific_count=0,goal_complete=False)
def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    result=run()
    output=HERE/"after1151_source_budget_results.json"
    data=json.dumps(result,ensure_ascii=False,indent=2)+"\n"
    if args.write_results:
        if output.exists(): assert output.read_text(encoding="utf-8")==data
        else: output.write_text(data,encoding="utf-8")
    else:
        old=json.loads(output.read_text(encoding="utf-8"))
        def compare(a,b):
            if isinstance(a,float): assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12),(a,b)
            elif isinstance(a,dict):
                assert a.keys()==b.keys()
                for k in a: compare(a[k],b[k])
            elif isinstance(a,list):
                assert len(a)==len(b)
                for x,y in zip(a,b): compare(x,y)
            else: assert a==b,(a,b)
        compare(result,old)
    print(json.dumps(dict(status="PASS",nonformal_groups=2,
      stationary_error=result["groups"][1]["stationary_source_error"])))
if __name__=="__main__": main()
