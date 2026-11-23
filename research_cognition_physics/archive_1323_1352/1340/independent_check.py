"""Round 1340 separate finite-dimensional source, drift and reference check."""
from pathlib import Path
from fractions import Fraction
import argparse
import json
import math
import numpy as np

def evolve(H,t):
    val,U=np.linalg.eigh(H)
    return (U*np.exp(-1j*t*val))@U.conj().T

def small_parent(directions,weights):
    slots=len(directions)+1
    eye=np.eye(2**slots)
    X=eye.copy()
    for j,w in enumerate(weights):
        for bits in np.ndindex(*(2,)*slots):
            index=np.ravel_multi_index(bits,(2,)*slots)
            moved=list(bits)
            moved[0],moved[j+1]=moved[j+1],moved[0]
            out=np.ravel_multi_index(moved,(2,)*slots)
            X[out,index]+=w
    return X

def full_probe():
    # f,ch,source, then data and all reference systems. Partners are passive.
    r=np.array([1,1j*1.7])/math.sqrt(1+1.7**2)
    w=np.array([1,1])/math.sqrt(2)
    z=np.array([0,1])
    banks=[r,w,z]
    seed=np.kron(r,np.kron(w,z))
    J0=np.kron(np.eye(2),seed[:,None])
    J=np.kron(np.eye(8),J0)
    X=np.array([[0.,1.],[1.,0.]])
    H=np.zeros((128,128),complex)
    Hd=np.zeros((16,16),complex)
    compression_error=0.
    eps=[]
    for i,(t,weights) in enumerate(zip([.2,.25],[[1,2,0],[0,1,3]])):
        X0=small_parent(banks,weights)
        eff=np.eye(2,dtype=complex)
        variance=np.zeros((2,2),complex)
        for v,amount in zip(banks,weights):
            P=np.outer(v,v.conj())
            eff+=amount*P
            variance+=amount**2*(np.eye(2)-P)
        diff=X0@J0-J0@eff
        compression_error=max(compression_error,np.linalg.norm(diff.conj().T@diff-variance,2))
        eps.append(t*math.sqrt(np.linalg.norm(variance,2)))
        block=np.diag([1-i,i])
        H+=np.kron(np.kron(np.kron(block,t*X),np.eye(2)),X0)
        Hd+=np.kron(np.kron(np.kron(block,t*X),np.eye(2)),eff)
    assert compression_error<4e-13
    V=np.zeros((8,8))
    for col in range(8):
        f,ch,src=np.unravel_index(col,(2,2,2))
        out=(src,ch,f) if ch==0 else (f,ch,src)
        V[np.ravel_multi_index(out,(2,2,2)),col]=1.
    assert np.array_equal(V@V,np.eye(8))
    fullV=np.kron(V,np.eye(16))
    effV=np.kron(V,np.eye(2))
    H_controls=math.pi*(np.eye(128)-fullV)/2
    assert np.allclose(evolve(H_controls,1.),fullV)
    a,b=.45,.65
    actual=evolve(H,b)@fullV@evolve(H,a)@J
    ideal=J@evolve(Hd,b)@effV@evolve(Hd,a)
    error=float(np.linalg.norm(actual-ideal,2))
    worst=(a+b)*max(eps)
    assert error<=worst+1e-12
    # Coherently controlled square operators prove the same bound tensored with identity.
    # Evaluate also all final record probabilities on a nonfactorized flavor-partner input.
    partner_seed=np.array([[1.,0.],[0.,1.]])/math.sqrt(2.)
    tiny_in=np.kron(np.eye(8)[:,1,None],partner_seed)
    def records(map_):
        state=(map_@tiny_in).reshape((2,2,2,16,2))
        return np.sum(abs(state)**2,axis=(3,4))
    pa,pb=records(actual),records(ideal)
    assert abs(pa.sum()-1)<1e-12 and abs(pb.sum()-1)<1e-12
    assert np.abs(pa-pb).sum()/2<=error+1e-12
    tau=.01
    actual_pulse=evolve(H*tau+H_controls,1.)
    pulse_error=float(np.linalg.norm(actual_pulse-fullV,2))
    pulse_limit=.25*(1+1+3)*tau
    assert pulse_error<=pulse_limit+1e-12
    return dict(flavor_dimension=2,full_dimension=128,passive_partner_dimension=2,
            variance_residual=float(compression_error),uniform_output_error=error,analytic_uniform_bound=worst,
            joint_final_record_total_variation=float(np.abs(pa-pb).sum()/2),finite_pulse_duration=tau,
            finite_pulse_error=pulse_error,simultaneous_drift_bound=pulse_limit,
            pulse_action_norm=math.pi)

def probability_check():
    # Independently evaluate three-flavor effective readouts by two matrix exponentials.
    r=np.sqrt(np.array([.5,.3,.2]))
    u=np.diag([1,1j,-1]);z=np.diag([1,-1,1])
    dirs=[r,u@r,u@u@r]
    X=np.array([[0,1],[1,0]])
    def H(directions,coef):
        return np.eye(3)+sum(a*np.outer(v,v.conj()) for a,v in zip(coef,directions))
    first=evolve(.2*np.kron(X,H(dirs,[1,2,0])),.45)@np.concatenate([r,0*r])
    left=np.concatenate([first[:3],0*r])
    probs=[]
    for d in (dirs,[z@v for v in dirs]):
        post=evolve(.25*np.kron(X,H(d,[0,1,3])),.65)@left
        probs.append(float(np.vdot(post[3:],post[3:]).real))
    bound=1.10*math.sqrt(10/76)/4
    for single,trials in [(1,(82,776)),(2,(86,808))]:
        n1,n2=trials
        miss=76*(single*Fraction(8,9)**n1+Fraction(80,81)**n2)
        assert miss<Fraction(1,100)
    gap=abs(probs[0]-probs[1])-2*(bound+.01+.0125)
    assert gap>.04
    return dict(shared_effective_readout=probs[0],separate_effective_readout=probs[1],copies_per_bank=76,
                preparation_fail_upper=.01,uniform_propagation_bound=bound,finite_pulse_gap_lower=gap,
                scope="extra finite pulse control assumed, tau=.01, no independent autonomous work-store construction")

def close(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:close(a[k],b[k])
    elif isinstance(a,float):
        assert math.isclose(a,b,abs_tol=2e-10,rel_tol=1e-8)
    else:assert a==b

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    r=dict(round=1340,status="PASS",full_test=full_probe(),source_holdout=probability_check())
    p=Path(__file__).with_name("independent_results.json")
    if args.write_new:
        with p.open("x",encoding="utf8",newline="\n") as stream:json.dump(r,stream,indent=2);stream.write("\n")
    else:close(r,json.loads(p.read_text(encoding="utf8")))
    print(json.dumps(r))

