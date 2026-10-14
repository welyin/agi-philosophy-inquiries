"""1194 internal non-blind independent review.
Does not import or execute 1194/check.py or the old Material implementation.
Uses exact spin-clock spectral moments, fixed noncommuting effects/entangled
references, analytic sharp filtering, and an independent rational endpoint bound.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def enc(q): return {"exact":str(q),"float":float(q)}
def read(p): return json.loads(p.read_text(encoding="utf-8-sig"))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def normhalf(a):
    a=(a+a.conj().T)/2
    return float(np.abs(np.linalg.eigvalsh(a)).sum()/2)
def trace_code(a):
    # Independent explicit block trace for code dimension 2, reference dimension 2.
    return a[:2,:2]+a[2:,2:]
def sqrt_positive(a):
    vals,vec=np.linalg.eigh(a)
    assert min(vals)>-1e-12
    return (vec*np.sqrt(np.maximum(vals,0)))@vec.conj().T

def clock_moments():
    # Hc=3I/2+Jx^(3/2), initial |0_z>; its spectral weights are binomial.
    weights=[F(1,8),F(3,8),F(3,8),F(1,8)]
    levels=list(map(F,range(4)))
    mean=sum(w*x for w,x in zip(weights,levels))
    second=sum(w*x*x for w,x in zip(weights,levels))
    variance=second-mean*mean
    assert mean==F(3,2) and variance==F(3,4)
    # Old full baseline is 17+2sqrt(19); no full-device eigenstate assumption.
    total_rational,total_radical=F(17)+mean,F(2)
    assert total_rational==F(37,2)
    second_rational=total_rational**2+19*total_radical**2+variance
    second_radical=2*total_rational*total_radical
    assert (second_rational,second_radical)==(F(419),F(74))
    data=read(HERE/"results.json")["clock_compression"]
    assert data["compressed_energy_exact"]=="37/2+2*sqrt(19)"
    assert data["compressed_variance_exact"]==str(variance)
    probs=[math.cos(t/2)**6 for t in (0.,.37,1.1,math.pi)]
    discrepancy=max(abs(a-b) for a,b in zip(probs,data["head_zero_probabilities"]))
    assert discrepancy<2e-12
    return {"method":"exact binomial spectral moments plus frozen FK intertwining",
            "head_mean":str(mean),"head_variance":str(variance),
            "compressed_energy":"37/2+2*sqrt(19)",
            "compressed_second_moment":"419+74*sqrt(19)",
            "head_probability_identity_max_error":discrepancy,
            "scope":"old device energy compression only; curved stress and physical filtering not identified"}

def sharp_filter():
    # E=diag(1,4). Initial Schmidt probabilities (2/3,1/3); success swaps them.
    pin=[F(2,3),F(1,3)]
    energy=[F(1),F(4)]
    h=sum(x*y for x,y in zip(pin,energy))
    pout=[x*y/h for x,y in zip(pin,energy)]
    assert pout==pin[::-1] and h/F(4)==F(1,2)
    # Squared overlap = (2sqrt(2)/3)^2; pure-state trace distance =1/3.
    overlap_squared=F(8,9)
    assert 1-overlap_squared==F(1,9)
    ref_distance=sum(abs(x-y) for x,y in zip(pin,pout))/2
    assert ref_distance==F(1,3)
    pfail=F(1)-h/4
    full_reference=[h/F(4)*pout[0]+pfail,h/F(4)*pout[1]]
    assert full_reference==pin
    # General elementary spectral inequality for sqrt(E) in [s,t]:
    # <E> <= (s+t)<sqrt(E)>-st; complete square proves the fidelity bound.
    # (s+t)^2 m^2 -4st[(s+t)m-st] = ((s+t)m-2st)^2 >=0.
    for s,t,m in ((F(1),F(2),F(4,3)),(F(999,1000),F(1001,1000),F(1))):
        left=(s+t)**2*m*m-4*s*t*((s+t)*m-s*t)
        assert left==((s+t)*m-2*s*t)**2
    return {"success_probability":"1/2","conditional_joint_distance":"1/3",
            "conditional_reference_distance":str(ref_distance),
            "full_flag_reference_preserved":True,
            "general_bound_proof":"purification overlap and completed-square spectral inequality"}

def fixed_noncommuting_cases():
    s,t=F(999,1000),F(1001,1000)
    delta=float((t-s)/(t+s)); epsilon=.0002
    rotation=np.array([[.6,-.8],[.8,.6]],complex)
    S=rotation@np.diag([float(s),float(t)])@rotation.conj().T
    E=S@S
    N=np.array([[.5,-.25j],[.25j,.5]],complex)
    A=S@N@S
    assert np.linalg.norm(E@N-N@E)>1e-4
    assert np.linalg.eigvalsh(A).min()>0
    assert np.linalg.eigvalsh(E-A).min()>0
    vectors=[np.array(v,complex) for v in
             ([1,1j,2,-1],[1,0,0,1],[0,2,1j,1],[3,-1,1,2j])]
    rhos=[]
    for psi in vectors:
        psi/=np.linalg.norm(psi)
        rhos.append(np.outer(psi,psi.conj()))
    rhos.append(.4*rhos[0]+.6*rhos[2])
    actual=N+epsilon*np.diag([1.,-1.])
    assert np.linalg.eigvalsh(actual).min()>0 and np.linalg.eigvalsh(actual).max()<1
    maxima={"identity":0.,"joint":0.,"cq":0.,"effect_excess":0.,"full_flag_R":0.}
    bigS=np.kron(S,np.eye(2)); bigN=np.kron(N,np.eye(2))
    Ks=S/float(t); Kf=sqrt_positive(np.eye(2)-Ks@Ks)
    assert np.linalg.norm(Ks.conj().T@Ks+Kf.conj().T@Kf-np.eye(2))<2e-12
    for rho in rhos:
        h=float(np.trace(np.kron(E,np.eye(2))@rho).real)
        weighted=bigS@rho@bigS.conj().T/h
        q1=float(np.trace(np.kron(A,np.eye(2))@rho).real)/h
        q2=float(np.trace(bigN@weighted).real)
        dist=normhalf(weighted-rho)
        block0=trace_code(bigN@(weighted-rho))
        block1=trace_code(np.kron(np.eye(2)-N,np.eye(2))@(weighted-rho))
        cq=normhalf(block0)+normhalf(block1)
        pco=float(np.trace(np.kron(actual,np.eye(2))@rho).real)
        full=(np.kron(Ks,np.eye(2))@rho@np.kron(Ks,np.eye(2)).conj().T
              +np.kron(Kf,np.eye(2))@rho@np.kron(Kf,np.eye(2)).conj().T)
        flag_error=np.linalg.norm(trace_code(full)-trace_code(rho))
        assert abs(q1-q2)<2e-12 and dist<=delta+2e-12 and cq<=delta+2e-12
        assert pco<=q1+delta+epsilon+2e-12 and flag_error<2e-12
        for key,val in (("identity",abs(q1-q2)),("joint",dist),("cq",cq),
                        ("effect_excess",pco-q1-delta-epsilon),("full_flag_R",flag_error)):
            maxima[key]=max(maxima[key],float(val))
    # Positive total does not imply positive regional operator measure.
    bad=np.diag([2.,-1.]); complement=np.eye(2)-bad
    assert np.array_equal(bad+complement,np.eye(2)) and np.linalg.eigvalsh(bad).min()<0
    return {"fixed_cases":len(rhos),"delta_exact":"1/1000","maxima":maxima,
            "effects_not_commuting_with_energy":True,
            "positive_total_insufficient":True,
            "scope":"finite mathematical effects and states; no actual CO hardware or common source"}

def endpoint():
    result=read(HERE/"results.json")["endpoint_certificate"]
    prior=read(HERE.parent/"1193/results.json")
    Z=F(prior["finite_budget"]["position_upper"]["exact"])
    s,t=F(999,1000),F(1001,1000)
    delta=(t-s)/(t+s)
    r,R,epsilon=F(1,10000),F(1,50),F(1,5000)
    bound=r/(R-Z)+delta+epsilon
    threshold=F(7,1000)
    assert bound<threshold
    assert str(bound)==result["tail_upper"]["exact"]
    assert str(threshold-bound)==result["strict_margin"]["exact"]
    # Independent coarser source bound, proven already in 1193 proof:
    Zcoarse=F(1207,10**8)
    assert Z<Zcoarse<R
    conservative=r/(R-Zcoarse)+delta+epsilon
    assert bound<=conservative<threshold
    # Closed-tail Markov example; outcome |x-Y| >= R includes the endpoints.
    masses=[F(1,400),F(199,200),F(1,400)]
    positions=[-R,F(0),R]
    first=sum(p*abs(x) for p,x in zip(masses,positions))
    tail=sum(p for p,x in zip(masses,positions) if abs(x)>=R)
    assert first==r and tail==r/R
    # Stronger joint-R bound can use any positive reference effect Omega<=I,
    # since Tr_R[Omega times positive outside cq block] <= total tail.
    return {"exact_author_budget":enc(bound),"exact_margin":enc(threshold-bound),
            "independent_conservative_bound":enc(conservative),"fixed_threshold":"7/1000",
            "closed_tail_example":{"r00":str(first),"tail":str(tail)},
            "fixed_center_and_same_state_required":True}

def bindings():
    paths=[HERE/name for name in
           ("proof.md","admission.md","sources.md","check.py","results.json","baseline.json","NEXT.md")]
    paths.append(HERE.parent/"research_note_1194.md")
    hashes={}
    for p in paths:
        assert p.is_file(),str(p)
        hashes[p.relative_to(HERE.parent).as_posix()]=sha(p)
    baseline=read(HERE/"baseline.json")
    for name,expected in baseline["reused_files_sha256"].items():
        assert sha(ROOT/name)==expected,name
        hashes["reused/"+name]=expected
    return hashes

def run():
    return {"status":"PASS_WITH_SCOPE","review_type":"internal non-blind independent mathematical review",
            "clock":clock_moments(),"sharp_filter":sharp_filter(),
            "noncommuting":fixed_noncommuting_cases(),"endpoint":endpoint(),
            "bound_files_sha256":bindings(),
            "scope":{"finite_probability_bridge":True,
                     "original_CO_instrument_implemented":False,
                     "proper_time_calibration_certified":False,
                     "physical_energy_filter_nondisturbing":False,
                     "curved_stress_equals_clock_H_certified":False,
                     "Einstein_dynamics_generated":False,"gr_goal_complete":False}}

def equivalent(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:equivalent(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):equivalent(x,y)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-10,abs_tol=2e-13),(a,b)
    else:
        assert a==b,(a,b)

if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--write",action="store_true")
    args=ap.parse_args();data=run();target=HERE/"review_checks.json"
    if args.write:
        with target.open("x",encoding="utf-8") as f:
            json.dump(data,f,ensure_ascii=False,indent=2);f.write("\n")
    else:equivalent(data,read(target))
    print(json.dumps({"status":data["status"],
                      "tail_upper":data["endpoint"]["exact_author_budget"]["float"],
                      "independent_upper":data["endpoint"]["independent_conservative_bound"]["float"],
                      "bound_files":len(data["bound_files_sha256"])},ensure_ascii=False))
