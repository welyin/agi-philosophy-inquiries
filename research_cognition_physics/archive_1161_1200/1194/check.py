"""1194: finite endpoint probabilities from encoded positive energy.
Default verifies saved output without mutation. No actual detector or Einstein
dynamics is inferred from the operator comparison."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,importlib.util,json,math
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def read(p):return json.loads(p.read_text(encoding="utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(q):return {"exact":str(q),"float":float(q)}
def sqrtm(a):
    w,v=np.linalg.eigh(a)
    assert min(w)>-1e-11
    return (v*np.sqrt(np.maximum(w,0)))@v.conj().T
def halfnorm(a):return float(np.sum(np.abs(np.linalg.eigvalsh((a+a.conj().T)/2)))/2)
def ptrace_q(rho,d,r):return np.einsum("iaib->ab",rho.reshape(d,r,d,r))
def history():
    b=read(HERE/"baseline.json"); p=ROOT/b["prior_progress_path"]
    assert sha(p)==b["prior_progress_sha256"]
    prior=read(p)
    for name,h in prior["science_files_sha256"].items():assert sha(ROOT/name)==h,name
    for n,h in prior["acceptance_sha256"].items():assert sha(HERE.parent/n/"acceptance.json")==h,n
    for name,h in b["reused_files_sha256"].items():assert sha(ROOT/name)==h,name
    return {"science_files_unchanged":len(prior["science_files_sha256"]),
            "acceptances_unchanged":len(prior["acceptance_sha256"]),
            "reused_files_unchanged":len(b["reused_files_sha256"])}

def clock_compression():
    p=ROOT/"research_cognition_physics/archive_1106_1119/1112/check.py"
    spec=importlib.util.spec_from_file_location("r1112_readonly",p)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    m=mod.Material(1);d=m.d;expected=37/2+2*math.sqrt(19)
    residuals=[];variance=[];head_zero=[]
    for t in [0.,.37,1.1,math.pi]:
        a=m.encode_joint(np.eye(d,dtype=complex).reshape(d,1,d),t)
        J=a.reshape(-1,d);HJ=m.h_action(a).reshape(-1,d)
        residuals.append(float(np.linalg.norm(J.conj().T@HJ-expected*np.eye(d),2)))
        variance.append(float(np.linalg.norm(HJ.conj().T@HJ-(expected**2+.75)*np.eye(d),2)))
        head_zero.append(float(np.vdot(a[0,...,0],a[0,...,0]).real))
        assert np.linalg.norm(J.conj().T@J-np.eye(d),2)<1e-12
    assert max(residuals)<2e-11 and max(variance)<3e-10
    assert head_zero[0]>.999999 and head_zero[-1]<1e-12
    return {"logical_dimension":d,"full_device_dimension":288,
            "compressed_energy_exact":"37/2+2*sqrt(19)","compressed_energy_float":expected,
            "compressed_variance_exact":"3/4","normalization_delta_exact":"0",
            "max_energy_residual":max(residuals),"max_second_moment_residual":max(variance),
            "head_zero_probabilities":head_zero,
            "scope":"original 1112 device energy only; neither complete curved body stress nor original CO matching"}

def normalization():
    # Exact sharp reference example; also tests full success/failure dilation.
    E=np.diag([1.,4.]);rootE=np.diag([1.,2.]);K=rootE/2
    Jf=sqrtm(np.eye(2)-K@K)
    psi=np.array([math.sqrt(2/3),0,0,math.sqrt(1/3)],complex)
    rho=np.outer(psi,psi.conj());bigK=np.kron(K,np.eye(2))
    un=bigK@rho@bigK.conj().T;ps=float(np.trace(un).real);success=un/ps
    failure=np.kron(Jf,np.eye(2))@rho@np.kron(Jf,np.eye(2))
    joint_distance=halfnorm(success-rho)
    ref_distance=halfnorm(ptrace_q(success,2,2)-ptrace_q(rho,2,2))
    assert abs(ps-.5)<1e-12 and abs(joint_distance-1/3)<1e-12
    assert abs(ref_distance-1/3)<1e-12
    assert np.linalg.norm(ptrace_q(un+failure,2,2)-ptrace_q(rho,2,2))<1e-12
    # Noncommuting effects on the encoded task; randomized arithmetic is diagnostic.
    rng=np.random.default_rng(1194);max_identity=0.;max_excess=0.;max_cq=0.
    aa=F(999,1000)**2;bb=F(1001,1000)**2
    E=np.diag([float(aa),float(bb)]);S=sqrtm(E);delta=.001
    for _ in range(40):
        z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        rho=z@z.conj().T;rho/=np.trace(rho)
        v=rng.normal(size=2)+1j*rng.normal(size=2);v/=np.linalg.norm(v)
        N=.2*np.eye(2)+.6*np.outer(v,v.conj());A=S@N@S
        en=float(np.trace(np.kron(E,np.eye(2))@rho).real)
        filt=np.kron(S,np.eye(2))@rho@np.kron(S,np.eye(2))/en
        q=float(np.trace(np.kron(A,np.eye(2))@rho).real)/en
        q2=float(np.trace(np.kron(N,np.eye(2))@filt).real)
        max_identity=max(max_identity,abs(q-q2))
        max_excess=max(max_excess,halfnorm(filt-rho)-delta)
        difference=ptrace_q(np.kron(N,np.eye(2))@(filt-rho),2,2)
        # Full two-outcome cq includes complementary outcome; the R marginal
        # of filt need not equal rho's, so the two blocks are computed separately.
        other=ptrace_q(np.kron(np.eye(2)-N,np.eye(2))@(filt-rho),2,2)
        cq=halfnorm(difference)+halfnorm(other)
        max_cq=max(max_cq,cq)
        assert cq<=delta+1e-12
        assert np.linalg.eigvalsh(A).min()>0
        assert np.linalg.norm(S@N@S+S@(np.eye(2)-N)@S-E)<1e-12
    assert max_identity<1e-12 and max_excess<1e-12
    # Total positivity alone does not define positive regional energy effects.
    bad=np.diag([2.,-1.]);bad_complement=np.eye(2)-bad
    assert np.linalg.norm(bad+bad_complement-np.eye(2))==0
    assert bad[1,1]<0
    return {"sharp_reference_example":{"energy":"diag(1,4)","success_probability":"1/2",
              "joint_trace_distance":"1/3","reference_trace_distance":"1/3",
              "full_flag_reference_preserved":True},
            "noncommuting_diagnostic":{"samples":40,"largest_identity_error":max_identity,
              "largest_joint_bound_excess":max_excess,"largest_cq_distance":max_cq,
              "delta_exact":"1/1000"},
            "positive_total_not_local":{"total_energy":"I","regional_minimum":"-1",
              "scope":"algebraic deletion example, not a complete physical stress construction"}}

def endpoint():
    prior=read(HERE.parent/"1193/results.json")
    Z=F(prior["finite_budget"]["position_upper"]["exact"])
    a=F(999,1000)**2;b=F(1001,1000)**2
    delta=F(1,1000);r=F(1,10000);R=F(1,50);eps=F(1,5000)
    tail=r/(R-Z)+delta+eps;threshold=F(7,1000)
    assert R>Z and tail<threshold
    # Existing scalar-compressed device would have zero normalization budget,
    # if the additional positive regional stress/matching hypotheses held.
    no_reweight=r/(R-Z)+eps
    # Positive symmetric three-point distribution: exact endpoint/first moment.
    p=r/(2*R)
    assert 2*p*R==r and 2*p==r/R
    return {"a":enc(a),"b":enc(b),"delta":enc(delta),"r00":enc(r),"R":enc(R),
            "Zx_from_1193":enc(Z),"epsilon_CO":enc(eps),"tail_upper":enc(tail),
            "threshold":enc(threshold),"strict_margin":enc(threshold-tail),
            "scalar_compression_conditional_upper":enc(no_reweight),
            "one_slice_markov_example":{"positions":["-1/50","0","1/50"],
                "probabilities":[str(p),str(1-2*p),str(p)],
                "mean":"0","first_absolute_moment":"1/10000",
                "tail_probability":str(2*p),
                "scope":"single-slice mathematical sharpness; not a dynamical quantum GR implementation"}}

def run():
    hist=history()
    return {"status":"PASS_WITH_SCOPE","check_groups":3,
            "clock_compression":clock_compression(),"normalization":normalization(),
            "endpoint_certificate":endpoint(),"frozen_history":hist,
            "new_adopted_axioms":0,"new_empirical_calibration":0,"gr_hard_deliverables_completed":0}
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-8,abs_tol=5e-12),(a,b)
    else:assert a==b,(a,b)
if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--write-results",action="store_true");args=ap.parse_args()
    d=run();p=HERE/"results.json"
    if args.write_results:
        with p.open("x",encoding="utf-8") as f:json.dump(d,f,ensure_ascii=False,indent=2);f.write("\n")
    else:compare(d,read(p))
    print(json.dumps({"status":d["status"],"groups":3,
      "tail_upper":d["endpoint_certificate"]["tail_upper"]["float"],
      "strict_margin":d["endpoint_certificate"]["strict_margin"]["float"],
      "history":d["frozen_history"]},ensure_ascii=False))
