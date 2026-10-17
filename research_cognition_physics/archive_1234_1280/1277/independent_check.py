"""Independent 1277: channel matrices and exact finite reachability bounds.
No author implementation is imported. Default is read-only; --write is exclusive.
Old channel identities are inherited through hashed proofs, not rederived as a
new carrier. Matrix samples diagnose algebra; the analytic bounds are general.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "research_cognition_physics").is_dir())

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def need(test, message):
    if not test:
        raise AssertionError(message)

def rat(value):
    v = F(value)
    return {"numerator": v.numerator, "denominator": v.denominator}

def clean(value):
    return float(f"{float(value):.12g}")

def trdist(a, b):
    d = (a-b + (a-b).conj().T)/2
    return float(np.sum(np.abs(np.linalg.eigvalsh(d)))/2)

def polyadd(*terms):
    n = max(map(len, terms))
    return [sum((p[j] if j < len(p) else F(0)) for p in terms) for j in range(n)]

def scale(p, a):
    return [a*x for x in p]

def q2mul(a, b):
    return (a[0]*b[0]+2*a[1]*b[1], a[0]*b[1]+a[1]*b[0])

def q2eval(poly, x):
    out = (F(0), F(0))
    for c in reversed(poly):
        out = q2mul(out, x)
        out = (out[0]+c, out[1])
    return out

def rational_checks():
    # Absolute-value quadratics obtained by summing the original sector spectra.
    A = [F(-1,4), F(0), F(1,2)]
    B = [F(1,10), F(3,5), F(-4,5)]
    C = [F(1,10), F(-2,5), F(1,5)]
    signs = [(-1,1,1), (-1,1,-1), (1,1,-1), (1,-1,-1)]
    pieces = [polyadd(scale(A,a),scale(B,b),scale(C,F(3,2)*c))
              for a,b,c in signs]
    need(pieces == [[F(1,2),F(0),F(-1)],
                    [F(1,5),F(6,5),F(-8,5)],
                    [F(-3,10),F(6,5),F(-3,5)],
                    [F(-1,2),F(0),F(1)]], "piecewise polynomials")
    sqrt2lo, sqrt2hi = F(1414213,10**6), F(1414214,10**6)
    sqrt17lo, sqrt17hi = F(4123105,10**6), F(4123106,10**6)
    need(sqrt2lo**2 < 2 < sqrt2hi**2, "sqrt2 enclosure")
    need(sqrt17lo**2 < 17 < sqrt17hi**2, "sqrt17 enclosure")
    # b=1-1/sqrt2, a=1/sqrt2, c=(3+sqrt17)/8.
    blo,bhi = 1-1/sqrt2lo,1-1/sqrt2hi
    alo,ahi = 1/sqrt2hi,1/sqrt2lo
    clo,chi = (3+sqrt17lo)/8,(3+sqrt17hi)/8
    need(0 < blo < bhi < alo < ahi < clo < chi < 1, "ordered roots")
    b,a = (F(1),F(-1,2)),(F(0),F(1,2))
    need(q2eval(pieces[0],b)==q2eval(pieces[1],b)==(F(-1),F(1)),
         "continuity at b")
    dstar = (F(-3,5),F(3,5))
    need(q2eval(pieces[1],a)==q2eval(pieces[2],a)==dstar,
         "exact minimizer value")
    # First piece decreases, second is concave (its minimum is at an endpoint),
    # third increases up to u=1, fourth increases on [c,1].
    need(pieces[0][2]<0 and pieces[1][2]<0 and
         pieces[2][1]==-2*pieces[2][2] and pieces[3][2]>0,
         "analytic monotonicity/concavity coefficients")
    need(q2eval(pieces[0],b)[0]==-1 and F(1)-F(3,5)>0,
         "D(b)-D(a)=2(sqrt2-1)/5 positive")
    dlow,dhigh = F(3,5)*(sqrt2lo-1),F(3,5)*(sqrt2hi-1)
    need(dlow > F(6,25), "positive fixed-time gap")
    # Pure Schmidt id-depolarizing difference: block trace=1/2,
    # determinant=-3p(1-p)/4, other eigenvalues -p/2 and -(1-p)/2.
    schmidt_cases=0
    for j in range(101):
        p=F(j,100)
        determinant=(p/2)*((1-p)/2)-p*(1-p)
        need(determinant == -F(3,4)*p*(1-p), "Schmidt block determinant")
        need(0<=p*(1-p)<=F(1,4), "Schmidt discriminant bound")
        need(1+12*p*(1-p)<=4, "upper bound 3/4")
        schmidt_cases+=1
    # Rank-three Bell-complement witness inside P4 x R.
    bell_witness=3*F(1,4)*F(1,2)
    need(bell_witness==F(3,8)==F(1,2)*F(3,4), "matching half-diamond bounds")
    return {
      "pieces_ascending_coefficients":[[rat(c) for c in p] for p in pieces],
      "ordered_breakpoints":["1-1/sqrt(2)","1/sqrt(2)","(3+sqrt(17))/8"],
      "minimizer":"u=1/sqrt(2)",
      "fixed_time_minimum":"3*(sqrt(2)-1)/5",
      "fixed_time_minimum_rational_interval":[rat(dlow),rat(dhigh)],
      "half_diamond_random_wait_lower":rat(bell_witness),
      "endpoint_mixture_half_diamond_upper":rat(F(3,8)),
      "Schmidt_exact_cases":schmidt_cases,
      "general_schmidt_half_trace_formula":"(1+sqrt(1+12*p*(1-p)))/4",
      "proof_scope":"piecewise signs/convexity and Schmidt reduction, not a grid proof"
    }

def basis_maps():
    # Choose a basis in which W's two-dimensional image is the first two H6
    # coordinates. This is a unitary coordinate choice for the inherited W,
    # not a new physical implementation.
    e4=np.zeros((7,2),complex);e4[:2]=np.eye(2)
    e6=np.zeros((7,2),complex);e6[2:4]=np.eye(2)
    target=np.diag([.25]*2+[.1]*5).astype(complex)
    return e4,e6,target

def kraus(weights):
    e4,e6,_=basis_maps()
    p0,p1,p2=weights
    result=[math.sqrt(max(0,p0))*e4,math.sqrt(max(0,p1))*e6]
    for j in range(5):
        for k in range(2):
            m=np.zeros((7,2),complex)
            m[2+j,k]=math.sqrt(max(0,p2)/5)
            result.append(m)
    return result

def apply(rho,weights,rdim=1):
    out=np.zeros((7*rdim,7*rdim),complex)
    for k in kraus(weights):
        k=np.kron(k,np.eye(rdim))
        out+=k@rho@k.conj().T
    return out

def state_at(weights,rho):
    return apply(rho,weights)

def entropy(rho):
    eig=np.linalg.eigvalsh((rho+rho.conj().T)/2)
    eig=eig[eig>1e-13]
    return float(-np.dot(eig,np.log(eig)))

def matrix_checks():
    e4,e6,target=basis_maps()
    rng=np.random.default_rng(1277001)
    chain=np.array([[0,1,0],[1,0,1],[0,1,0]],float)
    vals,vecs=np.linalg.eigh(chain)
    tstar=math.pi/math.sqrt(2)
    tmin=math.acos(math.sqrt(2)-1)/math.sqrt(2)
    times=list(np.linspace(0,2*tstar,31))+[tmin]
    max_channel_formula=0.;max_completeness=0.;max_chain=0.
    fixed_samples=[]
    for t in times:
        amp=vecs@(np.exp(-1j*t*vals)*vecs[0,:].conj())
        weights=np.abs(amp)**2
        u=(1+math.cos(math.sqrt(2)*t))/2
        exact=np.array([u*u,2*u*(1-u),(1-u)**2])
        max_chain=max(max_chain,float(np.max(np.abs(weights-exact))))
        ks=kraus(weights)
        max_completeness=max(max_completeness,float(np.max(np.abs(
            sum(k.conj().T@k for k in ks)-np.eye(2)))))
        rho=state_at(weights,np.eye(2)/2)
        value=trdist(rho,target)
        formula=abs(u*u/2-.25)+abs(.1+.6*u-.8*u*u)+1.5*abs(.1-.4*u+.2*u*u)
        max_channel_formula=max(max_channel_formula,abs(value-formula))
        fixed_samples.append([clean(t),clean(value)])
    need(max_chain<1e-13 and max_completeness<1e-13 and max_channel_formula<1e-13,
         "direct autonomous amplitudes/Kraus channel")
    dstar=3*(math.sqrt(2)-1)/5
    need(abs(fixed_samples[-1][1]-dstar)<1e-11,"direct minimizer matrix")

    bell=np.array([1,0,0,1],complex)/math.sqrt(2)
    rb=np.outer(bell,bell.conj())
    p4r=np.kron(e4@e4.conj().T,np.eye(2))
    bellout=np.kron(e4,np.eye(2))@bell
    witness=p4r-np.outer(bellout,bellout.conj())
    target_bell=np.kron(target,np.eye(2)/2)
    need(abs(np.trace(witness@target_bell).real-.375)<1e-14,"witness target")
    min_random_gap=1.;max_witness_actual=0.
    for _ in range(40):
        ts=rng.uniform(0,2*tstar,size=13)
        mix=rng.dirichlet(np.ones(13))
        weights=np.zeros(3)
        for w,t in zip(mix,ts):
            a=vecs@(np.exp(-1j*t*vals)*vecs[0,:].conj())
            weights+=w*np.abs(a)**2
        out=apply(rb,weights,2)
        max_witness_actual=max(max_witness_actual,abs(np.trace(witness@out)))
        gap=trdist(out,target_bell)
        min_random_gap=min(min_random_gap,gap)
        need(gap>=.375-1e-13,"random waiting lower bound diagnostic")
    need(max_witness_actual<1e-13,"Bell orthogonal effect annihilates P4 output")

    max_schmidt_formula_error=0.;max_endpoint_gap=0.; schmidt_matrix_cases=0
    for p in np.linspace(0,1,41):
        base=np.diag([math.sqrt(p),math.sqrt(1-p)])
        for _ in range(3):
            ua,_=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))
            ur,_=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))
            amplitude=ua@base@ur.T
            vector=amplitude.reshape(-1)
            rho=np.outer(vector,vector.conj())
            ref=amplitude.T@amplitude.conj()
            actual=apply(rho,[.5,0,.5],2)
            comparison=np.kron(target,ref)
            gap=trdist(actual,comparison)
            formula=(1+math.sqrt(1+12*p*(1-p)))/8
            max_schmidt_formula_error=max(max_schmidt_formula_error,abs(gap-formula))
            max_endpoint_gap=max(max_endpoint_gap,gap)
            schmidt_matrix_cases+=1
    need(max_schmidt_formula_error<1e-12 and abs(max_endpoint_gap-.375)<1e-12,
         "endpoint mixture upper bound diagnostic")

    # Arbitrary input eigenvalues: convexity/symmetry proves the general lower
    # bound; these matrices independently verify the claimed spectral use.
    unknown_cases=0;min_excess=1.
    for _ in range(120):
        z=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        rho=z@z.conj().T;rho/=np.trace(rho)
        u=float(rng.uniform())
        weights=[u*u,2*u*(1-u),(1-u)**2]
        gap=trdist(state_at(weights,rho),target)
        uniform=trdist(state_at(weights,np.eye(2)/2),target)
        min_excess=min(min_excess,gap-uniform)
        need(gap>=uniform-1e-13,"unknown input convex lower bound")
        unknown_cases+=1
    rho4=e4@e4.conj().T/2
    rho6=np.diag([0,0]+[.2]*5)
    cq=np.zeros((14,14),complex)
    cq[:7,:7]=rho4/2;cq[7:,7:]=rho6/2
    b=(rho4+rho6)/2
    mutual=math.log(2)+entropy(b)-entropy(cq)
    need(trdist(b,target)<1e-14 and abs(mutual-math.log(2))<1e-13,
         "known-input scheduler record correlation")
    return {
      "direct_chain_times":len(times),"max_chain_probability_error":clean(max_chain),
      "max_Kraus_completeness_error":clean(max_completeness),
      "max_trace_formula_error":clean(max_channel_formula),
      "fixed_time_samples":fixed_samples,
      "random_wait_Bell_checks":40,
      "max_Bell_witness_actual_probability":clean(max_witness_actual),
      "min_random_wait_Bell_distance":clean(min_random_gap),
      "endpoint_Schmidt_matrix_cases":schmidt_matrix_cases,
      "max_Schmidt_formula_error":clean(max_schmidt_formula_error),
      "endpoint_half_diamond_display":clean(max_endpoint_gap),
      "arbitrary_input_checks":unknown_cases,
      "min_arbitrary_input_excess_over_uniform":clean(min_excess),
      "known_input_mixture_reference_error":clean(trdist(b,target)),
      "scheduler_B_mutual_information_display":clean(mutual),
      "diagnostics_not_continuous_proof":True
    }

def scientific():
    return {"exact_certificates":rational_checks(),"independent_matrix_diagnostics":matrix_checks()}

def calculate():
    baseline=json.loads((HERE/"source_baseline.json").read_text(encoding="utf-8-sig"))
    for rel,expected in baseline["files_sha256"].items():
        need(sha(ROOT/rel)==expected,"old source binding changed: "+rel)
    names=["admission.md","source_baseline.json","proof.md","check.py","results.json",
           "NEXT.md","independent_check.py"]
    binding={str((HERE/n).relative_to(ROOT)).replace("\\","/"):sha(HERE/n) for n in names}
    note=HERE.parent/"research_note_1277.md"
    binding[str(note.relative_to(ROOT)).replace("\\","/")]=sha(note)
    return {
      "round":1277,"status":"PASS_WITH_SCOPE","independent_groups":1,
      "source_bindings":binding,"inherited_sources_verified":len(baseline["files_sha256"]),
      **scientific(),
      "scope":{
        "original_complete_process_and_input_code_retained":True,
        "fixed_time_lower_bound_covers_all_original_B_inputs":True,
        "random_scheduler_is_only_declared_mathematical_extension":True,
        "universal_reference_preparation_is_stronger_optional_contract":True,
        "reference_global_rank_used_in_place_of_B_bound":False,
        "half_diamond_norm_uses_arbitrary_passive_reference":True,
        "known_uniform_input_can_reach_reference_marginal_by_endpoint_mixture":True,
        "record_independence_inferred_from_B_marginal":False,
        "new_scheduler_measurement_or_physical_source_certified":False,
        "all_thermal_references_ruled_out":False,
        "geometric_or_Jacobson_identity_completed":False
      }
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--write",action="store_true")
    args=p.parse_args()
    actual=calculate()
    out=HERE/"independent_results.json"
    if args.write:
        with out.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(actual,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
        mode="EXCLUSIVE_WRITE"
    else:
        expected=json.loads(out.read_text(encoding="utf-8"))
        need(actual==expected,"READ_ONLY_MISMATCH")
        mode="READ_ONLY_MATCH"
    print("PASS_WITH_SCOPE "+mode)

if __name__=="__main__":
    main()

