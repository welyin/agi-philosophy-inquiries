"""Round 1261, one author parent-symbol group. Default replay is read-only."""
from fractions import Fraction as F
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

def null_contraction(P, kvec, n, mass):
    k = float(np.linalg.norm(kvec))
    if k == 0:
        raise ValueError("This screen has a positive-energy signal band.")
    E = np.sqrt(mass*mass + np.dot(P,P))
    u = kvec/k
    c = float(np.dot(n,u)); y = abs(c)
    a = float(np.dot(n,P)); d = float(np.dot(P,u))
    A = float(np.dot(P,P))
    Q = A-3*a*a+8*a*c*d+y*y*(a*a-3*A)+2*y*(2*d*d-A+a*a)
    return k/E * Q/(1+y)**2

def ordered_contraction(pa, pb, ma, mb, n):
    if mb == 0:
        if ma == 0:
            raise ValueError("Only one massless signal is in this screen.")
        return null_contraction(pa, pb, n, ma)
    pa2 = float(np.dot(pa,pa)); pb2 = float(np.dot(pb,pb))
    Ea = np.sqrt(ma*ma+pa2); Eb = np.sqrt(mb*mb+pb2)
    an = float(np.dot(pa,n)); bn = float(np.dot(pb,n))
    d = float(np.dot(pa,pb))
    y = np.sqrt(mb*mb+bn*bn)/Eb
    N = 2*(2*d*d*bn*bn-2*an*bn*d*pb2+an*an*pb2*pb2-d*d*pb2)/(Eb*Eb)
    N += 2*(-pa2*bn*bn+an*an*bn*bn+2*an*bn*d+d*d-an*an*pb2)
    N += (-3*pa2*bn*bn+an*an*bn*bn+8*an*bn*d+pa2*pb2-3*an*an*pb2)*y
    return N/(Ea*Eb*y*(1+y)**2)

def pair(pa,pb,ma,mb,n,r=1.0,G=1.0):
    Ea=np.sqrt(ma*ma+np.dot(pa,pa)); Eb=np.sqrt(mb*mb+np.dot(pb,pb))
    out=-G*Ea*Eb/r*(1+np.dot(pa,pa)/Ea**2+np.dot(pb,pb)/Eb**2)
    out+=G/(2*r)*(7*np.dot(pa,pb)+np.dot(pa,n)*np.dot(pb,n))
    out-=G/(4*r)*(ordered_contraction(pa,pb,ma,mb,n)
                         +ordered_contraction(pb,pa,mb,ma,-n))
    return float(out)


def build_diagnostic():
    # Pure polynomial cancellation: c^2=y^2, including c=0.
    count=0
    for c in [F(-1),F(-3,5),F(-1,7),F(0),F(1,7),F(3,5),F(1)]:
        y=abs(c)
        for A,a,d in [(F(7),F(1,3),F(-2,5)),(F(1),F(1),F(0)),(F(3),F(-1),F(1))]:
            N=2*(2*d*d*c*c-2*a*c*d+a*a-d*d)
            N+=2*(-A*c*c+a*a*c*c+2*a*c*d+d*d-a*a)
            N+=(-3*A*c*c+a*a*c*c+8*a*c*d+A-3*a*a)*y
            Q=A-3*a*a+8*a*c*d+y*y*(a*a-3*A)+2*y*(2*d*d-A+a*a)
            assert N==y*Q
            count+=1
    # Coefficient bounds, not a numerical substitute for the continuous proof.
    assert F(31)*64/147 < 14
    assert 3 + F(6+14,4) == 8
    assert F(4)+16*F(2,3)+4*F(2,3)**2 <= 6*(1+F(2,3))**2

    rng=np.random.default_rng(1260)
    null_ratio=reverse_ratio=pair_ratio=direct_error=0.0
    for _ in range(1200):
        m=float(rng.choice([27,32,1000]))
        v=float(rng.uniform(1e-4,.5))
        P=rng.normal(size=3); P=P/np.linalg.norm(P)*(m*v/np.sqrt(1-v*v))
        k=float(rng.uniform(.5,2))
        u=rng.normal(size=3);u/=np.linalg.norm(u)
        n=rng.normal(size=3);n/=np.linalg.norm(n)
        E=np.sqrt(m*m+np.dot(P,P)); A=np.dot(P,P)
        kvec=k*u
        nc=null_contraction(P,kvec,n,m)
        rc=ordered_contraction(kvec,P,0,m,-n)
        null_ratio=max(null_ratio,abs(nc)/(k*A/E))
        reverse_ratio=max(reverse_ratio,abs(rc)/(k*A/E))
        assert abs(nc)<=6*k*A/E*(1+1e-12)
        assert abs(rc)<=14*k*A/E*(1+1e-12)
        pv=pair(P,kvec,m,0,n)
        linear=k/2*(7*np.dot(P,u)+np.dot(P,n)*np.dot(u,n))
        residual=pv+2*m*k-linear
        pair_ratio=max(pair_ratio,abs(residual)/(E*k*v*v))
        assert abs(residual)<=8*E*k*v*v+1e-9
        # Direct unreduced massless expression away from the removable 1/y.
        c=float(np.dot(n,u)); y=abs(c)
        if y>1e-3:
            a=float(np.dot(P,n));d=float(np.dot(P,u))
            N=2*k*k*(2*d*d*c*c-2*a*c*d+a*a-d*d)
            N+=2*k*k*(-A*c*c+a*a*c*c+2*a*c*d+d*d-a*a)
            N+=k*k*(-3*A*c*c+a*a*c*c+8*a*c*d+A-3*a*a)*y
            direct=N/(E*k*y*(1+y)**2)
            direct_error=max(direct_error,abs(direct-nc)/max(1,abs(nc)))
            assert abs(direct-nc)<1e-8*max(1,abs(nc))

    # The full finite-mass ordered sum generally has a cusp, not a divergence.
    P=np.array([1.,0,0]);kvec=np.array([0.,0,2.]);m=3.
    def at(c):
        n=np.array([np.sqrt(1-c*c),0.,c])
        return pair(P,kvec,m,0,n)
    step=1e-6
    jump=(at(step)+at(-step)-2*at(0))/step
    expected=-2*2/np.sqrt(10)
    assert abs(jump-expected)<2e-5

    # Translational reaction, with all three positions/momenta present.
    # Values here are algebra diagnostics, not a weak-field physical task.
    q=np.array([[0.,0,0],[4.,2,-1],[-2.,1,3]])
    p=np.array([[.1,.2,-.1],[.2,-.1,.3],[1.,.4,.2]])
    masses=[27.,1000.,0.]
    def interaction(q):
        value=0.
        for a in range(3):
            for b in range(a+1,3):
                d=q[b]-q[a];r=np.linalg.norm(d)
                value+=pair(p[a],p[b],masses[a],masses[b],d/r,r)
        return value
    forces=np.zeros((3,3));h=1e-4
    for a in range(3):
        for j in range(3):
            qp=q.copy();qm=q.copy();qp[a,j]+=h;qm[a,j]-=h
            forces[a,j]=-(interaction(qp)-interaction(qm))/(2*h)
    reaction=float(np.linalg.norm(forces.sum(axis=0)))
    assert reaction<2e-6
    diagnostic = {
     "status":"PASS_1261_PARENT_SYMBOL_GROUP",
     "fraction_cancellations":count,
     "random_diagnostics":1200,
     "max_null_contraction_ratio_bound6":null_ratio,
     "max_reverse_contraction_ratio_bound14":reverse_ratio,
     "max_pair_remainder_ratio_bound8":pair_ratio,
     "direct_massless_formula_relative_error":direct_error,
     "angular_derivative_jump":jump,
     "predicted_jump":expected,
     "three_body_total_force_residual":reaction,
     "scope":"Finite symbol identities and diagnostic samples only; not source, apparatus, GR remainder, or actual arrival certification."
    }

    return diagnostic

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def build_results():
    here=Path(__file__).resolve().parent
    stage=here.parent
    return {
        "round":1261,
        "scientific_groups":1,
        "group_id":"parent_allvelocity_symbol",
        "bindings":{
            "parent_proof_sha256":digest(here/"parent.md"),
            "parent_code_sha256":digest(Path(__file__).resolve()),
            "frozen_signal_bridge_sha256":digest(stage/"_admission/same_source_signal_bridge_after1260/README.md"),
            "inherited_1260_proof_sha256":digest(stage/"1260/proof.md")
        },
        "diagnostic":build_diagnostic(),
        "continuous_proofs":"parent.md; random samples are diagnostics only",
        "actual_arrival_certified":False,
        "complete_source_insertion_certified":False,
        "GR_remainder_certified":False
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true",help="first exclusive result creation only")
    args=parser.parse_args()
    result=build_results()
    path=Path(__file__).resolve().with_name("parent_results.json")
    if args.write:
        with path.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write("\n")
        print("WROTE_ONCE_1261_PARENT")
    else:
        expected=json.loads(path.read_text(encoding="utf-8"))
        if result!=expected:
            raise AssertionError("1261 parent replay does not match saved result")
        print("READ_ONLY_MATCH_1261_PARENT; scientific_groups=1")

if __name__=="__main__":
    main()
