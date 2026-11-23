"""Independent round 1352. Mixed-tensor neutrality, full-group vs. Lie-algebra tests."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np


def su3_generators():
    out=[]
    for i,j in itertools.combinations(range(3),2):
        E=np.zeros((3,3),complex); E[i,j]=1
        out.extend([E+E.T,1j*(E-E.T)])
    return out+[np.diag([1.,-1.,0.]),np.diag([1.,1.,-2.])]


def so3_generators():
    return [T for T in su3_generators() if np.max(abs(T.imag))>0]


def invariant_projector(charges, reps, scalar=None):
    length=len(reps)
    dim=3**length
    quadratic=np.zeros((dim,dim),complex)
    for T in charges:
        total=np.zeros_like(quadratic)
        for pos in range(length):
            term=np.array([[1.]])
            for i,ref in enumerate(reps):
                U=T if ref==1 else -T.T
                term=np.kron(term,U if i==pos else np.eye(3))
            total+=term
        quadratic+=total@total
    if scalar is not None:
        phase=scalar**sum(reps)
        quadratic+=abs(phase-1)**2*np.eye(dim)
    eig,vec=np.linalg.eigh(quadratic)
    zero=vec[:,eig<1e-9]
    return zero@zero.conj().T


def ranks(group,scalar=None):
    return {
        pattern:int(round(np.trace(invariant_projector(group,reps,scalar)).real))
        for pattern,reps in [("RR",[1,1]),("RbarR",[1,-1]),("RRR",[1,1,1])]
    }


def valuation_test(N):
    # Character weight k has minimum possible total tensor degree sum_i |k_i|.
    # Positive or negative counts are supplied by a factor R or its conjugate.
    assert N%3
    low_bad=[]
    for k in itertools.product(range(-N,N+1),repeat=3):
        degree=sum(abs(j) for j in k)
        if degree>N or sum(k)%3==0: continue
        if (k[0]-k[2])%N==0 and (k[1]-k[2])%N==0:
            low_bad.append(degree)
    assert min(low_bad)==N
    eta=np.exp(2j*np.pi/N)
    # Exact order-N invariant sum of the three pure Nth tensor powers.
    for i,j in itertools.product(range(N),repeat=2):
        assert np.allclose(np.array([eta**i,eta**j,eta**(-i-j)])**N,np.ones(3))
    return dict(first_fractional_invariant_degree=N,group_order=3*N*N)


def run():
    root=np.exp(2j*np.pi/3)
    so=so3_generators();su=su3_generators()
    groups={"SO3":ranks(so),"SU3":ranks(su),"mu3SO3":ranks(so,root)}
    assert [groups[n]["RR"] for n in ("SO3","SU3","mu3SO3")]==[1,0,0]
    assert all(x["RRR"]==x["RbarR"]==1 for x in groups.values())
    # Order: active u,d; untouched partner1,partner2.
    full=np.eye(9).reshape(9,9)/3
    d=np.eye(3).reshape(9)/np.sqrt(3)
    Pd=np.outer(d,d)
    pinv=invariant_projector(so,[1,1])
    assert np.allclose(pinv,Pd)
    yes=pinv@full
    accept=np.linalg.norm(yes)**2
    assert np.allclose(yes/np.sqrt(accept),Pd)
    common_center=np.eye(9)*root**2
    assert np.allclose(common_center@Pd@common_center.conj().T,Pd)
    assert np.allclose(common_center@d,root**2*d)
    true_neutral={n:groups[n]["RR"]/9 for n in groups}
    dark=(1-np.cos(2*np.pi/3))/2
    assert np.isclose(accept*dark,1/12)
    # A left fundamental and a right conjugate with the same electric residue.
    LR=invariant_projector(su,[1,-1])
    assert np.allclose(LR,Pd)
    # Example instrument budget: route contrast >= .96, total phase deviation <= .06 rad.
    floor=(1/9-.002)*(.75-(.04+.06)/2)-.002
    assert floor>.07
    return dict(full_invariant_ranks=groups,
                true_neutral_and_dark={n:float(s*dark) for n,s in true_neutral.items()},
                only_quadratic_projector_joint_for_mu3SO3=float(accept*dark),
                chiral_graph_cross_singlets=int(round(np.trace(LR).real)),
                finite_irreducible_examples=[valuation_test(n) for n in (4,5,7,10)],
                SO3_predeclared_error_floor=floor)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a: compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    else:
        assert np.isclose(a,b,atol=1e-12,rtol=1e-12)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    data=run()
    out=Path(__file__).with_name("independent_results.json")
    if args.check:
        compare(data,json.loads(out.read_text(encoding="utf-8")))
    else:
        with out.open("x",encoding="utf-8") as stream:
            json.dump(data,stream,indent=2,allow_nan=False);stream.write("\n")
    print("PASS",data["true_neutral_and_dark"],data["SO3_predeclared_error_floor"])
