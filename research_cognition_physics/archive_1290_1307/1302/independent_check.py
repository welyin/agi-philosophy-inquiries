"""Round 1302: independent binomial-transform, anomaly and matrix diagnostics.
No author code is imported. The all-N proof is analytic; finite matrices diagnose it.
Default/--check read-only. --write exclusively creates independent_results.json.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from math import comb
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"


def n_from_flow(x):
    return [a+b for a,b in zip((0,)+tuple(x),tuple(x)+(0,))]


def moments(x):
    return [sum((-1)**i*i**r*a for i,a in enumerate(x)) for r in range(3)]


def shift_coefficients(x):
    # X(z) = sum_j b_j (z+1)^j; integer binomial translation, no division algorithm.
    b = [sum(x[i]*comb(i,j)*(-1)**(i-j) for i in range(j,len(x)))
         for j in range(len(x))]
    rebuilt = [sum(b[j]*comb(j,i) for j in range(i,len(b))) for i in range(len(b))]
    assert rebuilt == list(x)
    assert sum(x) == sum(a*2**j for j,a in enumerate(b))
    m0,m1,m2 = moments(x)
    assert b[0] == m0
    if len(b)>1: assert b[1] == -m1
    if len(b)>2: assert 2*b[2] == m2-m1
    return b


def anomaly_polynomials(x):
    """Exact ascending k coefficients computed directly from field charges."""
    n = n_from_flow(x)
    polynomials = {key:[0]*4 for key in ("g","g3","F","F3","gF2","g2F")}
    for i,count in enumerate(n):
        fi = (-1)**i
        sign = -fi
        for key,power,fpow in (("g",1,0),("g3",3,0),("F",0,1),
                               ("F3",0,3),("gF2",1,2),("g2F",2,1)):
            for j in range(power+1):
                polynomials[key][j] += count*sign**power*fi**fpow*comb(power,j)*i**(power-j)
    m0,m1,m2 = moments(x)
    expected = {
        "g":[m0,0,0,0],
        "g3":[m0+3*m1+3*m2,3*m0+6*m1,3*m0,0],
        "F":[0,0,0,0],"F3":[0,0,0,0],
        "gF2":[m0,0,0,0],
        "g2F":[-m0-2*m1,-2*m0,0,0]}
    assert polynomials == expected
    return polynomials


def evaluate(coeff,k):
    return sum(c*k**i for i,c in enumerate(coeff))


def algebra_certificates():
    cases = [
        ("linear_not_cancelled",(1,),1,False),
        ("old_N12",(2,2,0,1,1),2,False),
        ("old_N20",(3,5,2),1,False),
        ("N16_binomial",(1,3,3,1),1,True),
        ("N48_binomial",(3,9,9,3),7,True),
        # A different positive flow with a quotient containing a negative term.
        # X=(1+z)^3*(1-z+z^2); quotient positivity is not assumed.
        ("N16_nonunique_flow",(1,2,1,1,2,1),5,True)]
    rows = []
    for name,x,k,protected in cases:
        n = n_from_flow(x)
        assert x[0]>0 and x[-1]>0 and all(t>0 for t in n)
        b = shift_coefficients(x)
        poly = anomaly_polynomials(x)
        values = {key:evaluate(coeff,k) for key,coeff in poly.items()}
        if protected:
            assert all(value==0 for value in values.values())
            assert b[:3] == [0,0,0]
            c = sum(b[j]*2**(j-3) for j in range(3,len(b)))
            assert isinstance(c,int) and c>0 and sum(n)==16*c
        else:
            c = None
        if name in ("old_N12","old_N20"):
            assert values["g"]==values["g3"]==0 and values["g2F"]==2
            assert all(values[key]==0 for key in ("F","F3","gF2"))
        rows.append({"name":name,"flow":list(x),"k":k,"N":sum(n),"moments":moments(x),
                     "X_in_powers_of_z_plus_1":b,"local_anomaly_polynomials":poly,
                     "local_coefficients_at_k":values,"new_protection":protected,
                     "positive_multiple_c_when_protected":c})
    return rows


def matching_matrix(flow,k):
    n=n_from_flow(flow)
    roles=[(-1)**(i+1)*(k+i) for i in range(len(n))]
    charges=[q for q,count in zip(roles,n) for _ in range(count)]
    identity=[(-1)**i for i,count in enumerate(n) for _ in range(count)]
    groups=[];offset=0
    for count in n:
        groups.append(list(range(offset,offset+count)));offset+=count
    remaining=[g.copy() for g in groups]
    N=len(charges)
    base=np.zeros((N,N),complex)
    for edge,count in enumerate(flow):
        for _ in range(count):
            a,b=remaining[edge].pop(),remaining[edge+1].pop()
            base[a,b]=base[b,a]=1
    assert not any(remaining)
    assert np.array_equal(base@base,np.eye(N))
    q=np.array(charges)
    f=np.array(identity)
    perturbation=np.zeros_like(base)
    for a in range(N):
        for b in range(a+1,N):
            if abs(charges[a]+charges[b])==1:
                value=((1+(a+2*b)%5)+1j*((2*a+b)%3-1))/8
                assert abs(value)<1
                perturbation[a,b]=perturbation[b,a]=value
    actual=base+perturbation/(16*N)
    assert np.linalg.norm(perturbation,"fro")<N
    assert np.array_equal(actual,actual.T)
    assert np.array_equal((f[:,None]+f[None,:])*actual,np.zeros_like(actual))
    # Role-constant nonzero F is unique up to scale, since this connected path
    # forces f_{i+1}=-f_i. Check the edge incidence rank independently.
    incidence=np.zeros((len(n)-1,len(n)),int)
    for i in range(len(n)-1):incidence[i,i:i+2]=1
    assert np.linalg.matrix_rank(incidence)==len(n)-1
    assert np.array_equal(incidence@np.array([(-1)**i for i in range(len(n))]),
                          np.zeros(len(n)-1))
    allowed=(np.abs(q[:,None]+q[None,:])==1)
    assert np.array_equal(actual!=0,allowed)
    seen={0}
    while True:
        extended=seen|{j for i in seen for j in range(N) if actual[i,j]!=0}
        if extended==seen:break
        seen=extended
    assert len(seen)==N
    singular=np.linalg.svd(actual,compute_uv=False)
    assert min(singular)>=15/16-1e-12
    alpha=.419
    phases=np.exp(1j*alpha*f)
    assert np.linalg.norm(phases[:,None]*actual*phases[None,:]-actual)<1e-12
    assert np.linalg.norm(np.diag(f)@actual-actual@np.diag(f))>1
    assert all(v==0 for v in (sum(q),sum(q**3),sum(f),sum(f**3),sum(q*f**2),sum(q**2*f)))
    return {"N":N,"k":k,"flow":list(flow),"all_six_local_coefficients_zero":True,
            "role_constant_F_nullity":1,"all_allowed_couplings_nonzero":True,
            "actual_field_support_connected":True,
            "mass_identity_exact_F_transpose_M_plus_MF_zero":True,
            "ordinary_commutator_nonzero_as_required":True,
            "analytic_smallest_singular_value_bound":"15/16",
            "numerical_smallest_singular_value":round(float(min(singular)),12),
            "no_actual_instrument_certified":True}


# Tiny exact polynomial differential-form algebra, with integer coefficients.
ZERO=(0,0,0,0)
def padd(p,q):
    r=p.copy()
    for mon,c in q.items():r[mon]=r.get(mon,0)+c
    return {mon:c for mon,c in r.items() if c}
def pscale(p,c):return {mon:a*c for mon,a in p.items() if a*c}
def pmul(p,q):
    r={}
    for a,ca in p.items():
        for b,cb in q.items():
            mon=tuple(u+v for u,v in zip(a,b));r[mon]=r.get(mon,0)+ca*cb
    return {m:c for m,c in r.items() if c}
def pdiff(p,i):
    r={}
    for mon,c in p.items():
        if mon[i]:
            new=list(mon);new[i]-=1;r[tuple(new)]=c*mon[i]
    return r
def fadd(a,b,scale=1):
    out={key:p.copy() for key,p in a.items()}
    for key,p in b.items():
        out[key]=padd(out.get(key,{}),pscale(p,scale))
        if not out[key]:del out[key]
    return out
def wedge(a,b):
    out={}
    for ia,pa in a.items():
        for ib,pb in b.items():
            sequence=ia+ib
            if len(set(sequence))!=len(sequence):continue
            inv=sum(sequence[i]>sequence[j] for i in range(len(sequence)) for j in range(i+1,len(sequence)))
            key=tuple(sorted(sequence))
            out[key]=padd(out.get(key,{}),pscale(pmul(pa,pb),(-1)**inv))
            if not out[key]:del out[key]
    return out
def exterior(a):
    out={}
    for indices,p in a.items():
        for i in range(4):
            dp=pdiff(p,i)
            if dp:out=fadd(out,wedge({(i,):dp},{indices:{ZERO:1}}))
    return out


def compensation_patch_check():
    coords=[]
    for i in range(4):
        mon=[0]*4;mon[i]=1;coords.append({tuple(mon):1})
    x,y,z,t=coords
    A={(1,):x,(3,):z}  # A=x dy+z dt, G^G is nonzero.
    theta={():pmul(x,z)}
    lam={():pmul(y,t)}
    alpha={():padd(pmul(x,x),pmul(y,z))}
    C=fadd(A,exterior(theta),-1)
    G=exterior(A)
    transformed=fadd(fadd(A,exterior(lam)),exterior(fadd(theta,lam)),-1)
    assert transformed==C
    assert exterior(C)==G
    GG=wedge(G,G)
    assert GG=={(0,1,2,3):{ZERO:2}}
    current=wedge(C,G)
    assert exterior(current)==GG
    variation=wedge(exterior(alpha),current)
    total_derivative=exterior(wedge(alpha,current))
    # On compact-supported variations: integral variation = -integral alpha G^G.
    assert fadd(total_derivative,fadd(variation,wedge(alpha,GG)),-1)=={}
    return {"local_gauge_invariant_A_minus_dtheta":True,
            "dC_equals_G_exact":True,"G_wedge_G_coefficient":2,
            "d_current_shift_equals_G_wedge_G_exact":True,
            "counterterm_variation_mod_boundary":"-t * integral(alpha G wedge G)",
            "original_fermion_anomaly_coefficient_changed":False,
            "source_defined_total_identity_current_changed":True,
            "S_nonzero_patch_only":True,
            "regular_at_S_zero_or_global_patchability_proved_by_code":False}


def payload():
    sources=["../_shared/round1302_admission.md","../1301/proof.md","../1301/acceptance.json"]
    return {"status":"PASS_INDEPENDENT_1302_ALGEBRA_MATRIX_AND_PATCH_DIAGNOSTICS",
            "exact_binomial_transform_and_anomalies":algebra_certificates(),
            "finite_common_mass_objects":[matching_matrix((1,3,3,1),1),
                                         matching_matrix((3,9,9,3),7)],
            "Higgs_patch_compensation":compensation_patch_check(),
            "all_N_claim":"analytic: N=16c for every positive integer c in stated class",
            "extra_cognitive_axioms_adopted":0,"full_global_anomaly_or_O_claim":False,
            "source_bindings":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in sources},
            "independent_check_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group()
    flags.add_argument("--write",action="store_true")
    flags.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=payload()
    if args.write:
        with OUT.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2);stream.write("\n")
        print("CREATED",result["status"])
    else:
        if json.loads(OUT.read_text(encoding="utf-8"))!=result:
            raise SystemExit("MISMATCH: independent diagnostics or bindings changed")
        print("READ_ONLY_MATCH",result["status"])

if __name__=="__main__":
    main()

