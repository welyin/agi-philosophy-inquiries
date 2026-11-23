"""Independent resultant/Sturm count for the third-family spurion and the electron holdout."""
import argparse
import cmath
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def texture(rows):
    down=np.array(rows[:3]);t12,t13,t23,delta=rows[3:7]
    def rotate(i,j,a,phase=0.):
        v=np.eye(3,dtype=complex)
        v[i,i]=v[j,j]=math.cos(a)
        v[i,j]=math.sin(a)*cmath.exp(-1j*phase)
        v[j,i]=-v[i,j].conjugate()
        return v
    v=rotate(1,2,t23)@rotate(0,2,t13,delta)@rotate(0,1,t12)
    h=(v*(down/down[-1])**2)@v.conj().T
    left=np.diag([1,cmath.exp(1j*cmath.phase(h[0,1])),cmath.exp(1j*cmath.phase(h[0,2]))])
    h=left@h@left.conj().T
    A=math.sqrt(h[0,0].real)
    C,E=(h[0,1:]/A).real
    ff=math.sqrt(h[2,2].real-E**2)
    z=(h[1,2]-C*E)/ff
    B=-float(np.prod(down/down[-1]))/(A*ff)
    m=np.array([[0,A,0],[B,C,z],[0,E,ff]])
    np.testing.assert_allclose(m@m.conj().T,h,atol=3e-14)
    return [A,B,C,E,ff,z.real,z.imag]


def mtx(e,cl,t):
    A,B,C,E,ff,zr,zi=e
    ca,cb,cc=cl
    return np.array([[0,cb*A,0],[ca*B,cc*C,t*complex(zr,zi)],[0,t*E,t*ff]])


def char_resultant(e,cl,r):
    A,B,C,E,f,zr,zi=(F(f"{val:.12f}") for val in e)
    ca,cb,cc=map(lambda x:F(str(x)),cl)
    A*=cb;B*=ca;C*=cc
    z=zr*zr+zi*zi
    def p(*values): return Polynomial([F(x) for x in values])
    a=p(A*A+B*B+C*C,0,z+E*E+f*f)
    b=p(A*A*B*B,0,A*A*(f*f+z)+B*B*(E*E+f*f)+C*C*f*f,-2*C*f*E*zr,E*E*z)
    c=p(0,0,A*A*B*B*f*f)
    # For remaining eigenvalue x, the other roots y,z obey (y-r*z)*(z-r*y)=g(x).
    d=1+r+r*r
    h=-a*(1+r*r)
    j=b*(1+r)**2-r*a*a
    # det g(Q), Q the companion of x³-a x²+b x-c.
    mm=[[j,d*c,(h+d*a)*c],[h,j-d*b,-h*b+d*(c-a*b)],
        [p(d),h+d*a,j+h*a+d*(a*a-b)]]
    D=mm[0][0]*(mm[1][1]*mm[2][2]-mm[2][1]*mm[1][2])
    D-=mm[0][1]*(mm[1][0]*mm[2][2]-mm[2][0]*mm[1][2])
    D+=mm[0][2]*(mm[1][0]*mm[2][1]-mm[2][0]*mm[1][1])
    return list(D.coef)


def primitive(p):
    while p and p[-1]==0:p.pop()
    if not p:return []
    den=math.lcm(*(c.denominator for c in p))
    nums=[int(x*den) for x in p]
    gcd=math.gcd(*nums)
    return [x//gcd for x in nums]


def sturm(coeff):
    p=primitive(coeff)
    chain=[p,primitive([i*p[i] for i in range(1,len(p))])]
    while len(chain[-1])>1:
        rem=list(chain[-2]);div=chain[-1]
        while len(rem)>=len(div):
            k=len(rem)-len(div);s=F(rem[-1],div[-1])
            for i,d in enumerate(div):rem[i+k]-=s*d
            while rem and rem[-1]==0:rem.pop()
        rem=primitive([-x for x in rem])
        assert rem,"separate repeated roots before sign counting"
        chain.append(rem)
    return chain


def sign_at(p,x):
    den=x.denominator;n=x.numerator
    b=p[-1];power=1
    for q in p[-2::-1]:
        power*=den;b=b*n+q*power
    return (b>0)-(b<0)


def all_roots(coeff):
    chain=sturm(coeff)
    def changes(signs):
        s=[x for x in signs if x]
        return sum(a!=b for a,b in zip(s,s[1:]))
    def variations(x):
        return changes([sign_at(p,x) for p in chain])
    left=changes([(1 if p[-1]>0 else -1)*(-1)**(len(p)-1) for p in chain])
    right=changes([1 if p[-1]>0 else -1 for p in chain])
    radius=F(1)
    while variations(-radius)!=left or variations(radius)!=right:
        radius*=2
    pending=[(-radius,radius,left,right)]
    out=[]
    while pending:
        a,b,va,vb=pending.pop()
        if va==vb:continue
        if b-a<F(1,10**12) and va-vb==1:
            out.append(float((a+b)/2));continue
        mid=(a+b)/2
        assert sign_at(chain[0],mid)!=0
        d=variations(mid)
        pending.extend([(a,mid,va,d),(mid,b,d,vb)])
    assert len(out)==left-right
    return sorted(out)


def point_solution(edges,cl,t,mu,tau):
    m=mtx(edges,cl,t)
    sig=np.linalg.svd(m,compute_uv=False)[::-1]
    scale=tau/sig[-1]
    # Independent electron readout by determinant using only the fitted heavy pair.
    pred=scale**3*abs(np.linalg.det(m))/(mu*tau)
    return sig,scale,pred


def linked_root(rows,cl):
    edges=texture(rows)
    mu,tau=rows[7:9]
    def deviation(t):
        sig=np.linalg.svd(mtx(edges,cl,t),compute_uv=False)
        return sig[0]/sig[1]-tau/mu
    lo,hi=.3,2.3
    a,b=deviation(lo),deviation(hi)
    if a*b>=0:
        return None
    for _ in range(47):
        mid=(lo+hi)/2
        if deviation(mid)*a>0:lo=mid
        else:hi=mid
    t=(lo+hi)/2
    _,k,pred=point_solution(edges,cl,t,mu,tau)
    return t,k,pred


def solve():
    quark=read("../1397/inputs.json")["high_scale"]
    leptons=read("../1402/inputs.json")["scales"]["1e16"][3:]
    obs=leptons[0][0]
    bounds=quark["down"]+quark["angles_theta12_theta13_theta23_delta"]+leptons[1:]
    mids,errs=np.transpose(bounds)
    edges=texture(mids)
    e12=[float(f"{x:.12f}") for x in edges]
    ratio=(F(str(mids[7]))/F(str(mids[8])))**2
    table={}
    for key,cl in read("../1486/inputs.json")["clebsch"].items():
        roots=all_roots(char_resultant(edges,cl,ratio))
        paircounts=[0]*3;hits=[]
        for t in roots:
            y,k,pred=point_solution(e12,cl,t,*mids[7:])
            alt=np.array([y[1]/y[2],y[0]/y[1],y[0]/y[2]])**2
            rel=abs(alt/float(ratio)-1)
            i=np.argmin(rel);paircounts[i]+=1
            assert rel[i]<1e-6 and sorted(rel)[1]>.001
            if i==0:
                assert abs(k*y[0]/pred-1)<1e-7
                hits.append([t,pred/obs])
        light=np.linalg.svd(mtx(edges,cl,1)[:2,:2],compute_uv=False)
        table[key]={"Sturm_count_all_pairs":len(roots),"mu_tau_e_mu_e_tau_counts":paircounts,
                    "mu_tau_roots_and_e_over_measured":hits,
                    "drop_third_then_calibrate_mu_e_over_measured":mids[7]*light[1]/light[0]/obs}
    # A constructive test for whether the best central branch really fails the whole published box.
    GJ=[1,1,3]
    best=(0,None)
    for v in itertools.product((-1,1),repeat=len(mids)):
        ans=linked_root(mids+2*errs*np.array(v),GJ)
        if ans and ans[-1]>best[0]:best=(ans[-1],v)
    perturb=np.array(best[1])*2*errs
    h0=linked_root(mids,GJ)
    example={"largest_tested_corner_over_central_e":best[0]/obs,"direction_pm":list(best[1])}
    q_only=mids.copy();q_only[:7]+=perturb[:7]
    qt,qk,qe=linked_root(q_only,GJ)
    assert abs(qe-obs)<2*leptons[0][1]
    example["quark_corner_fixed_mu_tau"]={"t":qt,"k_times_yb":qk,"predicted_e":qe,
                                         "electron_marginal_pull":(qe-obs)/leptons[0][1]}
    if (best[0]-obs)*(h0[-1]-obs)<0:
        a,b=0.,1.
        for _ in range(42):
            mid=(a+b)/2
            if linked_root(mids+mid*perturb,GJ)[-1]<obs:a=mid
            else:b=mid
        vals=mids+((a+b)/2)*perturb
        t,k,pred=linked_root(vals,GJ)
        assert abs(pred/obs-1)<1e-8
        example.update({"max_abs_input_pull":a+b,"marginal_coordinates":vals.tolist(),
                        "t":t,"k_times_yb":k,"predicted_e":pred,
                        "warning":"experimental e selected interpolation; this point proves joint feasibility only"})
    return {"scale_GeV":1e16,"central_certificates":table,"extended_box_GJ_test":example}


def stable(x):
    if isinstance(x,dict):return {k:stable(v) for k,v in x.items()}
    if isinstance(x,list):return [stable(v) for v in x]
    if isinstance(x,str):return x
    return float(f"{x:.11g}")


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    results=stable(solve())
    path=HERE/"independent_results.json"
    if args.check:
        assert json.loads(path.read_bytes())==results
        print("1657 independent PASS: resultant/Sturm isolation, physical pair ordering, electron and joint input test")
    else:
        with path.open("x",encoding="utf8") as out:
            out.write(json.dumps(results,indent=2)+"\n")

