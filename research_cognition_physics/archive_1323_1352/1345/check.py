"""Round 1345. Complete four-route effects, including retained neutral partners."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def adj(a):
    return a.conj().T


def check(a,b,tol=5e-10):
    np.testing.assert_allclose(a,b,atol=tol,rtol=tol)


def su3():
    g=[]
    for j in (1,2):
        for i in range(j):
            m=np.zeros((3,3),complex)
            m[i,j]=1
            g.extend([(m+adj(m))/2,(m-adj(m))/(2j)])
        g.append(np.diag([1]*j+[-j]+[0]*(2-j))/math.sqrt(2*j*(j+1)))
    return g


def so3():
    g=[]
    for i,j in ((1,2),(2,0),(0,1)):
        m=np.zeros((3,3),complex)
        m[i,j],m[j,i]=-1j,1j
        g.append(m)
    return g


def rms(ops):
    return math.sqrt(sum(np.vdot(a,a).real for a in ops)/len(ops[0]))


def swap(n):
    return np.eye(n*n).reshape(n,n,n,n).swapaxes(0,1).reshape(n*n,n*n)


def wedge(n):
    assert n==3
    w=np.zeros((9,3),complex)
    for k,(i,j) in enumerate(((1,2),(2,0),(0,1))):
        w[3*i+j,k]=1/math.sqrt(2)
        w[3*j+i,k]=-1/math.sqrt(2)
    return w


def pair_on_four(q, legs):
    rest=[j for j in range(4) if j not in legs]
    inv=np.argsort(list(legs)+rest)
    axes=list(inv)+list(inv+4)
    return np.kron(q,np.eye(9)).reshape((3,)*8).transpose(axes).reshape(81,81)


def routes(q):
    left=[pair_on_four(x,(0,1))+pair_on_four(x,(2,3)) for x in q]
    right=[pair_on_four(x,(0,2))+pair_on_four(x,(1,3)) for x in q]
    return left,right


def partial(x):
    return np.trace(x.reshape(3,3,3,3),axis1=1,axis2=3)/3


def additive(a):
    return np.kron(a,np.eye(3))+np.kron(np.eye(3),a)


def responses(q,t):
    lifted=np.vstack(q)
    vals,v=np.linalg.eigh(adj(lifted)@lifted)
    lam=np.sqrt(np.maximum(vals,0))
    ground=(v*np.cos(t*lam))@adj(v)
    exc=-1j*lifted@((v*(t*np.sinc(t*lam/np.pi)))@adj(v))
    check(adj(ground)@ground+adj(exc)@exc,np.eye(len(lam)))
    return ground,exc


def four_records(q,t):
    left,right=routes(q)
    lr=responses(left,t)
    rr=responses(right,t)
    out={}
    check_full=np.zeros((81,81),complex)
    for sign,name in ((1,"plus"),(-1,"minus")):
        for a,b,code in zip(lr,rr,("ground","excited")):
            kraus=(a+sign*b)/2
            check_full+=adj(kraus)@kraus
            out[name+"_"+code]=float(np.vdot(kraus,kraus).real/81)
    check(check_full,np.eye(81))
    return out,rms([a-b for a,b in zip(left,right)])


def groups(g):
    c=sum(x@x for x in g)[0,0].real
    t=[x/math.sqrt(c) for x in g]
    canonical=[additive(x) for x in t]
    w=wedge(3)
    anti=w@adj(w)
    check(anti,(np.eye(9)-swap(3))/2)
    for a,x in zip(canonical,t):
        check(adj(w)@a@w,-x.T)
    rows={}
    for gain in (0.,.3,1.,1.4):
        q=[v@(anti+gain*(np.eye(9)-anti)) for v in canonical]
        for x,y in zip(q,t):
            check(adj(w)@x@w,-y.T)
        link=[partial(x) for x in q]
        z=rms([x-additive(y) for x,y in zip(q,link)])
        phi=.25
        four,defect=four_records(q,phi)
        check(defect,2*z)
        rms_error=rms([x-y for x,y in zip(q,canonical)])
        assert rms_error<=math.sqrt(7)*z+2e-13
        if gain==0:
            check(z**2,5/18)
            check_four=responses(q,math.pi/2)
            check(np.vdot(check_four[1],check_four[1]).real/9,1/3)
            for a,b in zip(link,t): check(a,b/6)
            assert four["minus_excited"] >= (phi*math.sqrt(5/18)-(2*phi)**3/6)**2
        if gain==1:
            check(four["minus_excited"],0.)
            check(four["minus_ground"],0.)
        rows[str(gain)]=dict(z=z,defect=defect,rms_from_canonical=rms_error,records=four)
    return rows


def generic_exchange(seed):
    # Does not select an irreducible-color basis. Verifies the entire
    # primitive decomposition for arbitrary swap symmetric pair operators.
    rng=np.random.default_rng(seed)
    sw=swap(3)
    inv=np.eye(9)
    worst=0.
    w=wedge(3)
    for _ in range(12):
        m=rng.normal(size=(9,9))+1j*rng.normal(size=(9,9))
        x=m+adj(m)
        x=(x+sw@x@sw)/2
        x-=inv*np.trace(x)/9
        a=partial(x)
        corr=x-additive(a)
        l,r=routes([x])
        diff=abs(rms([l[0]-r[0]])-2*rms([corr]))
        worst=max(worst,diff)
        check(diff,0.)
        # dual volume action of a general Hermitian matrix
        check(adj(w)@additive(a)@w,np.trace(a)*np.eye(3)-a.T)
    return worst


def run():
    angle=.25
    low=(angle*math.sqrt(5/18)-(2*angle)**3/6)**2
    assert low>.0123
    norms=[]
    for p in (1e-8,1e-6):
        # Only an illustration of the propagated finite witness, not a measured p.
        s=(3*math.sqrt(p)/8)**(1/3)
        norms.append(dict(t=s,upper_p=p,z_bound=math.sqrt(p)/s+2**3*s*s/6))
    return dict(round=1345,SU3=groups(su3()),SO3=groups(so3()),
                finite_theta=angle,competitor_probability_lower=low,
                generic_fourleg_max_error=generic_exchange(1345),fourleg_examples=norms)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    opts=parser.parse_args()
    path=Path(__file__).with_name("results.json")
    data=run()
    if opts.write:
        with path.open("x",encoding="utf-8") as f:
            json.dump(data,f,indent=2)
            f.write("\n")
    else:
        def match(x,y):
            if isinstance(x,dict):
                assert x.keys()==y.keys()
                for k in x: match(x[k],y[k])
            elif isinstance(x,list):
                for a,b in zip(x,y,strict=True): match(a,b)
            else: check(x,y)
        match(data,json.loads(path.read_text(encoding="utf-8")))
    print(json.dumps({"floor":data["competitor_probability_lower"],
                      "exact_su":data["SU3"]["0.0"]["records"]["minus_excited"],
                      "exact_so":data["SO3"]["0.0"]["records"]["minus_excited"]}))

