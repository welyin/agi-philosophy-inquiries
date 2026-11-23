"""Two independent neutral pairs and one retained scalar/adjoint source."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def adj(a):
    return a.conj().T


def su(n):
    ops = []
    for j in range(1, n):
        for i in range(j):
            m = np.zeros((n, n), complex)
            m[i, j] = 1
            ops.extend([(m + adj(m))/2, (m - adj(m))/(2j)])
        ops.append(np.diag([1]*j + [-j] + [0]*(n-j-1))/math.sqrt(2*j*(j+1)))
    return ops


def so3():
    ops = []
    for i, j in ((1, 2), (2, 0), (0, 1)):
        m = np.zeros((3, 3), complex)
        m[i, j], m[j, i] = -1j, 1j
        ops.append(m)
    return ops


def check(a,b,tol=3e-11):
    np.testing.assert_allclose(a,b,atol=tol,rtol=tol)


def blocks(gen, partner=False, one_arm=False):
    n=len(gen[0])
    eye=np.eye(n)
    left=[np.kron(g,eye) for g in gen]
    right=[np.kron(eye,-g.conj() if partner else g) for g in gen]
    return left if one_arm else [a+b for a,b in zip(left,right)]


def pulse_matrix(ops):
    lift=np.vstack(ops)
    data=lift.shape[1]
    h=np.zeros((len(lift)+data,len(lift)+data),complex)
    h[data:,:data]=lift
    h[:data,data:]=adj(lift)
    return h


def casimir(gen):
    n=len(gen[0])
    c=sum(g@g for g in gen)
    value=np.trace(c).real/n
    check(c,value*np.eye(n))
    return value


def pulse_stats(vals, t):
    return float(np.mean(np.sin(t*np.sqrt(np.maximum(vals,0)))**2))


def open_group(label, gen):
    n=len(gen[0])
    d=len(gen)
    c=casimir(gen)
    eye=np.eye(n*n)
    ops=blocks(gen)
    unit=[v/math.sqrt(c) for v in ops]
    l=sum(v@v for v in unit)
    lv=np.linalg.eigvalsh(l)
    check(np.trace(l)/(n*n),2.)
    m2=np.trace(l@l).real/(n*n)
    check(m2,4+4/d)
    singlet=np.eye(n).ravel()/np.sqrt(n)
    one=blocks(gen,one_arm=True)
    conjugate=blocks(gen,partner=True)
    check(np.linalg.norm(np.vstack(conjugate)@singlet),0)
    check(sum(adj(a)@a for a in one),c*eye)
    h=pulse_matrix(unit)
    energy,basis=np.linalg.eigh(h)
    # Each column is a reference-pair label. Empty initial source, equal
    # Schmidt weights give the two actual neutral pairs after mode regrouping.
    state=np.zeros((len(h),n*n),complex)
    state[:n*n]=eye/n
    initial=adj(basis)@state
    probs={}
    for t in (.1,.2,math.pi/2,math.pi):
        final=basis@(np.exp(-1j*energy*t)[:,None]*initial)
        seen=np.linalg.norm(final[n*n:])**2
        missed=np.linalg.norm(final[:n*n])**2
        check(seen+missed,1)
        # All original partners retained, completeness holds on arbitrary
        # attached spectators because the full pulse and records are isometric.
        check(adj(final)@final,eye/(n*n))
        check(seen,pulse_stats(lv,t))
        probs[str(t)]=float(seen)
    for a,t in enumerate(gen):
        # adjoint index convention O_ab in U T_b U† = sum_a O_ab T_a
        source=np.zeros((1+d,1+d),complex)
        for b,v in enumerate(gen):
            bracket=t@v-v@t
            for k,w in enumerate(gen):
                source[k+1,b+1]=np.trace(w@bracket)/np.trace(w@w)
        check(source,adj(source))
        total=np.kron(np.eye(1+d),ops[a])+np.kron(source,eye)
        check(total@h,h@total)
    expect=([(-.0,1),(1,3),(3,5)] if label=="SO3"
            else [(1,3),(2.5,6)])
    check(lv,np.repeat(*zip(*expect)))
    short=(4*probs["0.1"]-probs["0.2"])/(4*.1**4)
    return dict(n=n,adjoint_dimension=d,one_arm_casimir=c,second_moment=float(m2),
                two_time_estimate=short,pulse_probabilities=probs,
                normalized_total_casimir=[list(k) for k in expect])


def extra_groups():
    extra=[]
    for n in (2,4):
        g=su(n)
        vals=np.linalg.eigvalsh(sum(a@a for a in blocks(g))/casimir(g))
        check(np.mean(vals),2.)
        check(np.mean(vals**2),4+4/(n*n-1))
        extra.append(float(np.mean(vals**2)))
    # Irreducible tensor product of SU(2) doublet and SO(3) vector.
    a=[np.kron(x,np.eye(3)) for x in su(2)]
    b=[np.kron(np.eye(2),1.6*x) for x in so3()]
    p=casimir(a)/(casimir(a)+casimir(b))
    items=blocks(a)+blocks(b)
    l=sum(x@x for x in items)/(casimir(a)+casimir(b))
    expected=4+4*(p*p/3+(1-p)**2/3)
    check(np.trace(l@l)/len(l),expected)
    return dict(simple_second_moments=extra,first_weight=p,semisimple_second_moment=expected)


def result():
    so=open_group("SO3",so3())
    su3=open_group("SU3",su(3))
    gap=su3["pulse_probabilities"][str(math.pi)]-so["pulse_probabilities"][str(math.pi)]
    residual=(34/45)*(64/3)*(.1**2)
    tol=(5e-5)/(4*.1**4)
    assert 2*(residual+tol)<5/6
    for model in (so,su3):
        assert abs(model["two_time_estimate"]-model["second_moment"])<residual
    assert gap>.31
    return dict(round=1343,SO3=so,SU3=su3,other=extra_groups(),
                finite_window=dict(moment_model_gap=5/6,bias_bound=residual,
                                   input_record_error=1e-5,moment_error_from_two_reads=tol,
                                   shot_count_per_time=math.ceil(math.log(400)/(2*(5e-6)**2)),
                                   pi_absolute_gap=gap,pi_gap_after_two_005_errors=gap-.1))


if __name__=="__main__":
    arg=argparse.ArgumentParser()
    arg.add_argument("--write",action="store_true")
    opts=arg.parse_args()
    data=result()
    path=Path(__file__).with_name("results.json")
    if opts.write:
        with path.open("x",encoding="utf-8") as f:
            json.dump(data,f,indent=2)
            f.write("\n")
    else:
        def match(x,y):
            if isinstance(x,dict):
                assert x.keys()==y.keys()
                for k in x: match(x[k],y[k])
            else: check(x,y)
        match(data,json.loads(path.read_text(encoding="utf-8")))
    print(json.dumps(data["finite_window"]))

