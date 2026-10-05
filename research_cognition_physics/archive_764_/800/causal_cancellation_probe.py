"""Working 800: closed causal cancellation in one noncommutative BV algebra."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import importlib.util
import json
import sys
sys.dont_write_bytecode=True
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'798'))
import internal_register_probe as base
C,cmat,scale,dag,same=base.C,base.matrix,base.scale,base.dag,base.same
spec=importlib.util.spec_from_file_location('quartet800',HERE.parent/'784/local_charge_boundary.py')
quartet=importlib.util.module_from_spec(spec);spec.loader.exec_module(quartet)
m5,par5,q5,_,p5=quartet.quartet()
def to_c(a): return cmat(a.tolist())
I=base.eye(10);Z=scale(I,0)
metric=to_c(np.kron(m5,np.eye(2,dtype=object)))
charge=to_c(np.kron(q5,np.eye(2,dtype=object)))
parity=to_c(np.kron(par5,np.eye(2,dtype=object)))
sharp=lambda a:metric@dag(a)@metric
delta=lambda a,odd=False:charge@a-(-1 if odd else 1)*a@charge
def clean(s): return {k:a for k,a in s.items() if any(a.flat)}
def add(*args):
    o={}
    for s in args:
        for k,a in s.items():o[k]=o.get(k,Z)+a
    return clean(o)
def sc(s,v):return clean({k:scale(a,v) for k,a in s.items()})
def mul(s,t):
    out={}
    for k,a in s.items():
        for l,b in t.items():
            key=(k[0]+l[0],k[1]+l[1])
            if key[0]<=3 and key[1]<=2:out[key]=out.get(key,Z)+a@b
    return clean(out)
ONE={(0,0):I}
def adj(s):return {k:sharp(a) for k,a in s.items()}
def equal(s,t):return not add(s,sc(t,-1))
def exp(g):
    out,power=ONE,ONE
    for n in range(1,6):
        power=mul(power,g)
        out=add(out,sc(power,Q(1,factorial(n))))
    return out
def gen(mat,source=False):return {(0,1) if source else (1,0):scale(mat,C(0,1))}
def first(s):
    k=min(s)
    i,j=next((i,j) for i in range(10) for j in range(10) if s[k][i,j])
    return dict(epsilon=k[0],source=k[1],row=i,column=j,value=str(s[k][i,j]))

def run():
    def h(x,c):
        k=to_c(np.kron(quartet.elementary(5,0,4)+quartet.elementary(5,4,3),x))
        ex=delta(k,True);ex=ex+sharp(ex)
        return to_c(np.kron(np.eye(5,dtype=object),x))+scale(ex,c)
    hp=h(base.X,Q(1,5));hm=h(base.Z,Q(1,7))
    hf=h(base.X+base.Z,Q(1,11));obs=h(base.Y,Q(1,13))
    for x in (hp,hm,hf,obs):
        assert same(sharp(x),x) and not any(delta(x).flat)
        assert same(x@parity,parity@x)
    past,middle,future,source=[exp(gen(x,s)) for x,s in
                               ((hp,False),(hm,False),(hf,False),(obs,True))]
    original=mul(mul(future,middle),past)
    cut=mul(future,past)
    cut_source=mul(mul(future,source),past)
    rb=mul(adj(original),cut)
    rb_source=mul(adj(original),cut_source)
    double=mul(adj(rb),rb_source)
    common=mul(mul(adj(past),source),past)
    assert equal(double,common)
    assert equal(mul(adj(rb),rb),ONE)
    for family in (rb,rb_source,double,common):
        assert all(not any(delta(x).flat) for x in family.values())
    # Original late source retained, then represented with the same past dictionary.
    original_source=mul(mul(mul(future,source),middle),past)
    original_relative=mul(adj(original),original_source)
    shifted=mul(mul(adj(middle),source),middle)
    reconstructed=mul(mul(adj(past),shifted),past)
    assert equal(original_relative,reconstructed)
    naive=add(rb_source,sc(common,-1))
    assert naive
    # Replacing a completed compensator by a nonclosed one destroys Ward.
    raw=h(base.X,Q(1,6))+to_c(np.kron(quartet.elementary(5,2,4),base.I))
    raw=raw+sharp(raw)
    bad_compensator=exp(gen(raw))
    bad={k:delta(x) for k,x in bad_compensator.items() if any(delta(x).flat)}
    assert bad
    return dict(working_round=800,all_checks_passed=True,
                epsilon_order=3,source_order=2,full_BV_matrix_dimension=10,
                double_relative_cancellation_exact=True,
                same_past_unitary_independent_of_source=True,
                full_Ward_preserved=True,
                original_late_source_reconstructed_without_replacing_it=True,
                omission_of_compensator_inverse_defect=first(naive),
                uncompleted_compensator_Ward_defect=first(bad),
                original_continuum_zero_strip_completion_proven=False,
                physical_time_slice_full_theorem_signed=False,
                formal_round_completed=False,new_numbered_test_groups=0)
if __name__=='__main__':
    out=run()
    (HERE/'causal_cancellation_probe_results.json').write_text(
        json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,ensure_ascii=False,indent=2))
