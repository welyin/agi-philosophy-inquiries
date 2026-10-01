"""628: original quotient/matter bundle map and global anomaly interface.

Exact polynomial and bundle diagnostics only. Spin bordism and index theorems
are mathematical inputs, not conclusions of this finite computation.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_spinor_subgroup_mass as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_global_anomaly_bundle_results.json'


def clean(p):
    return {k:Fraction(v) for k,v in p.items() if v}


def add(*ps):
    out={}
    for p in ps:
        for k,v in p.items():out[k]=out.get(k,0)+v
    return clean(out)


def scale(p,a):
    return clean({k:a*v for k,v in p.items()})


def mul(p,q):
    out={}
    for a,x in p.items():
        for b,y in q.items():
            k=tuple(i+j for i,j in zip(a,b))
            out[k]=out.get(k,0)+x*y
    return clean(out)


def power(p,k,n=4):
    out={(0,)*n:Fraction(1)}
    for _ in range(k):out=mul(out,p)
    return out


def variables(n):
    return [{tuple(int(i==j) for i in range(n)):Fraction(1)} for j in range(n)]


def elementary(roots,k):
    n=len(next(iter(roots[0])))
    terms=[]
    for indices in itertools.combinations(range(len(roots)),k):
        p={(0,)*n:Fraction(1)}
        for i in indices:p=mul(p,roots[i])
        terms.append(p)
    return add(*terms)


def mod2(p):
    assert all(v.denominator==1 for v in p.values())
    return {k:int(v)%2 for k,v in p.items() if int(v)%2}


def sq2_roots(p):
    # Cartan formula on integral Chern roots: Sq^1 x=0, Sq^2 x=x^2.
    out={}
    for k,v in mod2(p).items():
        for i,m in enumerate(k):
            if m%2:
                new=list(k);new[i]+=1;new=tuple(new)
                out[new]=(out.get(new,0)+v)%2
    return clean(out)


def rows(p):
    return [dict(exponents=list(k),numerator=v.numerator,denominator=v.denominator)
            for k,v in sorted(p.items())]


def universal_invariant_check():
    x=variables(4)
    roots=x+[scale(add(*x),-1)]
    assert not add(*roots)
    weights=[add(*(roots[j] for j in s)) for s in original.STATES]
    ch1=add(*weights)
    tr3=add(*(power(w,3) for w in weights))
    ch2=scale(add(*(power(w,2) for w in weights)),Fraction(1,2))
    c2=elementary(roots,2)
    assert not ch1 and not tr3 and ch2==scale(c2,-4)
    wedge2_cubic=add(*(power(weights[i],3) for i,s in enumerate(original.STATES) if len(s)==2))
    wedge4_cubic=add(*(power(weights[i],3) for i,s in enumerate(original.STATES) if len(s)==4))
    fund3=add(*(power(w,3) for w in roots))
    assert wedge2_cubic==fund3 and wedge4_cubic==scale(fund3,-1) and fund3

    # Verify the characteristic-class identity as a polynomial, then restrict c1=0.
    allroots=variables(5)
    c1all=elementary(allroots,1);c2all=elementary(allroots,2);c3all=elementary(allroots,3)
    assert mod2(sq2_roots(c2all))==mod2(add(mul(c1all,c2all),c3all))
    c3=elementary(roots,3)
    assert mod2(sq2_roots(c2))==mod2(c3) and mod2(c3)
    # Restore the already-fixed original charge dictionary; no new embedding.
    J=original.dictionary()
    old_q=np.array([1]*6+[-4]*3+[2]*3+[-3]*2+[6,0])
    carrier_q=np.array([sum([-2,-2,-2,3,3][j] for j in s) for s in original.STATES])
    charge_error=float(np.max(abs(J.conj().T@np.diag(carrier_q)@J-np.diag(old_q))))
    assert charge_error==0.
    return dict(internal_dimension=len(weights),ch1_zero=not ch1,tr_F_cubed_zero=not tr3,
        ch2_equals_minus_four_c2=True,nonzero_wedge2_cubic=rows(wedge2_cubic),
        opposite_wedge4_cubic=rows(wedge4_cubic),ch2_polynomial=rows(ch2),
        Sq2_c2_equals_c1_c2_plus_c3_mod2=True,SU5_c3_mod2_nonzero=True,
        original_charge_dictionary_error=charge_error,
        topology_inputs='Spin AHSS low coefficient groups, H*(BSU5)=Z[c2,c3,c4,c5], and d2=Sq2_* after mod2 reduction are inherited theorems. Polynomial identities alone do not compute a bordism group.',
        derived_in_note='Only total-degree-five candidate E2_(4,1)=Z2 is killed by surjective d2 from (6,0).')


def nonliftable_bundle_check():
    # M=S2 x S2, L has c1=a+b and integral (a+b)^2=2.
    # E3=L+1+1, E2=L^-1+1; V=E3+E2 has determinant 1.
    q=[1,0,0,-1,0]
    weights=[sum(q[j] for j in s) for s in original.STATES]
    counts=Counter(weights)
    assert dict(counts)=={0:8,1:4,-1:4}
    J=original.dictionary()
    original_weights=[]
    for col in range(16):
        ids=np.flatnonzero(J[:,col])
        assert len(ids)==1 and abs(J[ids[0],col])==1
        original_weights.append(weights[int(ids[0])])
    names=('Q','uc','dc','L','ec','nuc');sizes=(6,3,3,2,1,1)
    modules={};start=0
    for name,size in zip(names,sizes):
        part=original_weights[start:start+size];start+=size
        # Ahat=1 in degree four here; index(L^q)=q^2.
        modules[name]=dict(line_powers=part,twisted_chiral_index=sum(w*w for w in part))
    c1_E2=[-1,-1];obstruction=[x%6 for x in c1_E2]
    assert obstruction==[5,5]
    c2_V=2*sum(q[i]*q[j] for i,j in itertools.combinations(range(5),2))
    index=sum(w*w for w in weights)
    assert c2_V==-2 and index==8 and sum(m['twisted_chiral_index'] for m in modules.values())==8
    assert index==-4*c2_V
    phi=np.zeros(5);h,d=original.old.mass_matrices(phi)
    assert np.max(abs(h))==0 and np.max(abs(d))==0
    F=float(original.old.original.F(phi))
    potential=float(original.old.original.node_potential(phi))
    assert F>0 and np.isfinite(potential)
    return dict(auxiliary_base='Closed Euclidean spin S2 x S2; not a derived Lorentz spacetime or a computed original graph continuum limit',
        cohomology='a^2=b^2=0, integral a*b=1, x=a+b',
        E3_line_powers=q[:3],E2_line_powers=q[3:],det_V_trivial=sum(q)==0,
        c1_E2=c1_E2,necessary_sixth_root_obstruction_mod6=obstruction,
        no_lift_to_SU3_times_SU2_times_U1=True,integral_c2_V=c2_V,
        integral_p1_TM=0,representation_line_multiplicities={str(k):v for k,v in sorted(counts.items())},
        module_indices=modules,full_twisted_chiral_index=index,
        lower_bound_on_massless_chiral_zero_modes=abs(index),
        original_scalar_zero=dict(F=F,potential=potential,all_original_mass_blocks_zero=True),
        scope='Index theorem forces massless continuum zero modes on this bundle. It does not determine the zero-mode count exactly, prove a zero of a massive/Yukawa-deformed Pfaffian, or construct a full continuum measure. Lack of anomaly does not imply a nowhere-zero determinant.')


def run():
    deps=('research_note_531.md','research_note_598.md','research_note_614.md',
          'research_note_615.md','research_note_616.md','research_note_617.md',
          'research_note_627.md','joint_spinor_subgroup_mass.py','joint_fermion_gauss_completion.py')
    return dict(round=628,tests_run=2,failures=0,errors=0,
        original_representation_invariants=universal_invariant_check(),
        original_quotient_bundle_and_index=nonliftable_bundle_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Conditional ordinary-spin four-dimensional continuum fermion anomaly interface for the original Z6 quotient and complete one-generation representation. Mature bordism/eta tools imply trivial closed-five-manifold anomaly phases, including non-liftable quotient bundles. No new spacetime dimension, gauge group, field content, local chiral regulator, regional CP process, nonbounding-sector phase choice, or GR dynamics is derived.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(round=628,tests_run=2,
        universal_invariants={k:v for k,v in result['original_representation_invariants'].items() if 'polynomial' not in k and 'cubic' not in k},
        bundle=result['original_quotient_bundle_and_index']),ensure_ascii=False,indent=2))

