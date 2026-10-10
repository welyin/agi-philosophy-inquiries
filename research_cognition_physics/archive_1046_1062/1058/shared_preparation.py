"""Shared-preparation process checks; finite matrices plus exact Pauli support."""
from pathlib import Path
from fractions import Fraction
from itertools import product
import importlib.util
import argparse
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('round1057',HERE.parent/'1057/composition_stability.py')
base=importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def allowed_after_grouping(words):
    if all(all(k==0 for k in w) for w in words): return True
    for input_leg,output_leg in ((0,1),(2,3)):
        if all(w[output_leg]==0 for w in words) and any(w[input_leg]!=0 for w in words):
            return True
    return False


def sharp_support(n,t):
    atoms=[(0,0,0,0),(0,3,3,0),(3,0,0,3)]
    entries=[]
    for labels in product(range(3),repeat=n):
        # Averaging the same sign s over all copies removes odd powers of Q.
        if sum(k==2 for k in labels)%2: continue
        coefficient=t**sum(k!=0 for k in labels)
        words=tuple(atoms[k] for k in labels)
        entries.append((words,coefficient,allowed_after_grouping(words)))
    return entries


def calculate():
    rng=np.random.default_rng(1058)
    words=base.R_TYPES+base.P_TYPES+base.Q_TYPES
    basis=np.array([base.pauli_word(w) for w in words])
    mixture_rows=[]
    max_born_error=max_four_moment_error=max_restriction_error=0.0
    dense_pairings=0
    for case in range(8):
        mu=np.array([0.2,0.3,0.5])
        states=[]; approximants=[]; cs=[]; ds=[]
        for _ in range(3):
            h=np.einsum('a,aij->ij',rng.normal(size=87),basis[1:])
            h*=0.8/np.max(abs(np.linalg.eigvalsh(h)))
            w=(np.eye(16)+h)/16
            assert np.linalg.eigvalsh(w).min()>0
            coeff=np.einsum('ij,aji->a',w,basis).real
            c,d=coeff[16:52],coeff[52:]
            v=base.twirl(w,(3,)) if np.dot(d,d)<=np.dot(c,c) else base.twirl(w,(1,))
            states.append(w); approximants.append(v); cs.append(c); ds.append(d)
            assert base.trace_distance(w,v) <= .5*(np.dot(c,c)*np.dot(d,d))**.25+1e-12
        states=np.array(states); approximants=np.array(approximants)
        cs=np.array(cs); ds=np.array(ds)
        gg=np.einsum('ki,kj->kij',cs,ds)
        first=np.einsum('k,kij->ij',mu,gg)
        second=np.einsum('k,kij->ij',mu,gg**2)
        four_plus=1+2*first+second; four_minus=1-2*first+second
        eps=float(max(np.max(abs(four_plus-1)),np.max(abs(four_minus-1))))
        assert np.max(second)<=eps+1e-12 and np.min(second)>=0
        mean_distance=sum(float(p)*base.trace_distance(w,v) for p,w,v in zip(mu,states,approximants))
        moment_bound=.5*float(second.sum())**.25
        assert mean_distance<=moment_bound+1e-12
        assert moment_bound<=3*eps**.25+1e-12
        marginal=np.einsum('k,kij->ij',mu,states)
        approximate=np.einsum('k,kij->ij',mu,approximants)
        first_distance=base.trace_distance(marginal,approximate)
        assert first_distance<=mean_distance+1e-12
        two=sum(float(p)*np.kron(w,w) for p,w in zip(mu,states))
        two_approx=sum(float(p)*np.kron(v,v) for p,v in zip(mu,approximants))
        second_distance=base.trace_distance(two,two_approx)
        assert second_distance<=2*mean_distance+1e-12
        reduced=np.einsum('iaja->ij',two_approx.reshape(16,16,16,16))
        max_restriction_error=max(max_restriction_error,float(np.max(abs(reduced-approximate))))
        assert np.linalg.eigvalsh(two_approx).min()>=-1e-12
        forbidden=two_approx-base.twirl(two_approx,(1,5))-base.twirl(two_approx,(3,7))+base.twirl(two_approx,(1,3,5,7))
        assert np.linalg.norm(forbidden)<1e-12
        for _ in range(6):
            i,j=int(rng.integers(36)),int(rng.integers(36))
            p,q=base.P_TYPES[i],base.Q_TYPES[j]
            ua,va=np.kron(base.PAULI[p[0]],base.PAULI[q[0]]),np.kron(base.PAULI[p[1]],base.PAULI[q[1]])
            ub,vb=np.kron(base.PAULI[p[2]],base.PAULI[q[2]]),np.kron(base.PAULI[p[3]],base.PAULI[q[3]])
            four=[]
            for sign in (-1,1):
                ma,_=base.channel(ua,va,1)
                mb,_=base.channel(ub,vb,sign)
                z=[]
                for w,g in zip(states,gg[:,i,j]):
                    grouped=16*base.reorder(np.kron(w,w),[0,4,1,5,2,6,3,7])
                    value=base.pairing(grouped,np.kron(ma,mb))
                    max_born_error=max(max_born_error,abs(value-(1+sign*g)))
                    z.append(value)
                    dense_pairings+=1
                # The four-copy operator is a product of these two blocks in
                # each fixed branch; average only after squaring the contraction.
                value=float(np.dot(mu,np.square(z)))
                expected=1+2*sign*first[i,j]+second[i,j]
                max_four_moment_error=max(max_four_moment_error,abs(value-expected))
                four.append(value)
            assert abs(sum(four)/2-1-second[i,j])<2e-12
        mixture_rows.append({'case':case,'four_copy_menu_defect':eps,
                             'mean_component_twirl_distance':mean_distance,
                             'actual_moment_bound':moment_bound,'uniform_bound':3*eps**.25,
                             'first_marginal_distance':first_distance,
                             'two_marginal_distance':second_distance})
    exact=[]
    P=base.pauli_word((0,3,3,0)); Q=base.pauli_word((3,0,0,3))
    for t in (Fraction(1,8),Fraction(1,4),Fraction(1,2)):
        ss=[(np.eye(16)+float(t)*P+s*float(t)*Q)/16 for s in (-1,1)]
        assert min(np.linalg.eigvalsh(w).min() for w in ss)>-1e-12
        two=sum(np.kron(w,w) for w in ss)/2
        forbidden=two-base.twirl(two,(1,5))-base.twirl(two,(3,7))+base.twirl(two,(1,3,5,7))
        assert np.linalg.norm(forbidden)<1e-12
        support2=sharp_support(2,t); support4=sharp_support(4,t)
        assert all(allowed for _,_,allowed in support2)
        bad4=[(w,c) for w,c,allowed in support4 if not allowed]
        assert bad4 and all(c>0 for _,c in bad4)
        target=((0,3,3,0),(3,0,0,3),(0,3,3,0),(3,0,0,3))
        assert dict(bad4)[target]==t**4
        e_g=sum((s*t*t for s in (-1,1)),Fraction(0))/2
        e_g2=sum(((s*t*t)**2 for s in (-1,1)),Fraction(0))/2
        assert e_g==0 and e_g2==t**4
        exact.append({'t':str(t),'two_copy_all_grouped_types_valid':True,
                      'two_copy_nonzero_types':len(support2),'four_copy_nonzero_types':len(support4),
                      'four_copy_forbidden_types':len(bad4),
                      'two_copy_both_signed_total_weights':'1',
                      'four_copy_both_signed_total_weights':str(1+t**4),
                      'first_marginal_already_fixed_order':True})
    assert max(max_born_error,max_four_moment_error,max_restriction_error)<2e-11
    return {'round':1058,'all_checks_passed':True,'fixed_four_copy_settings':2592,
            'mathematical_representation_is_an_extra_input':True,
            'mixtures':8,'components_per_mixture':3,'dense_two_copy_born_pairings':dense_pairings,
            'four_copy_checks_use_tensor_factorization_not_dense_65536_matrices':True,
            'max_dense_born_error':max_born_error,'max_four_copy_moment_error':max_four_moment_error,
            'max_approximant_restriction_error':max_restriction_error,
            'mixture_diagnostics':mixture_rows,'exact_hidden_sign_examples':exact,
            'quarter_power_optimality_claimed':False,'empirical_data':False,
            'new_cognitive_axioms':0,'independent_review_included':False}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--save',action='store_true')
    args=parser.parse_args()
    r=calculate(); path=HERE/'results.json'
    if args.save:
        with path.open('x',encoding='utf-8') as out:
            json.dump(r,out,ensure_ascii=False,indent=2); out.write('\n')
    else:
        base.compare(json.loads(path.read_text(encoding='utf-8')),r)
    print(json.dumps({k:v for k,v in r.items() if k not in ('mixture_diagnostics','exact_hidden_sign_examples')},ensure_ascii=False))
