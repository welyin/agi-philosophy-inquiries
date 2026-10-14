"""1208 exact permission/lattice certificates. Read-only unless --write."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent


def dec(x, upper=False, places=35):
    s = 10**places
    a = -((-x.numerator*s)//x.denominator) if upper else x.numerator*s//x.denominator
    return ('-' if a < 0 else '') + str(abs(a)//s)+'.'+str(abs(a)%s).zfill(places)


def native():
    old = json.loads((HERE.parent/'1202/results.json').read_text(encoding='utf-8'))
    p = [F(a) for a in old['one_marginal_rational']]
    boxes = [tuple(F(a, old['preparation']['root_denominator']) for a in b)
             for b in old['preparation']['root_boxes']]
    x, y = (0,30,0,20), (28,0,22,0)
    z = tuple(a-b for a,b in zip(x,y))
    assert sum(x) == sum(y) == 50 and z[2]%2 == 0
    a,b = (math.prod(p[i]**c[i] for i in range(4)) for c in (x,y))
    q, kappa = a+b, abs(a-b)/(a+b)
    m = (q.denominator+q.numerator-1)//q.numerator
    assert (m-1)*q < 1 <= m*q
    lo = sum(c*(box[0] if c >= 0 else box[1]) for c,box in zip(z,boxes))
    hi = sum(c*(box[1] if c >= 0 else box[0]) for c,box in zip(z,boxes))
    assert 0 < lo < hi < F('.0008474425022')
    # e > 163/60, hence 1-exp(-1) > 103/163. Strict, no huge exponent.
    lower = kappa*F(103,163)
    assert lower > F('.47')
    return {'occupation_difference':list(z), 'block_length':50,
        'original_single_z_parity_qualified':False,
        'doubled_z_parity_qualified':True, 'blocks':str(m),
        'total_materials':str(50*m), 'q_exact':str(q),
        'conditional_asymmetry_exact':str(kappa),
        'energy_interval':[dec(lo), dec(hi,True)],
        'trace_distance_strict_lower_exact':str(lower),
        'trace_distance_strict_lower_decimal':dec(lower),
        'physical_preparation_control_readout_certified':False}


def finite_charge_contrast():
    p = (F(1,2),F(1,4),F(1,4))
    x,y = (0,0,2,2),(1,1,1,1)
    a,b = math.prod(p[i] for i in x),math.prod(p[i] for i in y)
    assert a == F(1,64) and b == F(1,256)
    q,kappa=a+b,(a-b)/(a+b)
    rows=[]
    for blocks in (1,2):
        def perm(s):
            for start in range(0,len(s),4):
                part=s[start:start+4]
                if part in (x,y):
                    return s[:start]+(y if part==x else x)+s[start+4:]
            return s
        strings=list(product(range(3),repeat=4*blocks))
        dist=F(0); binary=F(0)
        for s in strings:
            t=perm(s)
            assert perm(t)==s and sum(s)==sum(t) and s.count(2)%2==t.count(2)%2
            ps,pt=(math.prod(p[i] for i in w) for w in (s,t))
            dist+=abs(ps-pt)/2
            first=next((s[j:j+4] for j in range(0,len(s),4) if s[j:j+4] in (x,y)),None)
            if first==x: binary+=ps-pt
        formula=(1-(1-q)**blocks)*kappa
        assert dist==binary==formula
        rows.append({'blocks':blocks,'basis_size':len(strings),'distance':str(dist),
            'energy_exactly_preserved':True,'parity_preserved':True,'first_x_attains_distance':True})
    d52=(1-(1-q)**52)*kappa
    assert d52>F('.38')
    # Exact generalized Gibbs weights with beta=gamma=log 2.
    weights=[F(2)**(-i+(1 if i==2 else 0)) for i in range(3)]
    assert tuple(w/sum(weights) for w in weights)==p
    # At fixed N2, probability ratios equal 2^(-energy difference), for all these strings.
    pairs=0
    for s in product(range(3),repeat=4):
        for t in product(range(3),repeat=4):
            if s.count(2)==t.count(2):
                ratio=math.prod(p[i] for i in t)/math.prod(p[i] for i in s)
                assert ratio==F(2)**(sum(s)-sum(t)); pairs+=1
    return {'probabilities':[str(v) for v in p], 'pair_probability':str(q),
        'conditional_asymmetry':str(kappa),'full_basis_checks':rows,
        '52_block_distance_exact':str(d52),'52_block_distance_lower':dec(d52),
        'exact_N2_generalized_Gibbs_beta_gamma':'log 2, log 2',
        'exact_N2_pair_checks':pairs, 'same_p_not_ordinary_Gibbs':p[1]**2!=p[0]*p[2]}


def permission_lattice():
    # z = z1*(e1-e0)+z3*(e3-e0)+(z2/2)*2(e2-e0).
    basis=(( -1,1,0,0),(-1,0,0,1),(-2,0,2,0))
    checked=0
    for b,c,d in product(range(-4,5),repeat=3):
        z=(-b-c-d,b,c,d)
        if c%2==0:
            coeff=(b,d,c//2)
            reconstructed=tuple(sum(coeff[j]*basis[j][i] for j in range(3)) for i in range(4))
            assert reconstructed==z; checked+=1
    # Original basis Q flips 1 and 2. h never connects 2 to another level.
    q={(1,2):1,(2,1):1}
    parity=(1,1,-1,1)
    assert all(parity[i]*parity[j]*v==-v for (i,j),v in q.items())
    edges=[((i,k),(j,l)) for i,j in q for k,l in q]
    assert all(parity[s[0]]*parity[s[1]]==parity[t[0]]*parity[t[1]] for s,t in edges)
    max_comm=max(abs(sum(v==2 for v in s)-sum(v==2 for v in t)) for s,t in edges)
    assert max_comm==2
    return {'parity_lattice_basis':[list(v) for v in basis],
        'index_in_sum_zero_lattice':2, 'bounded_integer_decompositions_checked':checked,
        'native_QQ_preserves_parity':True,'native_QQ_N2_commutator_nonzero_entry_abs':max_comm,
        'matrix_element_lattice_is_not_a_controllability_certificate':True,
        'preparation_count_preservation_requires_all_stocks_included':True}


def returned_finite_charge_reference():
    # Integer representation of W, and J=sqrt(2)*[psi -> psi tensor |+>].
    p=(1,-1,1); z=(1,-1)
    images={}
    for s,a in product(range(3),range(2)):
        t,b=(1-s,1-a) if s<2 else (s,a)
        assert p[s]*z[a]==p[t]*z[b]
        images[(s,a)]=(t,b)
    assert len(set(images.values()))==6
    for s in range(3):
        actual=sorted(images[(s,a)] for a in range(2))
        target=1-s if s<2 else s
        assert actual==[(target,0),(target,1)]
    assert {p[s]*p[1-s if s<2 else s] for s in range(3)}=={-1,1}
    # Two-level X tensor X also returns maximally mixed auxiliary, while X anticommutes Z.
    for s in range(2):
        assert sorted((1-s,1-a) for a in range(2))==[(1-s,0),(1-s,1)]
    return {'dimension':6,'W_is_parity_preserving_permutation':True,
        'WJ_equals_JT_exact_integer_columns':True,
        'auxiliary_plus_exact_return_for_all_unknown_system_R':True,
        'target_T_not_in_target_parity_preserving_class':True,
        'two_level_maximally_mixed_auxiliary_returns_under_X_tensor_X':True,
        'auxiliary_and_W_not_granted_to_native_QQ':True}


def main():
    names=['1201/proof.md','1202/results.json','1205/proof.md','1206/results.json',
           '1207/proof.md','1207/acceptance.json']
    result={'schema':'restricted_reference_1208_v1','status':'AUTHOR_CHECKS_PASS_NOT_ACCEPTANCE',
        'source_sha256':{name:hashlib.sha256((HERE.parent/name).read_bytes()).hexdigest() for name in names},
        'permission_lattice':permission_lattice(), 'finite_charge_contrast':finite_charge_contrast(),
        'finite_charge_return_counterexample':returned_finite_charge_reference(),
        'native_parity_witness':native(), 'all_checks_passed':True,
        'new_cognitive_axioms_adopted':0,
        'general_classification_requires_analytic_proof':True}
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true'); args=parser.parse_args()
    path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    assert json.loads(path.read_text(encoding='utf-8'))==result
    print(json.dumps({'status':'AUTHOR_CHECKS_PASS','all_checks_passed':True,'results_match':True,
        'native_parity_distance_lower':result['native_parity_witness']['trace_distance_strict_lower_decimal'],
        'native_parity_blocks':result['native_parity_witness']['blocks']}))


if __name__=='__main__':main()
