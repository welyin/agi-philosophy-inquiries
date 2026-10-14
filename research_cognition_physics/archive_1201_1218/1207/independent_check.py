"""Independent 1207 audit: exact first-hit permutations and rational bounds.
Uses only the Python standard library; imports no author or historical script.
Default is read-only and prints the complete matching JSON. --write creates once.
All-n, arbitrary-unitary and arithmetic-sequence claims are proved in the review;
finite enumerations below check the mechanisms, not these universal quantifiers.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from pathlib import Path
import argparse, hashlib, itertools, json, math
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
DEN=10**70

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(x,rounding=None):
    x=F(x)
    with localcontext() as ctx:
        ctx.prec=36
        if rounding: ctx.rounding=rounding
        return str(Decimal(x.numerator)/Decimal(x.denominator))
def product(seq):
    out=F(1)
    for x in seq: out*=x
    return out

def first_hit(word,x,y):
    out=list(word);k=len(x)
    assert len(y)==k and len(word)%k==0 and x!=y
    for start in range(0,len(word),k):
        block=tuple(word[start:start+k])
        if block==x:
            out[start:start+k]=y;break
        if block==y:
            out[start:start+k]=x;break
    return tuple(out)

def permutation_certificate(prob,energy,x,y,m):
    k=len(x);n=k*m
    words=list(itertools.product(range(len(prob)),repeat=n))
    prior={w:product(prob[i] for i in w) for w in words}
    mapping={w:first_hit(w,x,y) for w in words}
    assert len(set(mapping.values()))==len(words)
    assert all(mapping[mapping[w]]==w for w in words)
    after={mapping[w]:v for w,v in prior.items()}
    distance=sum(abs(after[w]-prior[w]) for w in words)/2
    norm=max(abs(sum(energy[i] for i in mapping[w])-sum(energy[i] for i in w)) for w in words)
    a=product(prob[i] for i in x);b=product(prob[i] for i in y);q=a+b
    formula=(1-(1-q)**m)*abs(a-b)/q
    expected_norm=abs(sum(energy[i] for i in x)-sum(energy[i] for i in y))
    assert distance==formula and norm==expected_norm
    return {'d':len(prob),'block_length':k,'blocks':m,'full_basis_dimension':len(words),
            'pair':[list(x),list(y)],'pair_probability':str(q),
            'exact_trace_distance':str(distance),'exact_operator_energy_error':str(norm),
            'formula_matches_full_permutation':True,'is_involution':True}

def f_int(x): return 32000*x**3-31200*x*x*DEN-2095*x*DEN**2-17*DEN**3

def isolate(a,b):
    lo=int(F(a)*DEN);hi=int(F(b)*DEN);flo=f_int(lo)
    assert flo*f_int(hi)<0
    while hi-lo>1:
        mid=(lo+hi)//2;fm=f_int(mid);assert fm
        if fm*flo>0:lo=mid;flo=fm
        else:hi=mid
    assert f_int(lo)*f_int(hi)<0
    return lo,hi

def signed_linear_interval(z,boxes):
    lo=sum(v*(b[0] if v>=0 else b[1]) for v,b in zip(z,boxes))
    hi=sum(v*(b[1] if v>=0 else b[0]) for v,b in zip(z,boxes))
    return F(lo,DEN),F(hi,DEN)

def coarse_hit_bound(q):
    assert 0<q<=1
    # m*q >= 5 and (1-q)^m <= exp(-m*q) <= exp(-5) < 1/100.
    # The last strict inequality follows from the positive exponential series.
    m=(5*q.denominator+q.numerator-1)//q.numerator
    assert m*q>=5
    series=sum(F(5**j,math.factorial(j)) for j in range(7))
    assert series>100
    return m,series

def main_result():
    p=[F(1,2),F(1,3),F(1,6)];E=[F(0),F(1),F(2)]
    tests=[permutation_certificate(p,E,(0,),(1,),m) for m in range(1,6)]
    tests += [permutation_certificate(p,E,(0,2),(1,1),m) for m in range(1,4)]
    # Equal energies with unequal populations: a genuine exact-kernel case.
    tests.append(permutation_certificate(p,[F(0),F(0),F(1)],(0,),(1,),4))
    # Swapping EVERY block accumulates operator energy error; first-hit does not.
    words=list(itertools.product(range(2),repeat=4))
    norm_all=max(abs(sum(1-i for i in w)-sum(w)) for w in words)
    norm_first=max(abs(sum(first_hit(w,(0,),(1,)))-sum(w)) for w in words)
    assert (norm_all,norm_first)==(4,1)
    # Without iid block inputs, the probability factorization need not hold.
    corr={(0,0):F(3,4),(1,1):F(1,4)}
    corr_after={first_hit(w,(0,),(1,)):v for w,v in corr.items()}
    support=set(corr)|set(corr_after)
    corrD=sum(abs(corr.get(w,0)-corr_after.get(w,0)) for w in support)/2
    iid_prediction=F(1,2)
    assert corrD==1 and corrD!=iid_prediction
    # Single-state mean energy does not certify operator energy error.
    pm=[F(1,2),F(1,6),F(1,3)];cycle=[1,2,0]
    mean_delta=sum(pm[i]*(E[cycle[i]]-E[i]) for i in range(3))
    after=[pm[cycle.index(i)] for i in range(3)]
    meanD=sum(abs(a-b) for a,b in zip(pm,after))/2
    mean_norm=max(abs(E[cycle[i]]-E[i]) for i in range(3))
    assert mean_delta==0 and meanD==F(1,3) and mean_norm==2
    # Positive/negative beta=+/-ln2 on E=(0,1,2), plus beta=0.
    # All tests satisfy the stronger rational bound 2*D^2/delta <= 2/3 < ln2.
    # ln2=2*(1/3+1/(3^3*3)+...) > 2/3.
    gibbs=[]
    for weights,beta in [([4,2,1],'+ln2'),([1,2,4],'-ln2'),([1,1,1],'0')]:
        pg=[F(x,sum(weights)) for x in weights]
        checked=0;max_ratio=F(0)
        for n in [1,2]:
            states=list(itertools.product(range(3),repeat=n))
            weightsn=[product(pg[i] for i in w) for w in states]
            energies=[sum(E[i] for i in w) for w in states]
            for i,j in itertools.combinations(range(len(states)),2):
                delta=abs(energies[i]-energies[j]);D=abs(weightsn[i]-weightsn[j])
                if delta:
                    ratio=2*D*D/delta
                    assert ratio<=F(2,3)
                    max_ratio=max(max_ratio,ratio)
                else: assert D==0
                if beta=='0': assert D==0
                checked+=1
        gibbs.append({'beta':beta,'population':[str(v) for v in pg],
                       'exact_transposition_checks':checked,
                       'max_2D_squared_over_delta':str(max_ratio),
                       'positive_and_negative_beta_retained':True})
    # The same Gibbs marginal can change its purifying reference correlations.
    overlap=F(2,3)-F(1,3)
    assert 1-overlap*overlap==F(8,9)
    # Native 1201 h and r*: same 1206 near-resonance, now first-hit amplified.
    boxes=[isolate('-.055','-.054'),isolate('-.010','-.009'),(DEN,DEN),isolate('1.038','1.039')]
    z=[-14,15,-11,10];eta=signed_linear_interval(z,boxes)
    assert sum(z)==0 and F('.0004237212509')<eta[0]<eta[1]<F('.0004237212510')
    rp=[F(3,8),F(1,8),F(1,8),F(3,8)]
    plus=[max(0,v) for v in z];minus=[max(0,-v) for v in z]
    k=sum(plus);assert k==sum(minus)==25
    a=product(r**c for r,c in zip(rp,plus));b=product(r**c for r,c in zip(rp,minus))
    q=a+b;asym=abs(a-b)/q
    assert b/a==81 and asym==F(40,41)
    m,series=coarse_hit_bound(q)
    lower=F(99,100)*asym;assert lower==F(198,205)
    # Main author claim uses the independently prepared original cold marginal.
    # Recompute it from frozen exact populations, not the new author results.
    colddata=json.loads((STAGE/'1202'/'results.json').read_text(encoding='utf-8'))
    pc=[F(t) for t in colddata['one_marginal_rational']]
    assert sum(pc)==1 and all(pc[i]>pc[i+1]>0 for i in range(3))
    activation_gap=signed_linear_interval([-1,1,1,-1],boxes)
    activation_population=pc[1]*pc[2]-pc[0]*pc[3]
    assert activation_gap[0]>F('.0060792930693608')
    assert activation_population>F('.0000188096985031')
    activation_work=(activation_gap[0]*activation_population,activation_gap[1]*activation_population)
    ca=product(r**c for r,c in zip(pc,plus))
    cb=product(r**c for r,c in zip(pc,minus))
    cq=ca+cb;contrast=abs(ca-cb)/cq
    cm=(cq.denominator+cq.numerator-1)//cq.numerator
    assert (cm-1)*cq<1<=cm*cq
    # A different conservative positive-term bound: e > sum_{j=0}^8 1/j!.
    eseries=sum(F(1,math.factorial(j)) for j in range(9))
    chit_lower=1-1/eseries
    cD_lower=chit_lower*contrast
    assert cm==4375693525386771 and k*cm==109392338134669275
    assert F('.4479105')<contrast<F('.4479107') and cD_lower>F('.283')
    # Exact Pell witnesses for a distinct illustrative spectrum (0,1,sqrt2).

    s2=math.isqrt(2*DEN*DEN)
    s2box=(F(s2,DEN),F(s2+1,DEN))
    assert s2*s2<2*DEN*DEN<(s2+1)**2
    pell=[];u,v=1,1
    for index in range(1,7):
        u,v=u+2*v,u+v
        assert abs(u*u-2*v*v)==1
        zz=[v-u,u,-v];assert sum(zz)==0
        raw=(F(u)-v*s2box[1],F(u)-v*s2box[0])
        absbox=(min(abs(t) for t in raw),max(abs(t) for t in raw))
        assert raw[0]*raw[1]>0 and absbox[1]<F(1,v)
        aa=F(1,3)**u;bb=F(1,2)**(u-v)*F(1,6)**v
        ratio=max(aa,bb)/min(aa,bb)
        qq=aa+bb;mm,_=coarse_hit_bound(qq)
        dd=F(99,100)*abs(aa-bb)/qq
        pell.append({'p':u,'q':v,'occupation_difference':zz,
            'block_length':u,'energy_error_lower':dec(absbox[0],ROUND_FLOOR),
            'energy_error_upper':dec(absbox[1],ROUND_CEILING),
            'probability_ratio_exact':str(ratio),'blocks_for_hit_gt_99_percent':str(mm),
            'trace_distance_strict_lower':dec(dd,ROUND_FLOOR),
            'different_spectrum_not_native_permission':True})
    assert all(F(a['energy_error_upper'])>F(b['energy_error_upper']) for a,b in zip(pell,pell[1:]))
    assert all(F(a['probability_ratio_exact'])<F(b['probability_ratio_exact']) for a,b in zip(pell,pell[1:]))
    files=['1201/proof.md','1201/results.json','1202/results.json','1206/proof.md','1206/results.json','1206/acceptance.json']
    return {'status':'INDEPENDENT_SCOPED_CHECKS_PASS_NOT_ACCEPTANCE','all_checks_passed':True,
        'method':'standard library exact enumeration, integer/rational root intervals; no author import',
        'source_sha256':{name:sha(STAGE/name) for name in files},
        'independent_script_sha256':sha(Path(__file__)),
        'first_hit_full_permutation_checks':tests,
        'operator_norm_accumulation_control':{'blocks':4,'swap_every_block_norm':norm_all,'first_hit_norm':norm_first},
        'iid_deletion':{'actual_correlated_trace_distance':str(corrD),'iid_formula_if_misapplied':str(iid_prediction)},
        'mean_energy_not_operator_budget':{'mean_change':str(mean_delta),'trace_distance':str(meanD),'operator_error':str(mean_norm)},
        'gibbs_sign_controls':gibbs,
        'arbitrary_reference_not_preserved':{'gibbs_population':['2/3','1/3'],'energy_error':'0','marginal_distance':'0','purification_distance_squared':'8/9'},
        'native_amplified_witness':{'population':[str(v) for v in rp],
            'occupation_difference':z,'positive_counts':plus,'negative_counts':minus,
            'block_length':k,'root_denominator':str(DEN),
            'root_integer_boxes':[[str(x) for x in box] for box in boxes],
            'energy_error_lower_exact':str(eta[0]),'energy_error_upper_exact':str(eta[1]),
            'declared_total_operator_energy_budget':'.000424',
            'probability_x':str(a),'probability_y':str(b),'hit_probability_one_block':str(q),
            'conditional_asymmetry':str(asym),'blocks':str(m),'total_material_copies':str(k*m),
            'm_times_q':str(m*q),'positive_exp5_series_partial_sum':str(series),
            'hit_probability_strict_lower':'99/100','trace_distance_strict_lower':str(lower),
            'trace_distance_lower_decimal':dec(lower,ROUND_FLOOR),
            'resources_enormous_finite_not_experimental_budget':True},
        'original_cold_amplified_witness':{
            'preparation':'independent cold pools for every material copy; not many correlated draws from one pool',
            'occupation_difference':z,'block_length':k,
            'population_exact':[str(t) for t in pc],
            'non_gibbs_two_copy_inversion':{
                'higher_energy_string':[1,2],'lower_energy_string':[0,3],
                'energy_gap_lower':dec(activation_gap[0],ROUND_FLOOR),
                'energy_gap_upper':dec(activation_gap[1],ROUND_CEILING),
                'higher_energy_probability_excess_exact':str(activation_population),
                'ideal_cyclic_work_lower':dec(activation_work[0],ROUND_FLOOR),
                'ideal_cyclic_work_upper':dec(activation_work[1],ROUND_CEILING),
                'not_a_battery_or_control_implementation':True},
            'pair_probability_exact':str(cq),'pair_probability_decimal':dec(cq),
            'conditional_asymmetry_exact':str(contrast),'conditional_asymmetry_decimal':dec(contrast),
            'blocks':str(cm),'total_material_copies':str(k*cm),
            'm_times_q_exact':str(cm*cq),
            'positive_e_series_partial_sum':str(eseries),
            'hit_probability_strict_lower':str(chit_lower),
            'trace_distance_strict_lower_exact':str(cD_lower),
            'trace_distance_strict_lower_decimal':dec(cD_lower,ROUND_FLOOR),
            'declared_trace_distance_lower':'.283',
            'same_energy_bound_as_native_r_star':True,
            'not_1205_resources_or_declared_thermometry_precision':True},
        'pell_sequence_illustration':pell,
        'analytic_scope':{'gibbs_modulus':'min(1,sqrt(abs(beta)*delta/2)) for all finite n and all U',
            'non_gibbs_uniform_response':'R_r(delta)=1 for each delta>0; supremum, not finite full-rank attainment',
            'native_r_star_exact_response':'R_r(0)=0 by frozen 1206',
            'finite_copy_cutoff':'For native r*, every fixed n cutoff admits a modulus by compactness; no finite test certifies all n.',
            'permissions':'All finite global unitaries form a mathematical upper permission class; CC does not grant these native operations.',
            'not_claimed':['reference preparation','thermalization','beta nonnegative','actual first-hit controller','autonomous work cycle','geometry identity','arbitrary purifying R return']}}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    r=main_result();dest=HERE/'independent_results.json'
    if args.write:
        if dest.exists(): raise SystemExit('refuse overwrite of independent_results.json')
        dest.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads(dest.read_text(encoding='utf-8'))
        if saved!=r: raise SystemExit('independent result mismatch')
    print(json.dumps(r,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
