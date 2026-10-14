"""1208 independent exact audit; no author or historical code import.
Default recomputes read-only and prints the complete saved JSON. --write creates once.
The general subgroup classification is an analytic theorem, not a finite test claim.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
import argparse, hashlib, itertools, json, math
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
DEN=10**70

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(x,up=False):
    x=F(x)
    with localcontext() as ctx:
        ctx.prec=42;ctx.rounding=ROUND_CEILING if up else ROUND_FLOOR
        return str(Decimal(x.numerator)/Decimal(x.denominator))
def prod(xs):
    a=F(1)
    for x in xs:a*=x
    return a

def counts(w,d):return tuple(w.count(i) for i in range(d))
def diff(a,b):return tuple(x-y for x,y in zip(a,b))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def first_hit(w,x,y):
    k=len(x);out=list(w)
    for j in range(0,len(w),k):
        block=tuple(w[j:j+k])
        if block==x:out[j:j+k]=y;break
        if block==y:out[j:j+k]=x;break
    return tuple(out)

def matrix_product(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def commutator(a,b):
    ab=matrix_product(a,b);ba=matrix_product(b,a)
    return [[x-y for x,y in zip(ar,br)] for ar,br in zip(ab,ba)]
def diag(x):return [[int(i==j)*x[i] for j in range(len(x))] for i in range(len(x))]
def kron(a,b):
    return [[a[i][j]*b[k][l] for j in range(len(a[0])) for l in range(len(b[0]))]
            for i in range(len(a)) for k in range(len(b))]
def zero(a):return all(v==0 for row in a for v in row)

def first_hit_cases():
    p=[F(1,2),F(1,4),F(1,4)];E=(0,1,2)
    x=(0,0,2,2);y=(1,1,1,1);z=diff(counts(x,3),counts(y,3))
    assert z==(2,-4,2) and dot(z,E)==0 and z[2]%2==0 and z[2]!=0
    a=prod(p[i] for i in x);b=prod(p[i] for i in y)
    assert a==F(1,64) and b==F(1,256)
    q=a+b;contrast=abs(a-b)/q
    assert q==F(5,256) and contrast==F(3,5)
    cases=[]
    for m in (1,2):
        words=list(itertools.product(range(3),repeat=4*m))
        perm={w:first_hit(w,x,y) for w in words}
        probabilities={w:prod(p[i] for i in w) for w in words}
        assert len(set(perm.values()))==len(words)
        assert all(perm[perm[w]]==w for w in words)
        assert all((counts(perm[w],3)[2]-counts(w,3)[2])%2==0 for w in words)
        norm=max(abs(sum(E[i] for i in perm[w])-sum(E[i] for i in w)) for w in words)
        D=sum(abs(probabilities[w]-probabilities[perm[w]]) for w in words)/2
        assert norm==0 and D==(1-(1-q)**m)*contrast
        cases.append({'blocks':m,'dimension':len(words),'exact_operator_energy_error':str(norm),
            'exact_trace_distance':str(D),'parity_preserved':True,'involution':True})
    exact52=(1-(1-q)**52)*contrast
    assert exact52>F(38,100)
    return {'energy':[0,1,2],'population':['1/2','1/4','1/4'],
            'x':list(x),'y':list(y),'occupation_difference':list(z),
            'probability_x':str(a),'probability_y':str(b),'q':str(q),'contrast':str(contrast),
            'small_full_basis_checks':cases,'blocks_main':52,
            'trace_distance_main_exact':str(exact52),'trace_distance_main_lower':dec(exact52),
            'declared_strict_lower':'.38','forbidden_by_exact_N2':True}

def additive_charge_controls():
    E=(0,1,2)
    controls=[([F(1,2),F(1,4),F(1,4)],[-1,-2,-2],1),
              ([F(1,4),F(1,2),F(1,4)],[-2,-1,-2],-1),
              ([F(1,4),F(1,4),F(1,2)],[-2,-2,-1],0)]
    result=[]
    for p,log2p,beta_sign in controls:
        c=[a+beta_sign*b for a,b in zip(log2p,E)]
        assert c[0]==c[1] # c is constant + gamma*N2, hence in ann(L_N2).
        pair_checks=0;rotation_checks=0;largest=F(0)
        for n in (1,2,3):
            words=list(itertools.product(range(3),repeat=n))
            probs=[prod(p[i] for i in w) for w in words]
            energies=[sum(E[i] for i in w) for w in words]
            charges=[w.count(2) for w in words]
            cs=[sum(c[i] for i in w) for w in words]
            for i,j in itertools.combinations(range(len(words)),2):
                if charges[i]!=charges[j]:continue
                assert cs[i]==cs[j]
                for sin in (F(1),F(4,5)):
                    delta=abs(energies[i]-energies[j])*sin
                    D=abs(probs[i]-probs[j])*sin
                    if beta_sign==0:assert D==0
                    if delta:
                        ratio=2*D*D/delta
                        assert ratio<=F(2,3) # ln 2 > 2/3 by its positive atanh series.
                        largest=max(largest,ratio)
                    else:assert D==0
                    if sin==1:pair_checks+=1
                    else:rotation_checks+=1
        result.append({'population':[str(t) for t in p],'logp_over_ln2':log2p,
            'beta_over_ln2':beta_sign,'conserved_C_coefficients_over_ln2':c,
            'allowed_transpositions_checked':pair_checks,'nontrivial_coherent_rotations_checked':rotation_checks,
            'max_2D_squared_over_delta':str(largest),'beta_zero_exact_invariance':beta_sign==0})
    # With L={0}, equal-count words have identical product probabilities for any p.
    p=[F(1,6),F(1,3),F(1,2)];bycount={};n=5
    for w in itertools.product(range(3),repeat=n):
        key=counts(w,3);value=prod(p[i] for i in w)
        if key in bycount:assert bycount[key]==value
        else:bycount[key]=value
    return {'exact_N2_controls':result,'L_zero_count_blocks':len(bycount),
            'L_zero_arbitrary_population':['1/6','1/3','1/2'],
            'finite_index_example_basis':[[1,-1,0],[0,-2,2]],
            'finite_index_relative_to_Lambda':2,
            'analytic_bound':'min(1,sqrt(abs(beta)*delta/2)); arbitrary U by conserved C, not these finite checks'}

def native_QQ_algebra():
    S=diag([1,1,-1,1]);Q=[[0]*4 for _ in range(4)];Q[1][2]=Q[2][1]=1
    N=diag([0,0,1,0]);I=diag([1]*4)
    # h support is enough for exact [S,h]=[N,h]=0; includes irrational h_03.
    support=[(0,1),(1,0),(0,3),(3,0),(1,1),(2,2),(3,3)]
    assert all(S[i][i]==S[j][j] and N[i][i]==N[j][j] for i,j in support)
    SQS=matrix_product(matrix_product(S,Q),S)
    assert all(SQS[i][j]==-Q[i][j] for i in range(4) for j in range(4))
    V=kron(Q,Q);parity=kron(S,S)
    ns=diag([int(i==2)+int(j==2) for i in range(4) for j in range(4)])
    assert zero(commutator(parity,V))
    C=commutator(ns,V)
    assert C[10][5]==2 and C[5][10]==-2
    assert max(sum(abs(x) for x in row) for row in C)==2
    assert sum(C[i][5]**2 for i in range(16))==4 # norm >=2; row/column bounds <=2.
    # Energy-basis Q is the certified c-star; actual nonzero matrix elements are not controls.
    differences=set()
    for s in itertools.product(range(4),repeat=2):
        for t in itertools.product(range(4),repeat=2):
            if all((a==2 and b!=2) or (a!=2 and b==2) for a,b in zip(s,t)):
                differences.add(diff(counts(t,4),counts(s,4)))
    basis=[(-1,1,0,0),(-1,0,0,1),(-2,0,2,0)]
    assert all(b in differences for b in basis)
    for z in differences:
        assert sum(z)==0 and z[2]%2==0
        coeff=(z[1],z[3],z[2]//2)
        reconstructed=tuple(sum(coeff[j]*basis[j][i] for j in range(3)) for i in range(4))
        assert reconstructed==z
    # Index in Lambda, in coordinates (z1,z3,z2), is det diag(1,1,2)=2.
    return {'S_diagonal':[1,1,-1,1],'h_commutes_S_and_N2':True,'SQS_equals_minus_Q':True,
        'QQ_commutes_total_parity':True,'QQ_commutator_with_total_N2_norm':'2',
        'full_H_point1_QQ_commutator_with_total_N2_norm':'1/5',
        'QQ_star_distinct_occupation_differences':[list(v) for v in sorted(differences)],
        'integer_basis':[list(b) for b in basis],'index_in_Lambda':2,'real_rank':3,
        'all_generators_integer_decomposed':True,'basis_in_actual_matrix_element_support':True,
        'full_parity_block_unitary_reachability_claimed':False,
        'stock_copy_boolean_only':'Conserve all material-plus-inventory N_i; a direction is not an isolated programmable gate.'}

def finite_parity_catalyst_control():
    X=[[0,1],[1,0]];Z=[[1,0],[0,-1]];I=diag([1,1])
    W=kron(X,X);P=kron(Z,Z)
    assert zero(commutator(W,P)) and not zero(commutator(X,Z))
    # J maps psi -> psi tensor (|0>+|1>); normalization 1/sqrt2 cancels.
    J=[[1,0],[1,0],[0,1],[0,1]]
    assert matrix_product(W,J)==matrix_product(J,X)
    mixed=diag([F(1,2),F(1,2)])
    assert matrix_product(matrix_product(X,mixed),X)==mixed
    # Stronger: T mixes even and odd parity characters, so is not merely a
    # parity-covariant unitary channel with a global sign on its implementing U.
    Ps=diag([1,-1,1]);T=[[0,1,0],[1,0,0],[0,0,1]]
    flip=[[0,1,0],[1,0,0],[0,0,0]];spectator=diag([0,0,1])
    A=kron(flip,X);B=kron(spectator,I)
    W3=[[a+b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
    J3=[[int(i//2==j) for j in range(3)] for i in range(6)]
    assert matrix_product(W3,W3)==diag([1]*6)
    assert zero(commutator(W3,kron(Ps,Z)))
    assert matrix_product(W3,J3)==matrix_product(J3,T)
    ptp=matrix_product(matrix_product(Ps,T),Ps)
    assert ptp!=T and ptp!=[[-x for x in row] for row in T]
    return {'binary_product_W_X_tensor_X':{
        'total_parity_conserved':True,'target_U_does_not_commute_parity':True,
        'joint_isometry_WJ_equals_JX':True,'pure_plus_aux_returned_exactly':True,
        'mixed_I_over_2_aux_marginal_also_returned':True,
        'arbitrary_purification_of_mixed_aux_not_claimed_returned':True},
        'qutrit_controlled_flip_W':{'target_parity':[1,-1,1],
            'target_T':'swap(0,1) plus identity on level 2',
            'total_parity_conserved':True,'unitary_involution':True,
            'joint_isometry_WJ_equals_JT':True,'P_T_P_is_neither_T_nor_minus_T':True,
            'auxiliary':'pure |+>, returned exactly and uncorrelated for every target SR'},
        'meaning':'Global discrete parity plus auxiliary return does not derive target-only U_L.',
        'native_preparation_or_realization_claimed':False,
        'restriction_needed':'Target U_L is an independent effective-menu contract; absence of a supplied plus reference is not a general return-cycle no-go.'}

def polynomial_int(x):return 32000*x**3-31200*x*x*DEN-2095*x*DEN**2-17*DEN**3
def isolate(a,b):
    lo=int(F(a)*DEN);hi=int(F(b)*DEN);flo=polynomial_int(lo)
    assert flo*polynomial_int(hi)<0
    while hi-lo>1:
        mid=(lo+hi)//2;fm=polynomial_int(mid);assert fm
        if fm*flo>0:lo=mid;flo=fm
        else:hi=mid
    return lo,hi

def native_doubled_witness():
    p=[F(v) for v in json.loads((STAGE/'1202'/'results.json').read_text(encoding='utf-8'))['one_marginal_rational']]
    assert sum(p)==1 and min(p)>0
    boxes=[isolate('-.055','-.054'),isolate('-.010','-.009'),(DEN,DEN),isolate('1.038','1.039')]
    oldz=(-14,15,-11,10);z=tuple(2*t for t in oldz)
    assert oldz[2]%2==1 and z[2]%2==0 and sum(z)==0
    xp=[max(t,0) for t in z];yp=[max(-t,0) for t in z]
    k=sum(xp);assert k==sum(yp)==50
    lo=F(sum(c*(a if c>=0 else b) for c,(a,b) in zip(z,boxes)),DEN)
    hi=F(sum(c*(b if c>=0 else a) for c,(a,b) in zip(z,boxes)),DEN)
    assert F('.0008474425018')<lo<hi<F('.0008474425020')
    a=prod(v**j for v,j in zip(p,xp));b=prod(v**j for v,j in zip(p,yp));q=a+b
    contrast=abs(a-b)/q
    assert contrast>F('.746')
    m=(q.denominator+q.numerator-1)//q.numerator
    assert (m-1)*q<1<=m*q
    # Independent stronger finite bound e > its positive series through order 8.
    series=sum(F(1,math.factorial(j)) for j in range(9))
    lower=(1-1/series)*contrast
    assert lower>F('.47')
    assert m==31894572775086737403576853338681
    return {'source':'independent one-material cold marginals from frozen 1202',
        'old_occupation_difference':list(oldz),'old_rejected_by_parity':True,
        'new_occupation_difference':list(z),'positive_counts':xp,'negative_counts':yp,'block_length':k,
        'root_denominator':str(DEN),'root_integer_boxes':[[str(v) for v in box] for box in boxes],
        'energy_error_lower':dec(lo),'energy_error_upper':dec(hi,True),
        'declared_energy_bound':'.0008474425022',
        'q_lower':dec(q),'q_upper':dec(q,True),'contrast_lower':dec(contrast),'contrast_upper':dec(contrast,True),
        'blocks':str(m),'total_material_copies':str(k*m),'ceil_certificate_exact':True,
        'positive_exp1_series':str(series),'trace_distance_lower':dec(lower),
        'trace_distance_lower_using_author_coarser_hit_bound':dec(contrast*F(103,163)),
        'declared_trace_distance_lower':'.47','parity_preserved':True,
        'resources_or_control_not_certified':True}

def result():
    names=['1201/proof.md','1202/results.json','1203/proof.md','1205/proof.md',
           '1206/proof.md','1206/results.json','1207/proof.md','1207/results.json','1207/acceptance.json']
    return {'status':'INDEPENDENT_SCOPED_CHECKS_PASS_NOT_ACCEPTANCE','all_checks_passed':True,
        'method':'standard-library exact fractions, complete small permutations, algebraic support, independent root bisection; no author import',
        'source_sha256':{name:sha(STAGE/name) for name in names},'independent_script_sha256':sha(Path(__file__)),
        'qutrit_parity_witness':first_hit_cases(),'additive_charge_and_sign_controls':additive_charge_controls(),
        'native_QQ_algebra':native_QQ_algebra(),'finite_parity_catalyst_control':finite_parity_catalyst_control(),
        'native_doubled_cold_witness':native_doubled_witness(),
        'scope':{'classification':'uniform robustness iff log(p)+beta E annihilates real span(L); otherwise every positive delta has supremum 1',
            'sufficiency':'all actual tested comparisons are certified material unitaries preserving the relevant real additive charge; generic CPTP requires entropy/resource terms',
            'necessity':'the required first-hit family and binary readout must be available with full resources',
            'finite_index':'same real span does not imply same implementation cost or full controllability',
            'arbitrary_R':'system marginal only; no joint return guarantee',
            'not_claimed':['actual first-hit implementation','exact Gibbs preparation','beta positivity','thermalization','geometry','automatic adoption of a cognitive axiom']}}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    r=result();dest=HERE/'independent_results.json'
    if args.write:
        if dest.exists():raise SystemExit('refuse overwrite of existing independent_results.json')
        dest.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        if json.loads(dest.read_text(encoding='utf-8'))!=r:raise SystemExit('independent result mismatch')
    print(json.dumps(r,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
