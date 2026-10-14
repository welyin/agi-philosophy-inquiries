"""1206 independent exact algebra; does not import the author check.
Default: recompute, compare the saved JSON, print it. --write creates once.
The all-n claims use the proof stated in independent_review.md. Finite gap
calculations below are diagnostics, not thermalization or device evidence.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
DEN = 10**70
COEFF = [-17, -2095, -31200, 32000]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def decimal(q, places=32, rounding=None):
    q=F(q)
    with localcontext() as ctx:
        ctx.prec=places
        if rounding is not None: ctx.rounding=rounding
        return str(Decimal(q.numerator)/Decimal(q.denominator))

def add(a,b):
    c=[F(0)]*max(len(a),len(b))
    for i,x in enumerate(a): c[i]+=x
    for i,x in enumerate(b): c[i]+=x
    return c

def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c

def scale(a,t): return [t*x for x in a]

def poly_integer(x):
    return 32000*x**3-31200*x*x*DEN-2095*x*DEN**2-17*DEN**3

def isolate(a,b):
    lo=int(F(a)*DEN); hi=int(F(b)*DEN)
    flo=poly_integer(lo)
    assert flo*poly_integer(hi)<0
    while hi-lo>1:
        mid=(hi+lo)//2; fm=poly_integer(mid)
        assert fm != 0
        if fm*flo>0: lo=mid; flo=fm
        else: hi=mid
    assert poly_integer(lo)*poly_integer(hi)<0
    return (lo,hi)

def row_rank(rows):
    a=[list(map(F,row)) for row in rows]
    if not a: return 0
    rank=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(rank,len(a)) if a[i][col]),None)
        if pivot is None: continue
        a[rank],a[pivot]=a[pivot],a[rank]
        q=a[rank][col]; a[rank]=[v/q for v in a[rank]]
        for i in range(len(a)):
            if i!=rank and a[i][col]:
                q=a[i][col]
                a[i]=[v-q*w for v,w in zip(a[i],a[rank])]
        rank+=1
    return rank

def types(n):
    for a in range(n+1):
        for b in range(n-a+1):
            for c in range(n-a-b+1):
                yield (a,b,c,n-a-b-c)

def gaps(n,boxes):
    values=[]
    for counts in types(n):
        lo=sum(k*b[0] for k,b in zip(counts,boxes))
        hi=sum(k*b[1] for k,b in zip(counts,boxes))
        values.append((lo+hi,lo,hi,counts))
    values.sort()
    assert len(values)==math.comb(n+3,3)
    candidates=[]
    for v,w in zip(values,values[1:]):
        lower=w[1]-v[2]; upper=w[2]-v[1]
        assert lower>0
        candidates.append((lower,upper,v[3],w[3]))
    lower=min(v[0] for v in candidates)
    upper=min(v[1] for v in candidates)
    best=min(candidates,key=lambda v:v[0]+v[1])
    width_upper=F(boxes[3][1]-boxes[0][0],DEN)
    pigeon=F(n)*width_upper/(len(values)-1)
    assert F(upper,DEN)<=pigeon
    return {
        'n':n,'distinct_occupation_energies':len(values),
        'all_adjacent_interval_gaps_strictly_positive':True,
        'minimum_gap_lower':decimal(F(lower,DEN),rounding=ROUND_FLOOR),
        'minimum_gap_upper':decimal(F(upper,DEN),rounding=ROUND_CEILING),
        'minimum_gap_integer_bounds':[str(lower),str(upper)],
        'one_minimum_midpoint_pair':[list(best[2]),list(best[3])],
        'pigeonhole_upper':decimal(pigeon,rounding=ROUND_CEILING),
        'scope':'finite spectrum diagnostic; no collision permission or observable response lower bound'
    }

def main_result():
    # Independent determinant of xI-h on native indices 0,1,3.
    cubic=add(add(mul(mul([0,1],[-1,1]),[F(1,40),1]),
                  scale([F(1,40),1],F(-1,25))),
              scale([-1,1],F(-3,6400)))
    assert cubic==[F(x,32000) for x in COEFF]
    quartic=mul(cubic,[-1,1])
    residues=[sum(c*pow(x,k,11) for k,c in enumerate(COEFF))%11
              for x in range(11)]
    assert all(residues) and 32000%11
    a,b,c,d=32000,-31200,-2095,-17
    disc=b*b*c*c-4*a*c**3-4*b**3*d-27*a*a*d*d+18*a*b*c*d
    sq=math.isqrt(disc)
    assert sq*sq<disc<(sq+1)*(sq+1)
    assert disc==2736128576000000
    boxes=[isolate('-0.055','-0.054'),isolate('-0.010','-0.009'),
           (DEN,DEN),isolate('1.038','1.039')]
    assert all(boxes[i][1]<boxes[i+1][0] for i in range(3))
    f1=sum(COEFF); assert f1==-1312
    relation=[40,40,-39,40]
    assert 40*F(39,40)-39==0 and sum(relation)==81
    determinant=F(3)-F(39,40); assert determinant==F(81,40)
    # QQ support in the energy basis; center c=2, three nonzero star edges.
    vertices=list(itertools.product(range(4),repeat=2)); center=2
    def qedge(i,j): return (i==center)!=(j==center)
    edges=[(u,v) for i,u in enumerate(vertices) for v in vertices[i+1:]
           if qedge(u[0],v[0]) and qedge(u[1],v[1])]
    adjacency={v:set() for v in vertices}
    rows=[]
    for u,v in edges:
        adjacency[u].add(v);adjacency[v].add(u)
        row=[0]*8
        row[u[0]]+=1;row[4+u[1]]+=1
        row[v[0]]-=1;row[4+v[1]]-=1
        rows.append(row)
    remaining=set(vertices);components=[]
    while remaining:
        front=[min(remaining)];comp=set(front)
        while front:
            x=front.pop()
            for y in adjacency[x]-comp: comp.add(y);front.append(y)
        remaining-=comp;components.append(sorted(comp))
    components.sort(key=len)
    assert [len(v) for v in components]==[6,10]
    assert len(edges)==18
    rank=row_rank(rows);assert rank==6
    null_generators=[[1]*4+[0]*4,[0]*4+[1]*4]
    assert all(sum(a*b for a,b in zip(row,v))==0
               for row in rows for v in null_generators)
    assert 8-rank==len(null_generators)
    uniform=[F(1,4)]*4
    test=[F(3,8),F(1,8),F(1,8),F(3,8)]
    def violates(r,s):
        return sum(r[u[0]]*s[u[1]] != r[v[0]]*s[v[1]] for u,v in edges)
    assert violates(uniform,uniform)==0 and violates(test,test)>0
    # Witness: a commuting phase changes a purification but not its marginal.
    overlap=1-2*test[0]
    joint_distance_squared=1-overlap*overlap
    assert overlap==F(1,4) and joint_distance_squared==F(15,16)
    # A separate equally spaced qutrit: declared alternative spectrum/permission.
    z=[1,-2,1]; e=[0,1,2]
    assert sum(z)==sum(a*b for a,b in zip(z,e))==0
    q=[F(1,2),F(1,3),F(1,6)]
    weight_difference=abs(q[0]*q[2]-q[1]*q[1])
    assert weight_difference==F(1,36)
    inverted=[F(1,7),F(2,7),F(4,7)]
    assert inverted[0]*inverted[2]==inverted[1]**2
    # Independent spectral-projector formula for the native |1> weight:
    # (P_lambda)11 = minor11(lambda I-h)/det'(lambda I-h).
    def imul(x,y):
        products=[a*b for a in x for b in y]
        return (min(products),max(products))
    x=tuple(F(v,DEN) for v in boxes[0])
    minor0=imul(x,(x[0]+F(1,40),x[1]+F(1,40)))
    minor=(minor0[0]-F(3,6400),minor0[1]-F(3,6400))
    x2=imul(x,x)
    deriv=(3*x2[0]-F(39,20)*x[1]-F(2095,32000),
           3*x2[1]-F(39,20)*x[0]-F(2095,32000))
    assert min(minor)>0 and min(deriv)>0
    q2=(minor[0]/deriv[1],minor[1]/deriv[0])
    coldpath=STAGE/'1202'/'results.json'
    coldjson=json.loads(coldpath.read_text(encoding='utf-8'))
    cold=[F(v) for v in coldjson['one_marginal_rational']]
    assert sum(cold)==1 and min(cold)>0
    pdiff=cold[0]**2-cold[2]**2
    coldcomm=(F(1,10)*q2[0]*pdiff,F(1,10)*q2[1]*pdiff)
    assert F('.0003202089318')<coldcomm[0]<coldcomm[1]<F('.0003202089319')
    crossrows=[]
    for left,right in [(0,1),(1,2)]:
        for level in [1,2]:
            row=[0]*9
            row[3*left+level]+=1;row[3*right]+=1
            row[3*left]-=1;row[3*right+level]-=1
            crossrows.append(row)
    conv=[1,-2,1]+[0]*6
    rank_exchange=row_rank(crossrows)
    rank_calibrated=row_rank(crossrows+[conv])
    rank_delete_edge=row_rank(crossrows[:2]+[conv])
    assert (rank_exchange,rank_calibrated,rank_delete_edge)==(4,5,3)
    files={
        '1201/proof.md':STAGE/'1201'/'proof.md',
        '1204/proof.md':STAGE/'1204'/'proof.md',
        '1202/results.json':coldpath,
        '_shared/reference_preparation_screen.md':STAGE/'_shared'/'reference_preparation_screen.md',
    }
    return {
        'status':'independent exact algebra and scoped finite diagnostics; not an implementation certificate',
        'all_checks_passed':True,
        'method':'Python standard library only; no import of any author program',
        'source_sha256':{name:digest(path) for name,path in files.items()},
        'independent_script_sha256':digest(Path(__file__)),
        'native_spectrum':{
            'cubic_coefficients_low_to_high':COEFF,
            'full_characteristic_coefficients_low_to_high':[str(v) for v in quartic],
            'mod11_values_x0_to10':residues,
            'irreducible_cubic_over_Q_from_mod11':True,
            'discriminant':disc,'discriminant_isqrt':sq,'discriminant_not_square':True,
            'galois_group':'S3 (irreducible cubic and nonsquare discriminant)',
            'f_at_exact_center_energy_1':f1,
            'root_denominator':str(DEN),'root_integer_boxes':[[str(x) for x in b] for b in boxes],
            'root_midpoints':[decimal(F(lo+hi,2*DEN)) for lo,hi in boxes],
        },
        'fixed_n_resonance':{
            'root_sum':'39/40',
            'two_equation_determinant':str(determinant),
            'integer_energy_relation_without_fixed_n':relation,
            'relation_particle_number_change':sum(relation),
            'fixed_n_integer_resonance_lattice':'{0}',
            'real_kernel_dimension':2,
            'invariance':'Every bare-energy-preserving unitary fixes the diagonal iid system marginal for every finite n.',
            'not_claimed':'No claim that arbitrary purifying-reference correlations remain fixed.',
            'purification_phase_witness_overlap':str(overlap),
            'purification_phase_witness_trace_distance_squared':str(joint_distance_squared),
        },
        'qq_full_product_stationarity':{
            'assumptions':['g != 0','r and s individually stationary for nondegenerate h','full rank','initial product r tensor s','unchanged for every t'],
            'q_star_nonzero_overlap_argument':'If outer eigenvector component v1=0, row1 forces v0=0 and row0 forces v3=0; contradiction.',
            'qq_undirected_edge_count':len(edges),
            'components':[[list(v) for v in comp] for comp in components],
            'log_product_constraint_rank':rank,
            'log_product_constraint_nullspace_basis':null_generators,
            'unique_normalized_full_rank_product':['I4/4','I4/4'],
            'test_rational_population':[str(v) for v in test],
            'test_self_product_violated_edges':violates(test,test),
            'cold_reference_commutator':{
                'preparation':'one material from each of two independent cold pools',
                'entry':'absolute <00|[H_AB,r_C tensor r_C]|cc>',
                'method':'spectral-projector cofactor divided by characteristic-polynomial derivative; independent rational interval roots',
                'strict_lower':'0.0003202089318',
                'strict_upper':'0.0003202089319',
                'interval_lower_exact':str(coldcomm[0]),
                'interval_upper_exact':str(coldcomm[1]),
                'finite_time_measurement_claimed':False,
            },
            'excluded_claims':['not a finite-time zero-flow criterion','not a claim about correlated full-H Gibbs states','not arbitrary auxiliary or battery cycles'],
        },
        'near_resonance_diagnostics':{
            'general_bound':'delta_n <= n*(e_max-e_min)/(binomial(n+3,3)-1), hence O(n^-2)',
            'nearest_configuration_swap_commutator_norm':'delta_n exactly on the two-state swap block',
            'scope':'No n-uniform positive gap. No lower bound on the iid-state change or physically available collision is implied.',
            'samples':[gaps(n,boxes) for n in [1,2,10,30]],
        },
        'permitted_conversion_criterion':{
            'full_rank_reference_condition':'For each permitted equal-energy type conversion z, z dot log(p) = 0.',
            'gibbs_selection_condition':'real span of permitted integer z equals real ker([1;E]); then log p = alpha - beta E',
            'inverse_temperature_sign':'beta may be any real value; noninversion is an additional condition',
            'three_subject_qutrit_exact_ranks':{
                'cross_subject_exchange_only':rank_exchange,
                'with_one_conversion':rank_calibrated,
                'with_conversion_but_delete_last_subject_edge':rank_delete_edge,
                'number_of_log_population_variables':9,
            },
            'qutrit_example_only':{
                'energies':e,'conversion':z,'strings':[[0,2],[1,1]],
                'nonstationary_example':[str(v) for v in q],
                'product_state_trace_distance_after_string_swap':str(weight_difference),
                'inverted_stationary_example':[str(v) for v in inverted],
                'not_native_h_or_native_QQ_permission':True,
            }
        }
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=main_result()
    dest=HERE/'independent_results.json'
    if args.write:
        if dest.exists(): raise SystemExit('refuse overwrite of existing independent_results.json')
        dest.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads(dest.read_text(encoding='utf-8'))
        if result != saved: raise SystemExit('independent recomputation mismatch')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__': main()
