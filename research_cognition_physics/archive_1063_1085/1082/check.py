"""1082 exact finite certificates: independent random representation.

The general theorem is analytic.  This checker verifies rational examples,
finite experimental margins, and the scope of lottery coarse-graining.
Default prints JSON only; --check recomputes and compares without writing.
"""
from fractions import Fraction as F
from itertools import combinations
from math import isqrt
from pathlib import Path
import argparse
import hashlib
import json

ZERO=(F(0),F(0),F(0))

def vec(*xs):
    return tuple(F(x) for x in xs)

def add(x,y):
    return tuple(a+b for a,b in zip(x,y))

def sub(x,y):
    return tuple(a-b for a,b in zip(x,y))

def mul(t,x):
    return tuple(F(t)*a for a in x)

def dot(x,y):
    return sum((a*b for a,b in zip(x,y)),F(0))

def det3(a):
    x,y,z=a
    return x[0]*(y[1]*z[2]-y[2]*z[1])-x[1]*(y[0]*z[2]-y[2]*z[0])+x[2]*(y[0]*z[1]-y[1]*z[0])

def sqrtrational(x):
    x=F(x)
    n,d=isqrt(x.numerator),isqrt(x.denominator)
    if n*n!=x.numerator or d*d!=x.denominator:
        raise ValueError('this finite example was expected to have a rational square root')
    return F(n,d)

def lottery(*xs):
    p=tuple(F(x) for x in xs)
    if min(p)<0 or sum(p,F(0))!=1:
        raise ValueError('normalized nonnegative lottery required')
    return p

def pure(n,i):
    return tuple(F(j==i) for j in range(n))

def directions(V, gain):
    return {(i,j):mul(gain(i,j),sub(V[i],V[j])) for i in range(len(V)) for j in range(len(V))}

def mixed_direction(p,q,d):
    value=ZERO
    for i,pi in enumerate(p):
        for j,qj in enumerate(q):
            value=add(value,mul(pi*qj,d[i,j]))
    return value

def moments(p,V,w):
    weight=sum((pi*wi for pi,wi in zip(p,w)),F(0))
    numerator=ZERO
    for pi,wi,vi in zip(p,w,V):
        numerator=add(numerator,mul(pi*wi,vi))
    return weight,numerator,mul(1/weight,numerator)

def outvec(x):
    return [str(a) for a in x]

class Audit:
    def __init__(self):
        self.require_calls=0
        self.groups=[]
    def require(self,condition,label):
        self.require_calls+=1
        if not condition:
            raise AssertionError(label)
    def group(self,name,data):
        self.groups.append({'name':name,'passed':True,'certificate':data})

def effect_routing(audit,p,q,d,scale,probe):
    # Four Pauli coefficients (I,X,Y,Z) exactly describe the effects.
    # Routing choices are classical, mutually independent, and independent
    # of the one unknown probe. Fine records (i,j) can remain accessible.
    coefficients=[F(0)]*4
    branch_probability=F(0)
    branch_total=F(0)
    support=0
    for i,pi in enumerate(p):
        for j,qj in enumerate(q):
            weight=pi*qj
            if not weight:
                continue
            support+=1
            E=(F(1,2),*mul(scale,d[i,j]))
            for k in range(4):
                coefficients[k]+=weight*E[k]
            prob=E[0]+dot(probe,E[1:])
            audit.require(0<=prob<=1,'selected branch probability is physical')
            branch_probability+=weight*prob
            branch_total+=weight
    vector=mixed_direction(p,q,d)
    expected=(F(1,2),*mul(scale,vector))
    audit.require(branch_total==1 and tuple(coefficients)==expected,
                  'independent classical routing gives the exact coarse effect')
    audit.require(branch_probability==F(1,2)+scale*dot(probe,vector),
                  'branch probabilities reproduce the same coarse readout')
    return {'left':outvec(p),'right':outvec(q),'nonzero_fine_branches':support,
            'fine_probability_sum':str(branch_total),
            'coarse_effect_Pauli_coefficients':outvec(coefficients),
            'probe':outvec(probe),'coarse_probability':str(branch_probability)}

def run():
    audit=Audit()
    ex,ey,ez=vec(1,0,0),vec(0,1,0),vec(0,0,1)
    V=[ZERO,ex,ey,ez,vec(2,0,0),ZERO]
    w=tuple(map(F,(1,2,3,4)))+(F(5,2),F(2))
    gain=lambda i,j:w[i]*w[j]
    d=directions(V,gain)
    products=(gain(0,1)*gain(2,3),gain(0,2)*gain(1,3),gain(0,3)*gain(1,2))
    audit.require(products==(F(24),)*3,'factorized four-anchor opposite products agree')
    w0=sqrtrational(gain(0,1)*gain(0,2)/gain(1,2))
    recovered=(w0,)+(tuple(gain(0,i)/w0 for i in (1,2,3)))
    audit.require(recovered==w[:4] and min(recovered)>0,'recover positive anchor weights')
    audit.require(det3([sub(V[i],V[0]) for i in (1,2,3)])==1,
                  'independent affine anchor directions')
    recovered_extra={}
    for i in (4,5):
        # Use a distinct-position anchor; the coincident edge i=5,j=0
        # carries no gain data and is deliberately not divided by.
        j=2
        delta=sub(V[i],V[j])
        k=next(k for k,x in enumerate(delta) if x)
        wi=d[i,j][k]/(w[j]*delta[k])
        audit.require(wi==w[i] and wi>0,'recover extra partner weight from a nonzero edge')
        recovered_extra[str(i)]=str(wi)
    audit.require(V[4]==mul(2,V[1]) and V[5]==V[0],
                  'positive example contains collinear and repeated positions')
    audit.require(mul(F(77,5),sub(V[0],V[5]))==d[0,5]==ZERO,
                  'arbitrary zero-edge coefficient is not observable gain data')
    for i,j in combinations(range(6),2):
        audit.require(d[i,j]==mul(w[i]*w[j],sub(V[i],V[j])),
                      'all extended edges use the recovered positive weights')
    audit.group('factorization_with_collinear_and_repeated_partners',{
        'zero_based_positions':[outvec(v) for v in V],
        'weights':outvec(w),'anchor_opposite_products':outvec(products),
        'recovered_anchor_weights':outvec(recovered),'recovered_extra_weights':recovered_extra,
        'edge_vectors':{f'{i},{j}':outvec(d[i,j]) for i,j in combinations(range(6),2)},
        'zero_edge':{'pair':[0,5],'arbitrary_literal_coefficient':'77/5',
                     'factorized_coefficient':str(w[0]*w[5]),'same_observable_vector':['0','0','0']}})

    P=lottery(F(1,2),F(1,2),0,0,0,0)
    Q=lottery(0,0,F(1,2),F(1,2),0,0)
    T=lottery(F(1,6),0,F(1,6),F(1,6),F(1,3),F(1,6))
    pool=[pure(6,i) for i in range(6)]+[P,Q,T]
    moment_checks=[]
    for p in pool:
        wp,bp,vp=moments(p,V,w)
        audit.require(wp>0,'every normalized finite lottery has positive total weight')
        for q in pool:
            wq,bq,vq=moments(q,V,w)
            direct=mixed_direction(p,q,d)
            bilinear=sub(mul(wq,bp),mul(wp,bq))
            normalized=mul(wp*wq,sub(vp,vq))
            audit.require(direct==bilinear==normalized,'exact weighted-moment comparison identity')
        moment_checks.append({'lottery':outvec(p),'W':str(wp),'B':outvec(bp),'normalized_V':outvec(vp)})
    mixing=F(1,3)
    mix=lottery(*((1-mixing)*p+mixing*q for p,q in zip(P,Q)))
    wp,bp,vp=moments(P,V,w)
    wq,bq,vq=moments(Q,V,w)
    wm,bm,vm=moments(mix,V,w)
    audit.require(wm==(1-mixing)*wp+mixing*wq and bm==add(mul(1-mixing,bp),mul(mixing,bq)),
                  'mixture updates unnormalized W and B linearly')
    geometric_fraction=mixing*wq/wm
    audit.require(vm==add(mul(1-geometric_fraction,vp),mul(geometric_fraction,vq))
                  and geometric_fraction!=mixing,'visible mixture fraction is weight adjusted')
    audit.group('exact_lottery_moments_and_composition',{
        'identity':'D(P,Q)=W_Q B_P-W_P B_Q=W_P W_Q(V_P-V_Q)',
        'normalization':'W_P=sum p_i w_i; B_P=sum p_i w_i v_i; V_P=B_P/W_P',
        'normalized_lotteries_checked':len(pool),'ordered_pair_identities_checked':len(pool)**2,
        'moments':moment_checks,
        'mixture':{'classical_Q_fraction':str(mixing),'probabilities':outvec(mix),
                   'W':str(wm),'B':outvec(bm),'normalized_V':outvec(vm),
                   'geometric_Q_fraction':str(geometric_fraction)},
        'scope':'This is a lottery/coarse-record composition, not an autonomous physical meeting of endpoints.'})

    positive_scale=F(1,128)
    p0=F(1,2)+positive_scale*dot(ex,d[0,1])
    p5=F(1,2)+positive_scale*dot(ex,d[5,1])
    audit.require(V[0]==V[5] and w[0]!=w[5] and p0-p5==F(1,64),
                  'same visible coordinate with different weight has different full probabilities')
    audit.group('weight_is_visible_to_the_original_probability_task',{
        'same_position_partners':[0,5],'fixed_other_partner':1,'probe':outvec(ex),
        'effect_scale':str(positive_scale),'probabilities':[str(p0),str(p5)],
        'finite_gap':str(p0-p5),
        'conclusion':'Only sign-order identity can discard this weight; full probabilities cannot.'})

    # Nonfactorized gains: indices 1,2 mean the e1/e2 edge, not 0/e1.
    Vn=V[:4]
    gn=lambda i,j:F(2) if {i,j}=={1,2} else F(1)
    dn=directions(Vn,gn)
    opposite=(gn(0,1)*gn(2,3),gn(0,2)*gn(1,3),gn(0,3)*gn(1,2))
    audit.require(opposite==(F(1),F(1),F(2)),'nonfactorized opposite products differ')
    r=vec(F(1,4),F(1,2),F(3,4))
    audit.require(dot(r,r)==F(7,8)<=1,'cycle witness is a physical Bloch preparation')
    cyc=[lottery(0,F(1,6),F(2,3),F(1,6)),
         lottery(0,F(1,3),F(1,6),F(1,2)),
         lottery(F(1,6),F(1,6),0,F(2,3))]
    cycle_vectors=[mixed_direction(cyc[i],cyc[(i+1)%3],dn) for i in range(3)]
    margins=[dot(r,x) for x in cycle_vectors]
    audit.require(margins==[F(1,144),F(1,144),F(1,72)] and min(margins)>0,
                  'three finite lotteries have a strict cycle')
    small_scale=F(1,8)
    probabilities=[F(1,2)+small_scale*x for x in margins]
    audit.group('pure_weak_order_but_strict_lottery_cycle',{
        'positions':[outvec(v) for v in Vn],
        'exceptional_gain_pair_zero_based':[1,2],'exceptional_gain':'2',
        'other_gains':'1','opposite_products':outvec(opposite),
        'pure_partner_order':'sign r.d_ij = sign(r.v_i-r.v_j), since all gains are positive',
        'probe':outvec(r),'probe_norm_squared':'7/8',
        'lotteries_P_Q_R':[outvec(p) for p in cyc],
        'D_PQ_D_QR_D_RP':[outvec(v) for v in cycle_vectors],
        'unscaled_cyclic_margins':outvec(margins),
        'effect_scale':str(small_scale),'coarse_probabilities':outvec(probabilities),
        'probability_margins':outvec([small_scale*x for x in margins]),
        'signs':[1,1,1]})

    tieP=lottery(0,F(1,3),0,F(2,3))
    tieQ=pure(4,2)
    tieR=pure(4,0)
    tieprobe=vec(F(1,2),0,F(-1,2))
    tie_vectors=[mixed_direction(tieP,tieQ,dn),mixed_direction(tieQ,tieR,dn),mixed_direction(tieP,tieR,dn)]
    tie_margins=[dot(tieprobe,x) for x in tie_vectors]
    audit.require(dot(tieprobe,tieprobe)==F(1,2),'tie witness is a physical preparation')
    audit.require(tie_margins==[F(0),F(0),F(-1,6)],'binary random representative breaks tie substitution')
    audit.require(small_scale*tie_margins[2]==F(-1,48),'binary tie failure has a finite actual margin')
    audit.group('binary_representative_tie_substitution_failure',{
        'probe':outvec(tieprobe),'P':outvec(tieP),'Q':outvec(tieQ),'R':outvec(tieR),
        'pair_order':['P,Q','Q,R','P,R'],'unscaled_margins':outvec(tie_margins),
        'effect_scale':str(small_scale),
        'coarse_probabilities':outvec([F(1,2)+small_scale*x for x in tie_margins]),
        'nonzero_probability_difference':'-1/48'})

    legal=[]
    for name,n,dd,s in [('factorized',6,d,positive_scale),('nonfactorized',4,dn,small_scale)]:
        maximum=max(s*s*dot(dd[i,j],dd[i,j]) for i in range(n) for j in range(n))
        audit.require(maximum<F(1,4),'all branch effects are strictly legal')
        audit.require(all(dd[j,i]==mul(-1,dd[i,j]) for i in range(n) for j in range(n)),
                      'reverse branches have complementary effects')
        legal.append({'menu':name,'effect_scale':str(s),'maximum_Bloch_norm_squared':str(maximum),
                      'positivity_threshold':'1/4','strict_margin':str(F(1,4)-maximum)})
    routing=[effect_routing(audit,cyc[i],cyc[(i+1)%3],dn,small_scale,r) for i in range(3)]
    routing.extend([effect_routing(audit,tieP,tieQ,dn,small_scale,tieprobe),
                    effect_routing(audit,tieQ,tieR,dn,small_scale,tieprobe),
                    effect_routing(audit,tieP,tieR,dn,small_scale,tieprobe),
                    effect_routing(audit,P,Q,d,positive_scale,ex)])
    audit.group('legal_qubit_effects_and_actual_classical_routing',{
        'convention':'rho_r=(I+r.sigma)/2; E_ij=I/2+scale*d_ij.sigma',
        'menus':legal,'selected_mixture_checks':routing,
        'protocol':{'selection':'Independent finite classical I~P and J~Q, independent of the unknown probe.',
                    'probe_uses':1,'unknown_probe_copying':False,
                    'branch':'Apply the already specified binary instrument for the selected pair (I,J).',
                    'coarse_effect':'sum p_i q_j E_ij',
                    'fine_records':'I,J can remain recorded in the complete system; only the comparison report is coarse.',
                    'reference_scope':'The effect identity is operator-level and hence holds with any passive reference; no equality of branch poststates is claimed.',
                    'source_scope':'Finite routing and fixed instruments are declared permissions, not derived autonomous-H controls.'}})

    # Full lottery compatibility constrains its comparison menu; it does not
    # specify an independently supplied physical endpoint-meeting operation.
    # Here all comparison weights are 1, while c is a distinct resource label.
    state0=(ZERO,F(1))
    state0prime=(ZERO,F(2))
    state1=(ex,F(1))
    def other_meeting(left,right):
        a,c=left
        b,e=right
        return mul(1/(c+e),add(mul(c,a),mul(e,b))), (c+e)/2
    vm0,c0=other_meeting(state0,state1)
    vm1,c1=other_meeting(state0prime,state1)
    audit.require(vm0==mul(F(1,2),ex) and vm1==mul(F(1,3),ex),
                  'separate actual meeting can depend on a hidden resource label')
    pm0=F(1,2)+small_scale*dot(ex,vm0)
    pm1=F(1,2)+small_scale*dot(ex,vm1)
    audit.require(pm0==F(9,16) and pm1==F(13,24) and pm0-pm1==F(1,48),
                  'separate meeting yields a finite final probability difference')
    audit.group('lottery_compatibility_does_not_fix_a_separate_meeting',{
        'comparison_contract':'d_xy=v_x-v_y, all endpoint weights w=1; every independent lottery obeys ordinary expected-coordinate order.',
        'additional_label':'c in [1,2], distinct from the comparison weight w',
        'actual_meeting':'m((v,c),(u,e))=((c*v+e*u)/(c+e),(c+e)/2)',
        'inputs':{'left_visible':['0','0','0'],'left_resources':['1','2'],
                  'right_visible':['1','0','0'],'right_resource':'1'},
        'meeting_positions':[outvec(vm0),outvec(vm1)],
        'fixed_final_comparison':'against visible origin, probe +X, effect scale 1/8',
        'probabilities':[str(pm0),str(pm1)],'finite_gap':'1/48',
        'scope':'Abstract countermodel of an inference, not a new apparatus or a demand that position predict every future.'})

    return {'round':1082,'status':'exact_fraction_certificates_passed',
            'scientific_calibration_groups':1,'validation_group_count':len(audit.groups),
            'require_call_count':audit.require_calls,
            'arithmetic':'fractions.Fraction; no floating-point scientific computation',
            'groups':audit.groups,
            'scope':{'general_factorization_theorem_proved_by_finite_examples':False,
                     'binary_random_representative_tie_substitution_is_an_additional_cognitive_hypothesis':True,
                     'full_lottery_weak_order_is_a_conditional_consequence':True,
                     'unknown_quantum_input_copied':False,
                     'random_lottery_is_actual_meeting_endpoint':False,
                     'weight_identified_as_mass_energy_or_time':False,
                     'autonomous_physical_implementation_derived':False,
                     'actual_spatial_dimension_derived':False},
            'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true',help='read-only exact comparison with saved results.json')
    args=parser.parse_args()
    result=run()
    if args.check:
        path=Path(__file__).with_name('results.json')
        if result!=json.loads(path.read_text(encoding='utf-8')):
            raise AssertionError('saved results differ from exact recomputation')
        print(json.dumps({'round':1082,'status':'exact_saved_results_match',
                          'scientific_calibration_groups':1,
                          'validation_group_count':result['validation_group_count'],
                          'require_call_count':result['require_call_count'],
                          'checker_sha256':result['checker_sha256'],
                          'results_sha256':hashlib.sha256(path.read_bytes()).hexdigest()},
                         ensure_ascii=False,indent=2))
    else:
        print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
