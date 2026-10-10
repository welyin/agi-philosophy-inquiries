"""Round 1081: exact finite certificates for comparison compatibility.

Standard library only; all scientific arithmetic uses Fraction.  A positive
triangle circuit certifies every preparation, by taking its scalar product
with an arbitrary Bloch vector.  Examples do not prove the general theorem.
Default: print deterministic JSON, write nothing.  --check: compare the saved
results.json to a fresh computation, also without writing any file.
"""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import argparse
import hashlib
import json

ZERO = (F(0), F(0), F(0))

def vec(*xs):
    return tuple(F(x) for x in xs)

def add(x, y):
    return tuple(a+b for a,b in zip(x,y))

def sub(x, y):
    return tuple(a-b for a,b in zip(x,y))

def scale(a, x):
    return tuple(F(a)*b for b in x)

def dot(x, y):
    return sum((a*b for a,b in zip(x,y)), F(0))

def rank(rows):
    a = [list(map(F,row)) for row in rows]
    if not a:
        return 0
    pivot = 0
    for col in range(len(a[0])):
        candidate = next((i for i in range(pivot,len(a)) if a[i][col]), None)
        if candidate is None:
            continue
        a[pivot],a[candidate] = a[candidate],a[pivot]
        divisor = a[pivot][col]
        a[pivot] = [x/divisor for x in a[pivot]]
        for i in range(len(a)):
            if i != pivot:
                m = a[i][col]
                a[i] = [x-m*y for x,y in zip(a[i],a[pivot])]
        pivot += 1
        if pivot == len(a):
            break
    return pivot

def determinant(rows):
    a = [list(map(F,row)) for row in rows]
    n = len(a)
    if any(len(row) != n for row in a):
        raise ValueError('square matrix required')
    out = F(1)
    for col in range(n):
        candidate = next((i for i in range(col,n) if a[i][col]),None)
        if candidate is None:
            return F(0)
        if candidate != col:
            a[col],a[candidate] = a[candidate],a[col]
            out = -out
        pivot = a[col][col]
        out *= pivot
        for i in range(col+1,n):
            m = a[i][col]/pivot
            for j in range(col+1,n):
                a[i][j] -= m*a[col][j]
    return out

def table(n, upper):
    d = {(i,i):ZERO for i in range(n)}
    for i,j in combinations(range(n),2):
        d[i,j] = upper[i,j]
        d[j,i] = scale(-1,upper[i,j])
    return d

def set_edge(d, i, j, value):
    d[i,j] = value
    d[j,i] = scale(-1,value)

def serial_table(n, d):
    return {f'{i+1},{j+1}':[str(x) for x in d[i,j]]
            for i,j in combinations(range(n),2)}

class Audit:
    def __init__(self):
        self.calls = 0
        self.groups = []
    def require(self, condition, label):
        self.calls += 1
        if not condition:
            raise AssertionError(label)
    def group(self, name, result):
        self.groups.append({'name':name,'passed':True,'certificate':result})

def circuit(audit, d, tri, weights, positive=True):
    audit.require(all(w > 0 if positive else w >= 0 for w in weights)
                  and any(w > 0 for w in weights), 'circuit coefficient signs')
    total = ZERO
    for k,w in enumerate(weights):
        total = add(total,scale(w,d[tri[k],tri[(k+1)%3]]))
    audit.require(total == ZERO, 'exact triangle circuit closure')
    return {'cycle':[i+1 for i in tri], 'weights':[str(F(w)) for w in weights],
            'strictly_positive':positive}

def positive_factor(audit, direction, delta):
    if delta == ZERO:
        audit.require(direction == ZERO, 'coincident positions require a zero edge')
        return F(1)  # Arbitrary positive multiplier of the zero vector.
    k = next(i for i,x in enumerate(delta) if x)
    factor = direction[k]/delta[k]
    audit.require(factor > 0 and scale(factor,delta) == direction,
                  'positive, not merely signed, difference representation')
    return factor

def radial_from_one_edge(audit, ux, vj, dxj):
    # t ux - h dxj = vj, h=1/a_xj.  Only two independent coordinates
    # are used to infer t; the third is then independently checked.
    pair = next(((k,l) for k,l in combinations(range(3),2)
                 if ux[k]*(-dxj[l])-ux[l]*(-dxj[k])),None)
    audit.require(pair is not None, 'independent radial inference equations')
    k,l = pair
    den = ux[k]*(-dxj[l])-ux[l]*(-dxj[k])
    t = (vj[k]*(-dxj[l])-vj[l]*(-dxj[k]))/den
    h = (ux[k]*vj[l]-ux[l]*vj[k])/den
    audit.require(t>0 and h>0 and sub(scale(t,ux),scale(h,dxj)) == vj,
                  'positive inferred radial coordinate')
    return t

def has_strict_cycle(n, d, probe):
    edges = {i:[j for j in range(n) if j!=i and dot(probe,d[i,j])>0]
             for i in range(n)}
    state = [0]*n
    def visit(i):
        state[i]=1
        for j in edges[i]:
            if state[j]==1 or (state[j]==0 and visit(j)):
                return True
        state[i]=2
        return False
    return any(state[i]==0 and visit(i) for i in range(n))

def run():
    audit = Audit()
    ex,ey,ez = vec(1,0,0),vec(0,1,0),vec(0,0,1)

    # Rank 2: every direction gives a full weak order, but no common
    # positively scaled difference representation exists.
    d2 = table(4,{(0,1):vec(-1,0,0),(0,2):vec(0,-1,0),
                  (0,3):vec(-2,-3,0),(1,2):vec(1,-1,0),
                  (1,3):vec(-1,-3,0),(2,3):vec(-2,F(-5,2),0)})
    cert2 = [circuit(audit,d2,tri,w) for tri,w in
             [((0,1,2),(1,1,1)),((0,1,3),(1,1,1)),
              ((0,2,3),(1,2,2)),((1,2,3),(7,8,9))]]
    audit.require(rank(list(d2.values())) == 2,'rank-two comparison span')
    parallel = []
    for i,j in combinations(range(4),2):
        dx,dy,_ = d2[i,j]
        row = [F(0)]*6
        for k,sign in ((i,1),(j,-1)):
            if k:
                row[2*(k-1)] += sign*(-dy)
                row[2*(k-1)+1] += sign*dx
        parallel.append(row)
    det2 = determinant(parallel)
    audit.require(det2 != 0 and rank(parallel)==6,'no nonconstant planar coordinates')
    # In R3 every edge direction has z=0, so all positions have a common
    # z coordinate.  The planar homogeneous system is therefore exhaustive.
    audit.group('rank2_full_weak_order_without_common_positions',{
        'one_based_edge_table':serial_table(4,d2),
        'all_triangle_positive_circuits':cert2,
        'parallelism_matrix':[[str(x) for x in row] for row in parallel],
        'parallelism_determinant':str(det2), 'parallelism_rank':6,
        'scope':'Positive circuits cover every Bloch preparation, not sampled preparations.'})

    # Rank 3 and no strict cycle are insufficient when ties are inconsistent.
    # A={1,3}, B={2,4}; all A-to-B vectors are ex.
    d3 = table(4,{(0,1):ex,(0,2):ey,(0,3):ex,
                  (1,2):scale(-1,ex),(1,3):ez,(2,3):ex})
    cert3 = [circuit(audit,d3,tri,w,False) for tri,w in
             [((0,1,2),(1,1,0)),((0,1,3),(1,0,1)),
              ((0,2,3),(0,1,1)),((1,2,3),(1,1,0))]]
    audit.require(rank(list(d3.values()))==3,'strict example comparison rank three')
    # For this axis-only example these 27 probes exhaust all sign patterns,
    # including the tie strata.  No nonphysical Bloch length is required:
    # positive rescaling into the unit ball preserves every sign.
    patterns = list(product((-1,0,1),repeat=3))
    audit.require(all(not has_strict_cycle(4,d3,vec(*r)) for r in patterns),
                  'all 27 realizable axis-sign patterns have no strict cycle')
    audit.require(dot(ey,d3[0,1])==dot(ey,d3[1,2])==0
                  and dot(ey,d3[0,2])>0,'tie nontransitivity at a physical pure preparation')
    audit.require(d3[0,1]==ex and d3[2,1]==ex and d3[0,2]==ey,
                  'two collinear comparisons force an incompatible third axis')
    audit.group('rank3_strict_acyclic_but_ties_fail',{
        'one_based_edge_table':serial_table(4,d3),
        'all_triangle_nonnegative_circuits':cert3,
        'exhaustive_axis_sign_patterns':len(patterns),
        'tie_failure':{'probe':['0','1','0'],'tied_pairs':[[1,2],[2,3]],
                       'strict_pair':[1,3]},
        'has_common_positive_difference_coordinates':False})

    # Positive input consists of comparisons, not pre-installed positions.
    d = table(4,{(0,1):scale(-1,ex),(0,2):scale(-1,ey),
                 (0,3):scale(-1,ez),(1,2):vec(2,-3,0),
                 (1,3):vec(10,0,-21),(2,3):vec(0,5,-7)})
    cert = [circuit(audit,d,tri,w) for tri,w in
            [((0,1,2),(2,1,3)),((0,1,3),(10,1,21)),
             ((0,2,3),(5,1,7)),((1,2,3),(5,3,1))]]
    u = {i:d[i,0] for i in (1,2,3)}
    audit.require(determinant([list(u[i]) for i in (1,2,3)])==1,
                  'three independent anchor rays')
    alpha12,beta12,alpha23,beta23,alpha31,beta31 = map(F,(2,3,5,7,21,10))
    ratio_product = (alpha12/beta12)*(alpha23/beta23)*(alpha31/beta31)
    audit.require(ratio_product==1,'face consistency fixes radial ratio product')
    r = {1:F(1),2:beta12/alpha12,3:(beta12/alpha12)*(beta23/alpha23)}
    v = {0:ZERO, **{i:scale(r[i],u[i]) for i in u}}
    factors = {(i,j):positive_factor(audit,d[i,j],sub(v[i],v[j]))
               for i,j in combinations(range(4),2)}
    audit.group('rank3_independent_anchor_ratio_closure',{
        'one_based_edge_table':serial_table(4,d),
        'all_triangle_positive_circuits':cert,
        'cyclic_ratio_product':str(ratio_product),
        'radial_scales_from_anchor_1':{str(i+1):str(x) for i,x in r.items()},
        'inferred_positions':{str(i+1):[str(x) for x in p] for i,p in v.items()},
        'positive_edge_factors':{f'{i+1},{j+1}':str(x) for (i,j),x in factors.items()},
        'scope':'Conditional exact lemma example; no autonomous physical realization is assumed.'})

    # A parallel anchor ray is not a coincident point: infer its different
    # radius from two independent coordinate equations.
    de = dict(d)
    for j,a in [(0,ex),(1,ex),(2,vec(4,-3,0)),(3,vec(20,0,-21))]:
        set_edge(de,4,j,a)
    de[4,4] = ZERO
    t4 = radial_from_one_edge(audit,de[4,0],v[2],de[4,2])
    ve = {**v,4:scale(t4,de[4,0])}
    audit.require(t4==2 and ve[4]!=ve[1],'parallel ray distinct radial position')
    # Same inferred position, different amplitudes: comparisons from the
    # added partner are twice those of partner 2, except the zero edge.
    for j in range(5):
        set_edge(de,5,j,ZERO if j==1 else scale(2,de[1,j]))
    de[5,5] = ZERO
    t5 = radial_from_one_edge(audit,de[5,0],ve[2],de[5,2])
    ve[5] = scale(t5,de[5,0])
    audit.require(t5==F(1,2) and ve[5]==ve[1] and de[5,1]==ZERO,
                  'nonzero-star coincident partners')
    # Zero star means tied to the origin in every preparation.  Other edges
    # need only be positive multiples, not equal vectors.
    for j in range(6):
        set_edge(de,6,j,scale(3,de[0,j]))
    de[6,6] = ZERO
    ve[6] = ve[0]
    audit.require(de[6,0]==ZERO and de[6,2]!=de[0,2],
                  'zero-star ties do not require equal comparison amplitudes')
    extended_factors = {(i,j):positive_factor(audit,de[i,j],sub(ve[i],ve[j]))
                        for i,j in combinations(range(7),2)}
    extended_cert = []
    for tri in combinations(range(7),3):
        weights=[]
        for k in range(3):
            i,j=tri[k],tri[(k+1)%3]
            weights.append(1/extended_factors[tuple(sorted((i,j)))])
        extended_cert.append(circuit(audit,de,tri,weights))
    audit.require(len(extended_cert)==35,'all extended triples checked')
    audit.group('parallel_rays_and_zero_comparison_classes',{
        'one_based_edge_table':serial_table(7,de),
        'inferred_positions':{str(i+1):[str(x) for x in p] for i,p in ve.items()},
        'parallel_partner_5_radial_scale':str(t4),
        'coincident_partner_6_radial_scale':str(t5),
        'zero_edges':[[i+1,j+1] for i,j in combinations(range(7),2) if de[i,j]==ZERO],
        'positive_edge_factors':{f'{i+1},{j+1}':str(x) for (i,j),x in extended_factors.items()},
        'all_triangle_positive_circuits':extended_cert,
        'zero_edge_multiplier_convention':'1; any strictly positive value represents 0 = a*0.'})

    # A finite tetrahedron and a one-parameter moment curve both affinely
    # span R3.  Their actual topological dimension is not inferred by rank.
    audit.require(rank([sub(v[i],v[0]) for i in (1,2,3)])==3,
                  'finite tetrahedron has affine span three')
    parameters = [F(0),F(1,3),F(2,3),F(1)]
    curve = [vec(t,t*t,t*t*t) for t in parameters]
    vandermonde = determinant([list(sub(p,curve[0])) for p in curve[1:]])
    audit.require(vandermonde!=0 and all(p[0]==t for p,t in zip(curve,parameters)),
                  'nonplanar curve sample and inverse first-coordinate projection')
    audit.group('affine_span_is_not_actual_spatial_dimension',{
        'finite_tetrahedron':{'point_count':4,'affine_span':3,
                             'topological_dimension_by_finite_set_theorem':0},
        'moment_curve':{'definition':'t -> (t,t^2,t^3), 0 <= t <= 1',
                        'sample_parameters':[str(t) for t in parameters],
                        'sample_affine_determinant':str(vandermonde),
                        'inverse':'projection onto the first coordinate',
                        'topological_dimension_by_homeomorphism_to_interval':1},
        'scope':'Rank and rational identities are computed; the dimension claims use the stated elementary theorems.'})

    # E_ij=I/2+s d_ij.sigma.  For a qubit, E is an effect exactly when
    # |s d|^2 <= 1/4; reversal is I-E and the neutral state is balanced.
    legal=[]
    for name,n,dd,s in [('rank2',4,d2,F(1,8)),('strict_rank3',4,d3,F(1,8)),
                        ('tetrahedral',4,d,F(1,128)),('degenerate_extension',7,de,F(1,128))]:
        squares = [s*s*dot(dd[i,j],dd[i,j]) for i,j in combinations(range(n),2)]
        audit.require(max(squares)<=F(1,4),'every finite comparison is a legal qubit effect')
        audit.require(all(dd[j,i]==scale(-1,dd[i,j]) for i,j in combinations(range(n),2)),
                      'reversed effects are complementary')
        audit.require(s>0,'common scaling preserves every preference sign')
        legal.append({'example':name,'common_positive_scale':str(s),
                      'maximum_effect_Bloch_norm_squared':str(max(squares)),
                      'allowed_upper':'1/4',
                      'strict_margin':str(F(1,4)-max(squares))})
    audit.group('legal_fixed_qubit_effects',{
        'convention':'E_ij = I/2 + s*d_ij.sigma; rho_r=(I+r.sigma)/2',
        'probability':'Tr(rho_r E_ij) = 1/2 + s*r.d_ij',
        'reversal':'E_ji = I-E_ij',
        'examples':legal,
        'instrument_existence':'Any such effect admits a binary quantum instrument; no old autonomous H implementation is claimed.'})
    return {'round':1081,'status':'exact_fraction_certificates_passed',
            'scientific_calibration_groups':1,
            'validation_group_count':len(audit.groups),
            'require_call_count':audit.calls,
            'arithmetic':'fractions.Fraction; no floating-point scientific arithmetic',
            'groups':audit.groups,
            'scope':{'general_representation_theorem_proved_by_examples':False,
                     'universal_preparation_scope_of_triangle_certificates':True,
                     'new_cognitive_contracts_automatically_derived':False,
                     'autonomous_physical_comparison_implementation_derived':False,
                     'actual_spatial_dimension_derived':False,
                     'numerical_tolerance_used':False},
            'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true',help='read-only exact comparison with saved results.json')
    args=parser.parse_args()
    result=run()
    if args.check:
        path=Path(__file__).with_name('results.json')
        stored=json.loads(path.read_text(encoding='utf-8'))
        if result!=stored:
            raise AssertionError('saved results.json differs from exact recomputation')
        print(json.dumps({'status':'exact_saved_results_match','round':1081,
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
